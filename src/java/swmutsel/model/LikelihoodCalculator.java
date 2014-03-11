package swmutsel.model;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.Constants;
import swmutsel.MatrixArrayPool;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.Pair;
import swmutsel.utils.PhyloUtils;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Felsenstein's pruning algorithm to calculate the likelihood for codon based models. Can deal with heterogenous models.
 *
 * @author Asif Tamuri (tamuri@ebi.ac.uk)
 */
public class LikelihoodCalculator {

    private final Tree tree;
    private final Map<String, Byte> states;
    private final Map<String, SubstitutionModel> cladeModels;
    private final Map<Node, String> nodeLabels = Maps.newHashMap();

    private final double[] tipConditional = new double[GeneticCode.CODON_STATES];
    private final double[] gapConditional = CoreUtils.repd(1.0, GeneticCode.CODON_STATES);

    private final String ROOT_MODEL_NAME;

    // TODO: this is can be a constant for the entire run...
    private final int[] senseCodons = GeneticCode.getInstance().getSenseCodons();

    private double[] probMatrix;

    private double logScaling = 0.0;

    // TODO: Improve memory usage of these conditionals...do we need to store everything??!
    // TODO: For example, strictly speaking, we only need as many arrays as we have concurrent threads running, and each can be reused
    // Could we have a reusable pool of 64*64 arrays? We would then only create 16 (or whatever), instead of 3598!
    //  private final double[][] tipConditionals;
    // private final double[][] conditionals;
    private double[][] conditionals;

    private final List<Partial> conditionalAtRoot = Lists.newArrayList();

    public LikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, LinkedHashMap<String, SubstitutionModel> models) {
        this.tree = tree;
        this.states = siteStates;

        this.cladeModels = models;

        // For non-homogeneous models, the first added entry is the "root" model (so you should use LinkedHashMap).
        this.ROOT_MODEL_NAME = cladeModels.keySet().iterator().next();

        // set up name lookup TODO: this is shared across all sites!
        for (Node n : PhyloUtils.externalNodes(tree)) {
            nodeLabels.put(n, n.getIdentifier().getName());
        }
        for (Node n : PhyloUtils.internalNodes(tree)) {
            nodeLabels.put(n, n.getIdentifier().getName());
            if (n.getIdentifier().getName().length() == 0) {
                nodeLabels.put(n, this.ROOT_MODEL_NAME);
            }
        }
    }

    @SuppressWarnings("unchecked")
    public LikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, SubstitutionModel model) {
        this(tree, siteStates, CoreUtils.getLinkedHashMap(Pair.of("ALL", model)));
    }

    public double getLogLikelihood() {
        boolean manageStorage = false;

        if (probMatrix == null) {
            getStorage();
            manageStorage = true;
        }

        conditionalAtRoot.clear();

        logScaling = 0.0;

        double[] conditionals = downTree();
        double[] f = cladeModels.get(ROOT_MODEL_NAME).getCodonFrequencies();

        double sum = 0.0;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) sum += conditionals[i] * f[i];

        if (sum < 0) sum = 0;

        if (manageStorage) {
            releaseStorage();
        }

        return Math.log(sum) + logScaling;
    }

    private String getNodeLabel(Node n) {
        return nodeLabels.get(n);
    }

    public void getStorage() {
        probMatrix = MatrixArrayPool.popCodonMatrix();
        conditionals = MatrixArrayPool.popConditionals();
    }

    public void releaseStorage() {
        // So we can keep the reference to the LikelihoodCalculator without keeping this baggage
        MatrixArrayPool.pushCodonMatrix(probMatrix);
        MatrixArrayPool.pushConditionals(conditionals);
        probMatrix = null;
        conditionals = null;
    }

    private double[] downTree() {
        //long start = CodeTimer.start();
        for (Node node : PhyloUtils.internalNodes(tree)) {

            double[] conditional = new double[GeneticCode.CODON_STATES];
            for (int j : senseCodons) conditional[j] = 1.0;

            for (int j = 0; j < node.getChildCount(); j++) {
                Node child = node.getChild(j);

                double[] lowerConditional;

                if (child.isLeaf()) {
                    if (GeneticCode.getInstance().isUnknownCodonState(states.get(getNodeLabel(child)))) {
                        lowerConditional = gapConditional;
                    } else {
                        Arrays.fill(tipConditional, 0.0);
                        tipConditional[states.get(getNodeLabel(child))] = 1.0;
                        lowerConditional = tipConditional;
                    }

                } else {
                    lowerConditional = conditionals[child.getNumber()];
                }

                if (child.getParent().isRoot()) {
                    conditionalAtRoot.add(new Partial(child, PhyloUtils.getNodeNumber(tree, child), Arrays.copyOf(lowerConditional, GeneticCode.CODON_STATES), child.getBranchLength()));
                }

                if (cladeModels.size() == 1) { // homogeneous model

                    // TODO: A probability matrix (of sorts) already exists in TdGCodonModel...do we need to create it again?
                    cladeModels.get(ROOT_MODEL_NAME).getTransitionProbabilities(probMatrix, child.getBranchLength());
                    updateIntraCladeConditionals(lowerConditional, conditional, probMatrix);

                } else { // non-homogeneous model

                    if (getNodeLabel(node).length() == 0 // the root of the tree is a parent without a label
                            || getNodeLabel(child).substring(0, 2).equals(getNodeLabel(node).substring(0, 2))) { // or we're not switching to a different model

                        // System.out.printf("getnodelabel = %s\n", getNodeLabel(child));



                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getTransitionProbabilities(probMatrix, child.getBranchLength());
                        updateIntraCladeConditionals(lowerConditional, conditional, probMatrix);

                    } else { // this is a hostshift!
                        final double[] probMatrix1 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
                        cladeModels.get(getNodeLabel(node).substring(0, 2)).getTransitionProbabilities(probMatrix, child.getBranchLength() * Constants.CLADE_BRANCH_SPLIT);
                        cladeModels.get(getNodeLabel(child).substring(0, 2)).getTransitionProbabilities(probMatrix1, child.getBranchLength() * (1 - Constants.CLADE_BRANCH_SPLIT));
                        updateInterCladeConditionals(lowerConditional, conditional, probMatrix, probMatrix1);
                    }
                }
            }

            if (!(node.getParent() == null) && !node.getParent().isRoot()) {
                scaleConditionals(node, conditional);
            }

            conditionals[node.getNumber()] = conditional;

        }
        //CodeTimer.store("downTree", start);
        return conditionals[tree.getRoot().getNumber()];
    }

    public double getLogLikelihoodForBranch(int node, double branchLength) {
        // TODO: store the true branchlengths here, so they can be updated with the distribtuedrunner
        probMatrix = MatrixArrayPool.popCodonMatrix();

        double[] sumPartial = new double[GeneticCode.CODON_STATES];
        Arrays.fill(sumPartial, 1.0);

        for (Partial p : conditionalAtRoot) {
            if (p.number == node) {
                cladeModels.get(ROOT_MODEL_NAME).getTransitionProbabilities(probMatrix, branchLength);
            } else {
                cladeModels.get(ROOT_MODEL_NAME).getTransitionProbabilities(probMatrix, p.branchLength);
            }
            updateIntraCladeConditionals(p.partial, sumPartial, probMatrix);
        }

        double lnL = 0;
        final double[] frequencies = cladeModels.get("ALL").getCodonFrequencies();
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) lnL += sumPartial[i] * frequencies[i];

        if (lnL < 0) lnL = 0;

        lnL = Math.log(lnL) + logScaling;

        MatrixArrayPool.pushCodonMatrix(probMatrix);

        return lnL;
    }

    public void setBranch(int node, double bl) {
        for (Partial p : conditionalAtRoot) {
            if (p.number == node) {
                p.branchLength = bl;
            }
        }
    }

    private void scaleConditionals(Node node, double[] conditionals) {
        if (node.getNumber() % Constants.SCALING_NODE_STEP == 0) {
            double scalingFactor = 0;
            for (double conditional : conditionals) {
                if (conditional > 0 && conditional > scalingFactor) {
                    scalingFactor = conditional;
                }
            }

            if (scalingFactor < Constants.SCALING_THRESHOLD) {
                for (int i = 0; i < conditionals.length; i++) {
                    conditionals[i] = conditionals[i] / scalingFactor;
                }
                logScaling += Math.log(scalingFactor);
            }
        }
    }

    private void updateInterCladeConditionals(double[] lowerConditional, double[] conditionals, double[] probMatrix0, double[] probMatrix1) {
        double[] probMatrix = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
        for (int i : senseCodons) {
            double branchProb = 0.0;
            for (int j : senseCodons) {
                for (int k : senseCodons) {
                    probMatrix[i * GeneticCode.CODON_STATES + j] += probMatrix0[i * GeneticCode.CODON_STATES + k] * probMatrix1[k * GeneticCode.CODON_STATES + j];
                }
                branchProb += lowerConditional[j] * probMatrix[i * GeneticCode.CODON_STATES + j];
            }
            conditionals[i] *= branchProb;
        }
    }

    private void updateIntraCladeConditionals(double[] lowerConditional, double[] conditionals, double[] probMatrix) {
        for (int i : senseCodons) {
            double branchProb = 0.0;
            for (int j : senseCodons) {
                branchProb += lowerConditional[j] * probMatrix[i * GeneticCode.CODON_STATES + j];
            }
            conditionals[i] *= branchProb;
        }
    }

    private class Partial {
        public int number;
        public double[] partial;
        public double branchLength;
        public Node node;

        Partial(Node node, int number, double[] partial, double branchLength) {
            this.node = node;
            this.number = number;
            this.partial = partial;
            this.branchLength = branchLength;
        }
    }


}
