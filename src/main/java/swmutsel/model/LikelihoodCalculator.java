package swmutsel.model;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.ArrayPool;
import swmutsel.Constants;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.Pair;
import swmutsel.utils.PhyloUtils;

import java.util.*;

/**
 * Felsenstein's pruning algorithm to calculate the likelihood for codon based models. Can deal with heterogenous models.
 *
 * @author Asif Tamuri (tamuri@ebi.ac.uk)
 */
public class LikelihoodCalculator {
    private Tree tree;
    private final Map<String, Byte> states;
    private final Map<String, SubstitutionModel> cladeModels;
    private final Map<Node, String> nodeLabels = Maps.newHashMap();
    private final Map<String, TransProbCalculator> cladeCalculators;
    private final List<Node> nodeTraversalOrder = Lists.newArrayList();

    private final double[] tipConditional = new double[GeneticCode.CODON_STATES];
    private final double[] gapConditional = CoreUtils.rep(1.0d, GeneticCode.CODON_STATES);

    private final String ROOT_MODEL_NAME;

    // TODO: this is can be a constant for the entire run...
    private final int[] senseCodons = GeneticCode.getInstance().getSenseCodons();

    private double[][] probMatrix;

    private List<Double> scaledPartial = Lists.newArrayList();
    private List<Integer> scaledNodes = Lists.newArrayList();

    // TODO: Improve memory usage of these conditionals...do we need to store everything??!
    // TODO: For example, strictly speaking, we only need as many arrays as we have concurrent threads running, and each can be reused
    // Could we have a reusable pool of 64*64 arrays? We would then only create 16 (or whatever), instead of 3598!
    //  private final double[][] tipConditionals;
    // private final double[][] conditionals;
    private double[][] conditionals;

    private final List<Partial> conditionalAtRoot = Lists.newArrayList();

    public LikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, LinkedHashMap<String, SubstitutionModel> models) {
        this.tree = tree;

        setDefaultNodeTraversal();

        this.states = siteStates;
        this.cladeModels = models;
        this.cladeCalculators = Maps.newConcurrentMap();
        // For non-homogeneous models, the first added entry is the "root" model (so you should use LinkedHashMap).
        this.ROOT_MODEL_NAME = cladeModels.keySet().iterator().next();

        setNodeLabels(tree);
    }

    private void setDefaultNodeTraversal() {
        nodeTraversalOrder.clear();
        for (Node n : PhyloUtils.internalNodes(tree)) {
            nodeTraversalOrder.add(n);
        }
    }

    @SuppressWarnings("unchecked")
    public LikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, SubstitutionModel model) {
        this(tree, siteStates, CoreUtils.getLinkedHashMap(Pair.of("ALL", model)));
    }

    private void setNodeLabels(Tree tree) {
        nodeLabels.clear();
        // set up name lookup TODO: this is shared across all sites! Do it once and then pass reference!!
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

    public double updateTree(Tree tree, Map<Integer, Integer> newToOldNumbering, Set<Integer> updateNodesOriginal) {
        this.tree = tree;
        setNodeLabels(tree);

        // The node numbering for the re-rooted tree has been changed
        // Copy to new conditionals and scaledNodes lists
        double[][] newConditionals = new double[tree.getInternalNodeCount()][];
        List<Integer> newScaledNodes = Lists.newArrayList();
        for (Map.Entry<Integer, Integer> entry : newToOldNumbering.entrySet()) {
            newConditionals[entry.getValue()] = conditionals[entry.getKey()];
            if (scaledNodes.contains(entry.getKey())) {
                newScaledNodes.add(entry.getValue());
            }
        }

        scaledNodes = newScaledNodes;
        conditionals = newConditionals;

        // The nodes we need to visit to get new partial likelihoods
        nodeTraversalOrder.clear();
        Set<Integer> updateNodes = Sets.newHashSet(updateNodesOriginal);
        updateNodes.add(tree.getRoot().getNumber());

        for (int node : updateNodes) {
            // If we need to update a node that has been scaled
            if (scaledNodes.contains(node)) {
                // Remove its entry from scaled nodes/partials
                int index = scaledNodes.indexOf(node);
                scaledNodes.remove(index);
                scaledPartial.remove(0);
            }
        }

        // Order the nodes, so we visit in the right order (descendants are visited before parents)
        List<Integer> sortedNodes = Lists.newArrayList(updateNodes);
        Collections.sort(sortedNodes);
        for (int i : sortedNodes) nodeTraversalOrder.add(tree.getInternalNode(i));

        // TODO: fix - this is the same as getLogLikelihood() but without clearing scaledNodes!
        conditionalAtRoot.clear();

        cladeCalculators.clear();
        for (Map.Entry<String, SubstitutionModel> e : cladeModels.entrySet()) {
            cladeCalculators.put(e.getKey(), e.getValue().getPtCalculator());
        }

        double[] conditionals = downTree();
        double[] f = cladeModels.get(ROOT_MODEL_NAME).getCodonFrequencies();
        double sum = 0.0;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) sum += conditionals[i] * f[i];

        if (sum < 0) sum = 0;

        if (weight > -1) sum *= weight;

        double scaled = 0;
        for (double d : scaledPartial) scaled += Math.log(d);

        double lnl = Math.log(sum) + scaled;

        setDefaultNodeTraversal();

        return lnl;
    }

    public double getLikelihood(double weight) {
        this.weight = weight;
        return getLogLikelihood();
    }

    private double weight = -1;

    public double getLogLikelihood() {
        if (probMatrix == null) {
            probMatrix = new double[64][64];
            conditionals = new double[tree.getInternalNodeCount()][64];
        }

        conditionalAtRoot.clear();

        scaledPartial.clear();
        scaledNodes.clear();

        cladeCalculators.clear();
        for (Map.Entry<String, SubstitutionModel> e : cladeModels.entrySet()) {
            cladeCalculators.put(e.getKey(), e.getValue().getPtCalculator());
        }

        double[] conditionals = downTree();
        double[] f = cladeModels.get(ROOT_MODEL_NAME).getCodonFrequencies();

        double sum = 0.0;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) sum += conditionals[i] * f[i];

        if (sum < 0) sum = 0;

        if (weight > -1) sum *= weight;

        double scaled = 0;
        for (double d : scaledPartial) scaled += Math.log(d);

        // TODO: hack - clean this up
        double penalty = 0;
        if (cladeModels.get(ROOT_MODEL_NAME) instanceof SwMutSel) {
            SwMutSel model = (SwMutSel) cladeModels.get(ROOT_MODEL_NAME);
            if (model.hasFitnessFDS()) {
                penalty = model.getFDSPenalty();
            }
        }

        return Math.log(sum) + scaled + penalty;
    }

    private String getNodeLabel(Node n) {
        return nodeLabels.get(n);
    }

    public void getStorage() {
        probMatrix = ArrayPool.popStatesByStatesMatrix();
        conditionals = ArrayPool.popNodesByStatesMatrix();
    }

    public void releaseStorage() {
        // So we can keep the reference to the LikelihoodCalculator without keeping this baggage
        ArrayPool.pushStatesByStatesMatrix(probMatrix);
        ArrayPool.pushNodesByStatesMatrix(conditionals);
        probMatrix = null;
        conditionals = null;
    }

    private double[] downTree() {
        //long start = CodeTimer.start();
        for (Node node : nodeTraversalOrder) {

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

                if (node.isRoot()) {
                    conditionalAtRoot.add(new Partial(child, PhyloUtils.getNodeNumber(tree, child), Arrays.copyOf(lowerConditional, GeneticCode.CODON_STATES), child.getBranchLength()));
                }

                if (cladeModels.size() == 1) { // homogeneous model

                    // TODO: A probability matrix (of sorts) already exists in TdGCodonModel...do we need to create it again?
                    // cladeCalculators.get(ROOT_MODEL_NAME).getTransitionProbabilities(probMatrix, child.getBranchLength());
                    // updateIntraCladeConditionals(lowerConditional, conditional, probMatrix);
                    cladeCalculators.get(ROOT_MODEL_NAME).calculatePartialLikelihood(lowerConditional, conditional, child.getBranchLength());

                } else { // non-homogeneous model

                    if (getNodeLabel(node).length() == 0 // the root of the tree is a parent without a label
                            || getNodeLabel(child).substring(0, 2).equals(getNodeLabel(node).substring(0, 2))) { // or we're not switching to a different model

                        cladeCalculators.get(getNodeLabel(child).substring(0, 2)).getTransitionProbabilities(probMatrix, child.getBranchLength());
                        updateIntraCladeConditionals(lowerConditional, conditional, probMatrix);

                    } else { // this is a hostshift!
                        final double[][] probMatrix1 = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
                        cladeCalculators.get(getNodeLabel(node).substring(0, 2)).getTransitionProbabilities(probMatrix, child.getBranchLength() * Constants.CLADE_BRANCH_SPLIT);
                        cladeCalculators.get(getNodeLabel(child).substring(0, 2)).getTransitionProbabilities(probMatrix1, child.getBranchLength() * (1 - Constants.CLADE_BRANCH_SPLIT));
                        updateInterCladeConditionals(lowerConditional, conditional, probMatrix, probMatrix1);
                    }
                }
            }

            // don't scale the root node
            if (!node.isRoot()) {
                scaleConditionals(node, conditional);
            }

            conditionals[node.getNumber()] = conditional;

        }

        return conditionals[tree.getRoot().getNumber()];
    }

    public double getLogLikelihoodForBranch(int node, double branchLength) {
        // TODO: store the true branchlengths here, so they can be updated with the distributed runner

        double[] sumPartial = new double[GeneticCode.CODON_STATES];
        Arrays.fill(sumPartial, 1.0);

        probMatrix = ArrayPool.popStatesByStatesMatrix();
        for (Partial p : conditionalAtRoot) {
            if (p.number == node) {
                // cladeCalculators.get(ROOT_MODEL_NAME).getTransitionProbabilities(probMatrix, branchLength);
                cladeCalculators.get(ROOT_MODEL_NAME).calculatePartialLikelihood(p.partial, sumPartial, branchLength);
            } else {
                // cladeCalculators.get(ROOT_MODEL_NAME).getTransitionProbabilities(probMatrix, p.branchLength);
                cladeCalculators.get(ROOT_MODEL_NAME).calculatePartialLikelihood(p.partial, sumPartial, p.branchLength);
            }
            // updateIntraCladeConditionals(p.partial, sumPartial, probMatrix);
        }

        ArrayPool.pushStatesByStatesMatrix(probMatrix);

        double lnL = 0;
        final double[] frequencies = cladeModels.get("ALL").getCodonFrequencies();
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) lnL += sumPartial[i] * frequencies[i];

        // if mixture model (i.e. weight was set when getlikelihood was called)
        if (weight > -1) lnL *= weight;

        double scaled = 0;
        for (double d : scaledPartial) scaled += Math.log(d);

        return Math.log(lnL) + scaled;
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

            // only scale if all entries in conditional are less than threshold
            for (double c : conditionals)
                if (c > Constants.SCALING_THRESHOLD)
                    return;

            for (int i = 0; i < conditionals.length; i++) {
                conditionals[i] /= Constants.SCALING_THRESHOLD;
            }
            scaledPartial.add(Constants.SCALING_THRESHOLD);
            scaledNodes.add(node.getNumber());
        }
    }

    private static void updateInterCladeConditionals(final double[] lowerConditional, final double[] conditionals, final double[][] probMatrix0, final double[][] probMatrix1) {
        final double[][] probMatrix = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            double branchProb = 0.0;
            final double[] probMatrix_r = probMatrix[i];
            final double[] probMatrix0_r = probMatrix0[i];
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                for (int k = 0; k < GeneticCode.CODON_STATES; k++) {
                    probMatrix_r[j] += probMatrix0_r[k] * probMatrix1[k][j];
                }
                branchProb += lowerConditional[j] * probMatrix_r[j];
            }
            conditionals[i] *= branchProb;
        }
    }

    private static void updateIntraCladeConditionals(final double[] lowerConditional, final double[] conditionals, final double[][] probMatrix) {
        double b1, b2, b3, b4;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            final double[] probMatrix_r = probMatrix[i];
            {
                b1 = b2 = b3 = b4 = 0;
                b1 += lowerConditional[0] * probMatrix_r[0];
                b1 += lowerConditional[1] * probMatrix_r[1];
                b1 += lowerConditional[2] * probMatrix_r[2];
                b1 += lowerConditional[3] * probMatrix_r[3];
                b1 += lowerConditional[4] * probMatrix_r[4];
                b1 += lowerConditional[5] * probMatrix_r[5];
                b1 += lowerConditional[6] * probMatrix_r[6];
                b1 += lowerConditional[7] * probMatrix_r[7];
                b1 += lowerConditional[8] * probMatrix_r[8];
                b1 += lowerConditional[9] * probMatrix_r[9];
                b1 += lowerConditional[10] * probMatrix_r[10];
                b1 += lowerConditional[11] * probMatrix_r[11];
                b1 += lowerConditional[12] * probMatrix_r[12];
                b1 += lowerConditional[13] * probMatrix_r[13];
                b1 += lowerConditional[14] * probMatrix_r[14];
                b1 += lowerConditional[15] * probMatrix_r[15];
                b2 += lowerConditional[16] * probMatrix_r[16];
                b2 += lowerConditional[17] * probMatrix_r[17];
                b2 += lowerConditional[18] * probMatrix_r[18];
                b2 += lowerConditional[19] * probMatrix_r[19];
                b2 += lowerConditional[20] * probMatrix_r[20];
                b2 += lowerConditional[21] * probMatrix_r[21];
                b2 += lowerConditional[22] * probMatrix_r[22];
                b2 += lowerConditional[23] * probMatrix_r[23];
                b2 += lowerConditional[24] * probMatrix_r[24];
                b2 += lowerConditional[25] * probMatrix_r[25];
                b2 += lowerConditional[26] * probMatrix_r[26];
                b2 += lowerConditional[27] * probMatrix_r[27];
                b2 += lowerConditional[28] * probMatrix_r[28];
                b2 += lowerConditional[29] * probMatrix_r[29];
                b2 += lowerConditional[30] * probMatrix_r[30];
                b2 += lowerConditional[31] * probMatrix_r[31];
                b3 += lowerConditional[32] * probMatrix_r[32];
                b3 += lowerConditional[33] * probMatrix_r[33];
                b3 += lowerConditional[34] * probMatrix_r[34];
                b3 += lowerConditional[35] * probMatrix_r[35];
                b3 += lowerConditional[36] * probMatrix_r[36];
                b3 += lowerConditional[37] * probMatrix_r[37];
                b3 += lowerConditional[38] * probMatrix_r[38];
                b3 += lowerConditional[39] * probMatrix_r[39];
                b3 += lowerConditional[40] * probMatrix_r[40];
                b3 += lowerConditional[41] * probMatrix_r[41];
                b3 += lowerConditional[42] * probMatrix_r[42];
                b3 += lowerConditional[43] * probMatrix_r[43];
                b3 += lowerConditional[44] * probMatrix_r[44];
                b3 += lowerConditional[45] * probMatrix_r[45];
                b3 += lowerConditional[46] * probMatrix_r[46];
                b3 += lowerConditional[47] * probMatrix_r[47];
                b4 += lowerConditional[48] * probMatrix_r[48];
                b4 += lowerConditional[49] * probMatrix_r[49];
                b4 += lowerConditional[50] * probMatrix_r[50];
                b4 += lowerConditional[51] * probMatrix_r[51];
                b4 += lowerConditional[52] * probMatrix_r[52];
                b4 += lowerConditional[53] * probMatrix_r[53];
                b4 += lowerConditional[54] * probMatrix_r[54];
                b4 += lowerConditional[55] * probMatrix_r[55];
                b4 += lowerConditional[56] * probMatrix_r[56];
                b4 += lowerConditional[57] * probMatrix_r[57];
                b4 += lowerConditional[58] * probMatrix_r[58];
                b4 += lowerConditional[59] * probMatrix_r[59];
                b4 += lowerConditional[60] * probMatrix_r[60];
                b4 += lowerConditional[61] * probMatrix_r[61];
                b4 += lowerConditional[62] * probMatrix_r[62];
                b4 += lowerConditional[63] * probMatrix_r[63];
            }
            conditionals[i] *= b1 + b2 + b3 + b4;
        }
    }

    static private class Partial {
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

