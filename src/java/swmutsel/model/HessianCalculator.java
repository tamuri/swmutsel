package swmutsel.model;

import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.model.parameters.Mapper;
import swmutsel.runner.Runner;
import swmutsel.utils.CoreUtils;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 10/03/2014 13:24
 */
public class HessianCalculator {

    private final Tree tree;
    private final int siteCount;
    private final Runner runner;

    public HessianCalculator(Tree tree, int siteCount, Runner runner) {
        this.tree = tree;
        this.siteCount = siteCount;
        this.runner = runner;
    }

    public void get(final SubstitutionModel model) {
        final List<Node> nodes = Lists.newArrayList();
        traverse(tree.getRoot(), nodes);

        final double[] branchMLE = getTreeBranchLengths(nodes);
        double[] xmle = Doubles.concat(branchMLE, Mapper.getOptimisable(model.getParameters()));

        double[][] H = new double[xmle.length][xmle.length];
        double[] g = new double[xmle.length];

        calculateHessian(xmle, g, H,
                new UpdateModelCallback() {
                    @Override
                    public void update(double[] point) {
                        double[] branchLengths = Arrays.copyOfRange(point, 0, branchMLE.length);
                        setBranchLengths(branchLengths, nodes);
                        double[] modelParameters = Arrays.copyOfRange(point, branchMLE.length, point.length);
                        Mapper.setOptimisable(model.getParameters(), modelParameters);
                        model.build();
                    }
                },
                new LogLikelihoodCallback() {
                    @Override
                    public Map<Integer, Double> get() {
                        return runner.getLogLikelihood(model);
                    }
                }
        );

        printRST2(branchMLE, g, H);

    }

    public void get(final FitnessStore fitnesses, final SwMut mutation, final Penalty penalty) {
        final List<Node> nodes = Lists.newArrayList();
        traverse(tree.getRoot(), nodes);

        final double[] branchMLE = getTreeBranchLengths(nodes);
        double[] xmle = Doubles.concat(branchMLE, Mapper.getOptimisable(mutation.getParameters()));


        // Convert the MLE of branch lengths from neutral substitutions/sites to substitutions/codon
        double averageSubstitutionRate = runner.getAverageSubstitutionRate(mutation, fitnesses);
        CoreUtils.msg("Average substitution rate per site = %.7f\n", averageSubstitutionRate);

        double scalingMLE = mutation.getBranchScaling().get();
        double[] subsPerSiteBranchLengths = getTreeBranchLengths(nodes);
        for (int i = 0; i < subsPerSiteBranchLengths.length; i++) {
            subsPerSiteBranchLengths[i] *= scalingMLE;
            subsPerSiteBranchLengths[i] *= averageSubstitutionRate;
        }
        setBranchLengths(subsPerSiteBranchLengths, nodes);
        System.out.println(tree.toString());

        setBranchLengths(branchMLE, nodes);

        double[][] H = new double[xmle.length][xmle.length];
        double[] g = new double[xmle.length];

        calculateHessian(xmle, g, H,
                new UpdateModelCallback() {
                    @Override
                    public void update(double[] point) {
                        double[] branchLengths = Arrays.copyOfRange(point, 0, branchMLE.length);
                        setBranchLengths(branchLengths, nodes);
                        double[] modelParameters = Arrays.copyOfRange(point, branchMLE.length, point.length);
                        Mapper.setOptimisable(mutation.getParameters(), modelParameters);
                        mutation.build();
                    }
                },
                new LogLikelihoodCallback() {
                    @Override
                    public Map<Integer, Double> get() {
                        return runner.getLogLikelihood(tree, mutation, fitnesses, penalty);
                    }
                }
        );



        setBranchLengths(subsPerSiteBranchLengths, nodes);

        for (int i = 0; i < g.length; i++) {
            g[i] /= averageSubstitutionRate;
        }

        for (int i = 0; i < subsPerSiteBranchLengths.length; i++) {
            for (int j = 0; j < subsPerSiteBranchLengths.length; j++) {
                H[i][j] /= Math.pow(averageSubstitutionRate, 2);
            }
        }

        printRST2(subsPerSiteBranchLengths, g, H);

    }

    private void printRST2(double[] branchLengths, double[] g, double[][] H) {
        System.out.println();
        System.out.printf(" %s\n", tree.getExternalNodeCount());
        System.out.println();

        System.out.println(tree.toString());
        System.out.println();

        for (double b : branchLengths) {
            System.out.printf(" %9.6f", b);
        }
        System.out.println();
        System.out.println();

        for (int i = 0; i < branchLengths.length; i++) {
            if (branchLengths[i] > 0.0004 && Math.abs(g[i]) < 0.005) g[i] = 0;
            System.out.printf(" %9.6f", g[i]);
        }
        System.out.println();
        System.out.println();


        System.out.println("Hessian");
        System.out.println();

        for (int i = 0; i < branchLengths.length; i++) {
            for (int j = 0; j < branchLengths.length; j++) {
                System.out.printf("%10.4g\t", H[i][j]);
            }
            System.out.println();
        }
    }

    /**
     * Based on Ziheng's implementation of Seo et al. (2004).
     * See function HessianSKT2004() in treesub.c of PAML.
     */
    private void calculateHessian(double xmle[], double g[], double H[][], UpdateModelCallback updateModelCallback, LogLikelihoodCallback logLikelihoodCallback) {
        // MLEs -> maximum likelihood estimates for each parameter, branch lengths and model parameters concatenated
        // g[] -> array to hold gradient
        // H[]  -> array to hold Hessian
        // UpdateModelCallback -> to update model with new parameter values

        double[] x = new double[xmle.length];
        double[][] lnL = new double[2][xmle.length];
        double[][][] df = new double[2][xmle.length][siteCount];

        double EPSILON = 1e-6;
        double eh0 = EPSILON * 2;
        double SMALL = 4e-6;
        double eh;

        int nearZero = 0;

        for (int backforth = 0; backforth < 2; backforth++) {
            for (int i = 0; i < xmle.length; i++) {

                System.arraycopy(xmle, 0, x, 0, xmle.length);

                eh = eh0 * (Math.abs(xmle[i]) + 1);

                if (backforth == 0) x[i] = xmle[i] - eh;
                else x[i] = xmle[i] + eh;

                if (Math.abs(x[i]) <= SMALL) nearZero++;

                updateModelCallback.update(x);

                // Calculate the log-likelihood. Save the site-wise likelihoods, as well as the total
                //Map<Integer, Double> out = runner.get(tree, mutation, fitnesses, penalty);
                Map<Integer, Double> out = logLikelihoodCallback.get();
                for (int j = 0; j < siteCount; j++) {
                    df[backforth][i][j] = -out.get(j + 1); // only using negative because that's what ZY has. I don't think it matters!
                }

                lnL[backforth][i] = CoreUtils.sum(out.values());
            }
        }

        for (int i = 0; i < xmle.length; i++) {
            eh = eh0 * (Math.abs(xmle[i]) + 1);
            g[i] = (lnL[1][i] - lnL[0][i]) / (eh * 2);
        }


        for (int i = 0; i < xmle.length; i++) {
            eh = eh0 * (Math.abs(xmle[i]) + 1);

            for (int h = 0; h < siteCount; h++) {
                df[0][i][h] = (df[1][i][h] - df[0][i][h]) / (eh * 2);
            }
        }

        for (int i = 0; i < xmle.length; i++) {
            for (int j = 0; j < xmle.length; j++) {
                for (int h = 0; h < siteCount; h++) {
                    H[i][j] -= df[0][i][h] * df[0][j][h];
                }
            }
        }

        if (nearZero > 0) System.out.println("\nWarning: Hessian matrix may be unreliable for zero branch lengths\n");

    }

    private void setBranchLengths(double[] branchLengths, List<Node> nodes) {
        // Branch lengths should be in same order as given by getTreeBranchLengths().
        int i = 0;
        for (Node n : nodes) n.setBranchLength(branchLengths[i++]);
    }

    private double[] getTreeBranchLengths(List<Node> nodes) {
        double[] branchLengths = new double[nodes.size()];
        int i = 0;
        for (Node n : nodes) branchLengths[i++] = n.getBranchLength();
        return branchLengths;
    }

    private void traverse(Node node, List<Node> storage) {
        for (int i = 0; i < node.getChildCount(); i++) {
            storage.add(node.getChild(i));
            traverse(node.getChild(i), storage);
        }
    }

    private static interface UpdateModelCallback {
        public void update(double[] point);
    }

    private static interface LogLikelihoodCallback {
        public Map<Integer, Double> get();
    }
}
