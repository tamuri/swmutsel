package swmutsel;

import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;
import pal.alignment.Alignment;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.model.FMutSel0;
import swmutsel.model.SubstitutionModel;
import swmutsel.model.parameters.*;
import swmutsel.runner.MultiThreadedRunner;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.PhyloUtils;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 07/12/2013 17:08
 */
public class TestHessianSKT2004 {




    public static void main(String[] args) {
        TestHessianSKT2004 t = new TestHessianSKT2004();
        t.run();
    }

    private void pamlToSwMutSel() {


        Pattern p2 = Pattern.compile("\\s+");
        Iterable<String> f = Splitter.on(p2).split("-1.196811 -1.196811 -1.933164 -1.933164 -2.704113 -2.704113 -2.704113 -2.704113 -2.905748 -2.905748 -2.449835 -2.449835 -1.049355 -1.049355 -1.933164 -1.933164 -1.933164 -1.933164 -3.540852 -3.540852 -3.540852 -3.540852 -3.918200 -3.918200 -4.168093 -4.168093 -2.937931 -2.937931 -2.937931 -2.937931 -1.957810 -1.957810 -2.429621 -2.429621 -3.351957 -3.351957 -3.351957 -3.351957 -3.800958 -3.800958 -4.478663 -4.478663 -2.704113 -2.704113 -1.368296 -1.368296 -1.368296 -1.368296 -1.562600 -1.562600 -1.562600 -1.562600 -2.591327 -2.591327 -2.472289 -2.472289  0.000000  0.000000  0.000000  0.000000");
        List<String> ff = Lists.newArrayList(f);


        Map<Integer, Double> fits = Maps.newHashMap();

        int pos = 0;
        for (int k = 0; k < GeneticCode.CODON_STATES; k++) {

            if (GeneticCode.getInstance().isSenseCodon(k)) {

                System.out.printf("%s %s\n",
                        GeneticCode.getInstance().getAminoAcidCharByIndex(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(k)),
                        // GeneticCode.getInstance().getCodonTLA(k),
                        ff.get(pos)
                );

                fits.put(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(k), Double.parseDouble(ff.get(pos)));


                pos++;
            }

        }

        System.out.println();

        System.out.printf("%s", fits.get(0));
        for (int k = 1; k < 20; k++) {
            // System.out.printf("%s %s %s\n", k, GeneticCode.getInstance().getAminoAcidCharByIndex(k), fits.get(k));
            System.out.printf(",%s", fits.get(k));
        }
        System.out.println();

        System.exit(0);


    }

    private void run() {

        GeneticCode.setCode(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);


        final SubstitutionModel model;
        MultiThreadedRunner r;
        Tree t;

        // t = PhyloUtils.readTree("/Users/atamuri/Documents/2013/mutselbranch/data/small4s/small4s_mtCDNA_500.tree");
        t = PhyloUtils.readTree("/Users/atamuri/Documents/2013/mutselbranch/data/mdr.0226/tree.fmutsel0.txt");
        MatrixArrayPool.treeSize = t.getInternalNodeCount();
        // Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2013/mutselbranch/data/small4s/small4s_mtCDNA_500.txt");
        Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2013/mutselbranch/data/mdr.0226/aln.phy");

        ArgumentsProcessor p = new ArgumentsProcessor();
        Table<String, Integer, Byte> table = p.getCleanedSitesTable(a);

        r = new MultiThreadedRunner(3);
        r.setSites(table);

/*        model = new FMutSel0(
                new TsTvRatio(7.96527),
                new Omega(0.05327),
                new BaseFrequencies(new double[]{0.17572, 0.28418, 0.45954, 0.08056}),
                new Fitness(
                        new int[]{7, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
                        new double[]{-0.978738, -2.741638, -3.007054, -2.449005, -25.90481, -3.087357, -2.614215, 0.0, -3.078092, -1.279761, -1.310963, -3.866942, -1.893621, -0.793491, -2.537767, -2.17783, -2.577635, -1.070589, -3.024342, -0.95129}))*/;

        model = new FMutSel0(
                new TsTvRatio(3.24782),
                new Omega(0.09293),
                new BaseFrequencies(new double[]{0.15494, 0.31435, 0.47675, 0.05396}),
                new Fitness(
                        new int[]{7, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
                        new double[]{-1.5626,-2.937931,-3.800958,-2.591327,-2.449835,-4.168093,-2.472289,0.0,-3.9182,-1.95781,-1.933164,-4.478663,-2.429621,-1.196811,-3.540852,-2.704113,-3.351957,-1.049355,-2.905748,-1.368296}
                )

        );
        //


        model.build();

        final List<Node> nodes = Lists.newArrayList();
        traverse(t.getRoot(), nodes);

        final double[] branchMLE = getTreeBranchLengths(nodes);
        double[] paramMLE = Doubles.concat(branchMLE, Mapper.getOptimisable(model.getParameters()));
        double[] H = new double[paramMLE.length * paramMLE.length];
        double[] g = new double[paramMLE.length * paramMLE.length];

        calculateHessian(paramMLE, g, H, new UpdateCallback() {
            @Override
            public void update(double[] point) {
                double[] branchLengths = Arrays.copyOfRange(point, 0, branchMLE.length);
                setBranchLengths(branchLengths, nodes);
                double[] modelParameters = Arrays.copyOfRange(point, branchMLE.length, point.length);
                Mapper.setOptimisable(model.getParameters(), modelParameters);
                model.build();
            }
        }, model, t, r);

        System.out.println();
        System.out.printf(" %s\n", t.getExternalNodeCount());
        System.out.println();
        System.out.println(t.toString());
        System.out.println();


        for (double b : branchMLE) {
            System.out.printf(" %9.6f", b);
        }
        System.out.println();
        System.out.println();

        for (int i = 0; i < branchMLE.length; i++) {
            if (branchMLE[i] > 0.0004 && Math.abs(g[i]) < 0.005) g[i] = 0;
            System.out.printf(" %9.6f", g[i]);
        }
        System.out.println();
        System.out.println();


        System.out.println("Hessian");
        System.out.println();

        for (int i = 0; i < branchMLE.length; i++) {
            for (int j = 0; j < branchMLE.length; j++) {
                System.out.printf("%10.4g\t", H[i * paramMLE.length + j]);
            }
            System.out.println();
        }

        // See Yang (2006) CME, pg. 24, eq. 1.46:
        // \hat \theta \~ {N_k}(\theta , - H{}^{ - 1})
        // H = \{ {d^2}\ell /d{\theta _i}d{\theta _j}) \}

        // DoubleMatrix2D HInv = Algebra.DEFAULT.inverse(DoubleFactory2D.dense.make(H, paramMLE.length).assign(Functions.neg));

        /*System.out.println();
        System.out.println("SEs for parameters:");
        for (int i = 0; i < paramMLE.length; i++)
            System.out.printf(" %8.6f", (HInv.get(i, i) > 0. ? Math.sqrt(HInv.get(i, i)) : -1));
        System.out.println();*/


        r.shutdown();


        /*
        4

        (Diceros_bicornis:0.3234100,(Dicerorhinus_sumatrensis:0.4157900,Rhinoceros_unicornis:0.3019100):0.0319500,Equus_caballus:0.7360200);

        0.323410  0.415790  0.301910  0.736020  0.031950

        Hessian

        -465.1	    -40.65	    -122.5	    -35.90	    -279.2
        -40.65	    -393.0	    -90.16	    -31.80	    -154.2
        -122.5	    -90.16	    -590.5	    -40.11	    -185.6
        -35.90	    -31.80	    -40.11	    -181.7	    -71.47
        -279.2	    -154.2	    -185.6	    -71.47	     -1567

        SEs for parameters:
        0.053400 0.057508 0.044186 0.091907 0.027904 0.253167 0.417119 0.615201 1.005567 0.006400 0.281648 0.414777 0.358341 0.410100 10100692.062774 0.365875 0.408004 0.401932 0.329370 0.314518 0.427016 0.328143 0.367040 0.359651 0.322330 0.325553 0.351849 0.446997 0.285477
        */
    }

    private void setBranchLengths(double[] branchLengths, List<Node> nodes) {
        // branchLengths variable expects branch lengths in same order as given by getTreeBranchLengths().
        int pos = 0;

        for (Node n : nodes) {
            n.setBranchLength(branchLengths[pos++]);
        }
    }

    private double[] getTreeBranchLengths(List<Node> nodes) {
        double[] branchLengths = new double[nodes.size()];

        int i = 0;
        for (Node n : nodes) {
            branchLengths[i++] = n.getBranchLength();
        }

        return branchLengths;
    }


    private void traverse(Node node, List<Node> storage) {
        for (int i = 0; i < node.getChildCount(); i++) {
            storage.add(node.getChild(i));
            traverse(node.getChild(i), storage);
        }
    }


    /**
     * Based on Ziheng's implementation of Seo et al. (2004).
     * See function HessianSKT2004() in treesub.c of PAML.
     */
    private void calculateHessian(double MLEs[], double g[], double H[], UpdateCallback callback, SubstitutionModel model, Tree t, MultiThreadedRunner r) {
        // MLEs -> maximum likelihood estimates for each parameter, branch lengths and model parameters concatenated
        // g[] -> array to hold gradient
        // H[]  -> array to hold Hessian
        // UpdateCallback -> to update model with new parameter values

        int siteCount = r.getSites().columnKeySet().size(); // number of site patterns (or sites)

        double[] x = new double[MLEs.length];
        double[][] lnL = new double[2][MLEs.length];
        double[][][] df = new double[2][MLEs.length][siteCount];

        double EPSILON = 1e-6;
        double eh0 = EPSILON * 2;
        double SMALL = 4e-6;
        double eh;

        int nearZero = 0;

        Arrays.fill(H, 0);

        for (int bracket = 0; bracket < 2; bracket++) {
            for (int i = 0; i < MLEs.length; i++) {

                System.arraycopy(MLEs, 0, x, 0, MLEs.length);

                eh = eh0 * (Math.abs(MLEs[i]) + 1);

                if (bracket == 0) x[i] = MLEs[i] - eh;
                else x[i] = MLEs[i] + eh;

                if (Math.abs(x[i]) <= SMALL) nearZero++;

                callback.update(x);

                // Calculate the log-likelihood. Save the site-wise likelihoods, as well as the total
                Map<Integer, Double> out = r.getLogLikelihood(t, model);
                for (int j = 0; j < siteCount; j++) df[bracket][i][j] = out.get(j + 1);
                lnL[bracket][i] = CoreUtils.sum(out.values());
            }
        }

        for (int i = 0; i < MLEs.length; i++) {
            eh = eh0 * (Math.abs(MLEs[i]) + 1);
            g[i] = (lnL[1][i] - lnL[0][i]) / (eh * 2);
        }


        for (int i = 0; i < MLEs.length; i++) {
            eh = eh0 * (Math.abs(MLEs[i]) + 1);

            for (int k = 0; k < siteCount; k++) {
                df[0][i][k] = (df[1][i][k] - df[0][i][k]) / (eh * 2);
            }
        }


        for (int i = 0; i < MLEs.length; i++) {
            for (int j = 0; j < MLEs.length; j++) {
                for (int k = 0; k < siteCount; k++) {
                    H[i * MLEs.length + j] -= df[0][i][k] * df[0][j][k];
                }
            }
        }

        if (nearZero > 0) System.out.println("\nWarning: Hessian matrix may be unreliable for zero branch lengths\n");

    }


    interface UpdateCallback {
        public void update(double[] point);
    }

}
