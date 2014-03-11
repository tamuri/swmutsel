package swmutsel;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.internal.Maps;
import com.google.common.base.Splitter;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.collect.Table;
import com.google.common.io.Files;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleValueChecker;
import pal.alignment.Alignment;
import pal.alignment.DataTranslator;
import pal.datatype.Codons;
import pal.datatype.SpecificAminoAcids;
import pal.tree.*;
import swmutsel.model.SwMut;
import swmutsel.model.SwMutSel;
import swmutsel.model.parameters.*;
import swmutsel.runner.MultiThreadedRunner;
import swmutsel.trees.RerootedTreeIterator;
import swmutsel.utils.*;

import java.io.File;
import java.io.PushbackReader;
import java.io.StringReader;
import java.nio.charset.Charset;
import java.util.*;
import java.util.regex.Pattern;

/**
 * Author: Asif Tamuri (atamuri@ebi.ac.uk)
 * Date: 14/10/2013 16:50
 */
public class Test {
    public static void main(String[] args) throws Exception{


        calculatelambda();

        // System.out.println("\u2190 \u2191 \u2192 \u2193 \u21e1");
        //makephylip();

        //System.out.println(CoreUtils.calculateMachineEpsilonDouble());
        //testcolt();
        //testtree();
        //testbrent();

        //testTreeRerooting();
        //testJcommander();
        // singleSubstitutionCodons();
        //alignmentTranslate();
        //testfmutsel();
// testpatterns();
//         testS();

        //testswmutsel();


//        treewithsequence();


/*

          tree = PhyloUtils.readTree("/Users/atamuri/Documents/2013/mutselbranch/data/small4s/codeml/tree.txt");
        //Tree t = PhyloUtils.readTree("/Users/atamuri/Documents/2013/tdg12/etc/all.but.ND6.AUT.FMutSel0.tree");

        MatrixArrayPool.treeSize = tree.getInternalNodeCount();
        final Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2013/mutselbranch/data/small4s/small4s_mtCDNA_500.txt");
        //Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2013/tdg12/etc/all.but.ND6.AUT.phys");



        ArgumentsProcessor p = new ArgumentsProcessor();
        Table<String, Integer, Byte> rawtable = p.getCleanedSitesTable(a);

        table = GuavaUtils.getColumnSubset(rawtable, Ints.asList(CoreUtils.seqi(1, 311)));



        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            //System.out.printf("%s\t", i);

            tipConditionals[i * GeneticCode.CODON_STATES + i] = 1;
            tipConditionals[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES + i] = 1; // gap is row 64

        }


        // only needs to be done once
        internalNodeChildCount = new byte[tree.getInternalNodeCount()];

        int total = 0;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            internalNodeChildCount[i] = (byte) tree.getInternalNode(i).getChildCount();
            total += internalNodeChildCount[i];
        }

        internalNodeChildNumber = new int[total];
        internalNodeChildBranchLength = new double[total];
        childLeafOffset = new int[total];

        int pos = 0;

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            for (int j = 0; j < tree.getInternalNode(i).getChildCount(); j++) {
                if (!tree.getInternalNode(i).getChild(j).isLeaf()) { // not leaf, internal
                    internalNodeChildNumber[pos] = tree.getInternalNode(i).getChild(j).getNumber();
                    internalNodeChildBranchLength[pos] = tree.getInternalNode(i).getChild(j).getBranchLength();
                    childLeafOffset[pos] = tree.getInternalNode(i).getChild(j).getNumber() * 64;
                }

                pos++;
            }
        }

        newConditionals= new double[tree.getInternalNodeCount() * GeneticCode.CODON_STATES];




        model = new FMutSel0(
                new TsTvRatio(7.6),
                new Omega(0.001),
                new BaseFrequencies(new double[]{1,1,1,1}), // initial pi are observed frequencies
                new Fitness(
                        new int[]{7, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
                        new double[]{-0.978738, -2.741638, -3.007054, -2.449005, -25.90481, -3.087357, -2.614215, 0.0, -3.078092, -1.279761, -1.310963, -3.866942, -1.893621, -0.793491, -2.537767, -2.17783, -2.577635, -1.070589, -3.024342, -0.95129}));
        model.build();


        long start, end;

        for (int i = 0; i < 100000; i++) {
            start = System.currentTimeMillis();
            fastcalculator();
            System.out.printf("in %s\n", System.currentTimeMillis() - start);
        }
*/

    }

    private static void calculatelambda() {
        GeneticCode.setCode(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);
        SwMut mut = new SwMut(new TsTvRatio(TsTvRatio.getDefault()), new BranchScaling(BranchScaling.getDefault()), new BaseFrequencies(BaseFrequencies.getDefault()));
        mut.build();
        double[] f = CoreUtils.repd(0, 20);
        SwMutSel mutsel = new SwMutSel(mut, new Fitness(f));
        System.out.println(mutsel.getExpectedSubsPerSite());
    }

    private static void makephylip() throws Exception {
        String path = "/Users/atamuri/bin/paml4.7/examples/DatingSoftBound/mtCDNApri123_toparse.txt";

        Map<String, String> one = Maps.newHashMap();
        Map<String, String> two = Maps.newHashMap();
        Map<String, String> three = Maps.newHashMap();

        List<String> lines = Files.readLines(new File(path), Charset.defaultCharset());

        Pattern p = Pattern.compile("\\s+");

        for (String line : lines) {

            Iterable<String> split = Splitter.on(p).split(line);

            Iterator<String> it = split.iterator();

            String name = it.next();
            String sequence = it.next();


            if (!one.containsKey(name)) {
                one.put(name, sequence);
                continue;
            }

            if (!two.containsKey(name)) {
                two.put(name, sequence);
                continue;

            }

            if (!three.containsKey(name)) {
                three.put(name, sequence);
            }

        }


        System.out.printf("%s\n", one.size());
        System.out.printf("%s\n", two.size());
        System.out.printf("%s\n", three.size());



        /*
        for (String n : one.keySet()) {

            System.out.printf("%s      ", n);

            char[] c1 = one.get(n).toCharArray();
            char[] c2 = two.get(n).toCharArray();
            char[] c3 = three.get(n).toCharArray();

            for (int i = 0; i < c1.length; i++) {
                System.out.printf("%s%s%s", c1[i], c2[i], c3[i]);
            }

            System.out.println();

        }
*/


        Set<String> keys = one.keySet();
        List<String> orderKeys = Lists.newArrayList();
        for (String k : keys) orderKeys.add(k);


        Set<String> allPatterns = Sets.newHashSet();

        for (int i = 0; i < 3331; i++) {
            StringBuilder pattern = new StringBuilder();

            for (String k : orderKeys) {



                pattern.append(one.get(k).charAt(i));
                pattern.append(two.get(k).charAt(i));
                pattern.append(three.get(k).charAt(i));

            }

            allPatterns.add(pattern.toString());


        }


        System.out.println(allPatterns.size());





    }

    private static void testcolt() throws Exception {

        DoubleMatrix1D tipT = ColtMatrixUtils.make(1.0, 0.0, 0.0, 0.0);
        DoubleMatrix1D tipC = ColtMatrixUtils.make(0.0, 1.0, 0.0, 0.0);
        DoubleMatrix1D tipA = ColtMatrixUtils.make(0.0, 0.0, 1.0, 0.0);
        DoubleMatrix1D tipG = ColtMatrixUtils.make(0.0, 0.0, 0.0, 1.0);

        DoubleMatrix2D PtLeft, PtRight;


        double kappa = 2.0;

        // K80 from Ziheng's book pg. 103, t=0.2
        PtLeft = getK80(0.2, kappa);
        PtRight = getK80(0.1, kappa);

        System.out.println(PtLeft);
        System.out.println(tipT);
        DoubleMatrix1D x = Algebra.DEFAULT.mult(PtLeft, tipT);

        System.out.println(x);
        DoubleMatrix1D y = Algebra.DEFAULT.mult(PtLeft, tipC);




        DoubleMatrix1D n7 = mult(x, y);
        //System.out.println(n7);


        DoubleMatrix1D n6 = mult(Algebra.DEFAULT.mult(PtRight, n7), Algebra.DEFAULT.mult(PtLeft, tipA));
        System.out.println(n6);

        DoubleMatrix1D n8 = mult(Algebra.DEFAULT.mult(PtLeft, tipC), Algebra.DEFAULT.mult(PtLeft, tipC));
        System.out.println(n8);

        DoubleMatrix1D n0 = mult(
                Algebra.DEFAULT.mult(PtRight, n6),
                Algebra.DEFAULT.mult(PtRight, n8)
        );
        System.out.println(n0);

        n0.assign(Functions.mult(0.25));
        System.out.println(n0.zSum());
        System.out.println(Math.log(n0.zSum()));


        Tree tree = new ReadTree(new PushbackReader(new StringReader("(((A:0.2, B:0.2):0.1, C:0.2):0.1 , (D:0.2, E:0.2):0.1);   ")));



        for (Node n : PhyloUtils.externalNodes(tree)) {
            System.out.printf("%s\t%s\t%s\n", n.getIdentifier().getName(), n.getBranchLength(), n.getParent());
        }


        for (Node n : PhyloUtils.internalNodes(tree)) {
            System.out.printf("%s\t%s\t%s\n", n.getIdentifier().getName(), n.getBranchLength(), n.getParent());
        }


        for (Tree tr : TreeManipulator.getEveryRoot(tree)) {
            System.out.println(tr);

            Map<Node, DoubleMatrix1D> partials = Maps.newHashMap();
            for (int i = 0; i < tr.getInternalNodeCount(); i++) {
                partials.put(tr.getInternalNode(i), DoubleFactory1D.dense.make(4, 0));
            }


            for (int i = 0; i < tr.getInternalNodeCount(); i++) {
                Node node = tr.getInternalNode(i);

                Node left = node.getChild(0);
                Node right = node.getChild(1);

                DoubleMatrix1D r, l;
                DoubleMatrix2D lPt, rPt;

                if (left.isLeaf()) {
                    char id = left.getIdentifier().getName().charAt(0);
                    switch (id) {
                        case 'A':
                            l = tipT;
                            break;
                        case 'C':
                            l = tipA;
                            break;
                        default:
                            l = tipC;
                    }

                } else {
                    l = partials.get(left);
                }
                lPt = getK80(left.getBranchLength(), kappa);


                if (right.isLeaf()) {
                    char id = right.getIdentifier().getName().charAt(0);
                    switch (id) {
                        case 'A':
                            r = tipT;
                            break;
                        case 'C':
                            r = tipA;
                            break;
                        default:
                            r = tipC;
                    }
                } else {
                    r = partials.get(right);
                }


                rPt = getK80(right.getBranchLength(), kappa);

                DoubleMatrix1D p = mult(Algebra.DEFAULT.mult(lPt, l), Algebra.DEFAULT.mult(rPt, r));
                partials.get(node).assign(p);


                // System.out.printf("%s\t%s\n", i, partials.get(n));


            }


            DoubleMatrix1D root = partials.get(tr.getRoot());
            root.assign(Functions.mult(0.25));
            System.out.println(Math.log(root.zSum()));


            System.out.println();
            System.out.println();

            //break;
        }


    }

    private static DoubleMatrix2D getJC(double d) {
        // Jukes-Cantor
        double p0 = 0.25 + (0.75 * Math.exp(-4 * d));
        double p1 = 0.25 - (0.25 * Math.exp(-4 * d));

        DoubleMatrix2D Pt = DoubleFactory2D.dense.make(4, 4);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j) {
                    Pt.set(i, j, p0);
                } else {
                    Pt.set(i, j, p1);
                }
            }
        }

        return Pt;
    }

    /**
     * @param d distance
     * @param k kappa - ts/tv ratio
     * @return transition probabilities
     */
    private static DoubleMatrix2D getK80(double d, double k) {
        DoubleMatrix2D Pt = DoubleFactory2D.dense.make(4, 4);

        double p0 = 0.25 + (0.25 * Math.exp(-4 * d / (k + 2))) + (0.5 * Math.exp(-2 * d * (k + 1) / (k + 2)));
        double p1 = 0.25 + (0.25 * Math.exp(-4 * d / (k + 2))) - (0.5 * Math.exp(-2 * d * (k + 1) / (k + 2)));
        double p2 = 0.25 - (0.25 * Math.exp(-4 * d / (k + 2)));

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {

                if (i == j) {
                    Pt.set(i, j, p0);
                } else if (i == 0 && j == 1 || i == 1 && j == 0 || i == 2 && j == 3 || i == 3 && j == 2) {
                    Pt.set(i, j, p1);
                } else {
                    Pt.set(i, j, p2);

                }
            }
        }

        return Pt;
    }


    private static DoubleMatrix1D mult(DoubleMatrix1D x, DoubleMatrix1D y) {
        DoubleMatrix1D z = DoubleFactory1D.dense.make(x.size());
        z.assign(x);
        z.assign(y, Functions.mult);
        return z;
    }

    private static void testbrent() throws Exception {
        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);

        double one = -56318.16538729038;
        double two = -56318.15065647167;

        PointValuePair previous = new PointValuePair(null, one);
        PointValuePair next = new PointValuePair(null, two);




        System.out.println(convergenceChecker.converged(100, previous, next));


    }

    private static void testtree() throws Exception {
        GeneticCode.setCode(GeneticCode.STANDARD_CODE);


        //      (((x,y),(z,s)),(t,u))

        Tree t = new ReadTree(new PushbackReader(new StringReader("( ((x,y),(z,s)), ((t,u),(v,w)) );")));
        t = TreeManipulator.getUnrooted(t);
/*

        System.out.println(t.getInternalNodeCount());
        System.out.println(t.getExternalNodeCount());
        System.out.println();

        for (Node n : PhyloUtils.internalNodes(t)) {
            System.out.println(n.getNumber());
            System.out.println(n);
            System.out.println(n.getParent());

            List<String> nodes = Lists.newArrayList();
            allChildren(n, nodes);
            System.out.println(nodes);
            System.out.println();
        }

*/


        Set<List<String>> one = Sets.newHashSet();
        Set<List<String>> two = Sets.newHashSet();

        RerootedTreeIterator rti = new RerootedTreeIterator(t);

        for (Tree t1 : rti) {

            System.out.println(t1.toString().replaceAll(":0.0000000", ""));
            for (Node n : PhyloUtils.internalNodes(t1)) {
                List<String> nodes = Lists.newArrayList();
                allChildren(n, nodes);
                two.add(nodes);
                if (n.getParent() != null && n.getParent().isRoot()) System.out.printf(">> ");
                if (n.isRoot()) {
                    for (Node n2 : PhyloUtils.externalNodes(t1)) {
                        if (n2.getParent().isRoot()) {
                            System.out.printf(">> ");
                            System.out.println("[" + n2.getIdentifier().getName() + "]");
                        }
                    }
                    System.out.printf("** ");
                }
                System.out.println(nodes);
            }

            System.out.println();

            System.out.println(one);
            System.out.println(two);
            System.out.println();

            if (one.size() > 0 && two.size() > 0) {
                System.out.printf("Different: %s\n", Sets.difference(one, two));
                System.out.printf("Same: %s\n", Sets.intersection(one, two));
                System.out.println();
            }

            one = two;
            two = Sets.newHashSet();

        }



    }

    private static void allChildren(Node node, List<String> nodes) {

        for (int i = 0; i < node.getChildCount(); i++) {
            if (node.getChild(i).isLeaf()) {
                nodes.add(node.getChild(i).getIdentifier().getName());
            }
            allChildren(node.getChild(i), nodes);
        }

    }

    private static void treewithsequence() {
        GeneticCode.setCode(GeneticCode.STANDARD_CODE);

        Tree t = PhyloUtils.readTree("/Users/atamuri/Downloads/tmp/hiv_vpu_tdg_stuff/vpu_fasttree_codeml_rooted_Hu_nH.nwk.tre");
        Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Downloads/tmp/hiv_vpu_tdg_stuff/vpu_final_codon.phy");


        ArgumentsProcessor ap = new ArgumentsProcessor();
        Table<String, Integer, Byte> codonTable = ap.getCleanedSitesTable(a);


        boolean translate = true;
        int siteOfInterest = 5;

        for (Node n : PhyloUtils.externalNodes(t)) {

            String id = n.getIdentifier().getName();

            byte codon = codonTable.get(id, siteOfInterest);

            String label;

            if (GeneticCode.getInstance().isSenseCodon(codon)) {

                if (translate) {

                    label = GeneticCode.getInstance().getCodonTLA(codon) + " " + GeneticCode.getInstance().getAminoAcidCharByIndex(GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(codon));

                } else {

                    label = GeneticCode.getInstance().getCodonTLA(codon);

                }

            } else {

                label = "-";

            }

            System.out.printf("%s\t%s\t%s\n", codon, GeneticCode.getInstance().getCodonTLA(codon), label);

            n.getIdentifier().setName(label + " " + n.getIdentifier().getName());


        }


        System.out.printf("%s\n", t.toString());


    }

    /*
    static final double[] probMatrix = new double[64 * 64];
    static final double[] tipConditionals = new double[(GeneticCode.CODON_STATES + 1) * GeneticCode.CODON_STATES];


    static byte[] internalNodeChildCount;
    static int[] internalNodeChildNumber;
    static double[] internalNodeChildBranchLength;
    static double[] newConditionals;
    static int[] codons_ = GeneticCode.getInstance().getSenseCodons();
    static int[] childLeafOffset;
    static Table<String, Integer, Byte> table;
static Tree tree;
    static int pos;
    static SubstitutionModel model;*/

   /* private static void fastcalculator() {
        GeneticCode.setCode(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);

        double alignmentsum = 0;



        for (int site : table.columnKeySet()) {
            pos = 0;
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                for (int j = 0; j < tree.getInternalNode(i).getChildCount(); j++) {
                    if (tree.getInternalNode(i).getChild(j).isLeaf()) { // not leaf, internal
                        internalNodeChildNumber[pos] = -1 - tree.getInternalNode(i).getChild(j).getNumber();
                        internalNodeChildBranchLength[pos] = tree.getInternalNode(i).getChild(j).getBranchLength();


                        byte codon = table.column(site).get(tree.getInternalNode(i).getChild(j).getIdentifier().getName());

                        if (codon < 64 && codon >= 0) {
                            childLeafOffset[pos] = codon * 64;
                        } else {
                            childLeafOffset[pos] = 64 * 64;
                        }
                    }

                    pos++;
                }
            }


            Arrays.fill(newConditionals, 1);
            pos = 0;
            double branchProb;
            double[] lowerToUse;
            int lowerOffset, childNumber;
            int i, j, childcount, left, right, leftoffset, rightoffset;


            for (i = 0; i < internalNodeChildCount.length; i++) {

                //System.out.printf("internal node = %s\n", i);


                for (j = 0; j < internalNodeChildCount[i]; j++) {

                    childNumber = internalNodeChildNumber[pos];

                    //  System.out.printf("child %s\n", childNumber);


                    model.getTransitionProbabilities(probMatrix, internalNodeChildBranchLength[pos]);


                    if (childNumber < 0) {
                        //    System.out.printf("it's a leaf\n");
                        // is a leaf
                        //lowerToUse = tipConditionals;
                        lowerOffset = childLeafOffset[pos];

                        for (int f : codons_) {
                            branchProb = 0;
                            for (int t : codons_) {
                                branchProb += tipConditionals[lowerOffset + t] * probMatrix[f * GeneticCode.CODON_STATES + t]; // dot product of conditional and P(t) row
                            }
                            newConditionals[i * GeneticCode.CODON_STATES + f] *= branchProb;
                            // System.out.printf(" %9.6f", newConditionals[i * GeneticCode.CODON_STATES + f]);
                        }


                    } else {
                        //lowerToUse = newConditionals;
                        lowerOffset = childLeafOffset[pos];

                        for (int t : codons_) {
                            branchProb = 0;
                            for (int f : codons_) {
                                branchProb += newConditionals[lowerOffset + f] * probMatrix[t * GeneticCode.CODON_STATES + f];
                            }
                            newConditionals[i * GeneticCode.CODON_STATES + t] *= branchProb;
                            // System.out.printf(" %9.6f", newConditionals[i * GeneticCode.CODON_STATES + f]);
                        }
                    }


                    // System.out.printf("lower to use:");

          *//*      System.out.println();
                for (int iii = 0; iii < 64; iii++) {
                    System.out.printf("%9.6f", lowerToUse[lowerOffset + iii]);
                }
                System.out.println();
*//*


                    //System.out.println();


                    pos++;
                }


            }

        *//*

        System.out.println();
        for (int i = 0; i < 64; i++) {
            System.out.printf("%9.6f", newConditionals[tree.getRoot().getNumber() * GeneticCode.CODON_STATES + i]);
        }
        System.out.println();*//*


            double[] f = model.getCodonFrequencies();

            double sum = 0.0;
            for (i = 0; i < GeneticCode.CODON_STATES; i++)
                sum += newConditionals[tree.getRoot().getNumber() * GeneticCode.CODON_STATES + i] * f[i];

            if (sum < 0) sum = 0;
            sum = Math.log(sum);

            alignmentsum +=sum;
        }


        System.out.printf("%s\n", alignmentsum);


    }*/


    // System.out.printf("%s\n", Ints.join(" ", internalNodeChildNumber));


        /*
        {
                    internalNodeChildNumber[pos] = -1 - tree.getInternalNode(i).getChild(j).getNumber();
                    internalNodeChildBranchLength[pos] = tree.getInternalNode(i).getChild(j).getBranchLength();


                    byte codon = states.get(tree.getInternalNode(i).getChild(j).getIdentifier().getName());

                    if (codon < 64 && codon >= 0) {
                        childLeafOffset[pos] = codon * 64;
                    } else {
                        childLeafOffset[pos] = 64 * 64;
                    }
                }
         */




    private static void testswmutsel() {
        GeneticCode.setCode(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);


        SwMut m = new SwMut(5, 1, new double[]{1, 1, 1, 1});
        m.build();


        Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2013/tdg12/etc/all.but.ND6.AUT.phys");
        Table<String, Integer, Byte> sites = new ArgumentsProcessor().getCleanedSitesTable(a);


        Tree t = PhyloUtils.readTree("/Users/atamuri/Documents/2013/tdg12/etc/all.but.ND6.AUT.FMutSel0.tree");
        MatrixArrayPool.setTreeSize(t.getInternalNodeCount());


        MultiThreadedRunner r = new MultiThreadedRunner(3);
        r.setSites(sites);


        FitnessStore f = new FitnessStore();
        for (int i = 1; i < 3599; i++) {
            f.put(i, new Fitness(new double[20]));
        }

        System.out.printf("%s\n", r.getLogLikelihood(t, m, f, null));

    }

    private static void testpatterns() {
        Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2013/tdg12/etc/all.but.ND6.AUT.phys");

        Table<String, Integer, Byte> sites = new ArgumentsProcessor().getCleanedSitesTable(a);


        Set<Collection<Byte>> set = Sets.newHashSet();

        Set<String> set2 = Sets.newHashSet();

        for (int site : sites.columnKeySet()) {
            String s = "";
            Collection<Byte> b = sites.column(site).values();
            for (byte bb : b) {
                s = s.concat(new String(new byte[]{bb}));
            }
            set2.add(s);
        }


        System.out.printf("%s -> %s\n", sites.columnKeySet().size(), set2.size());


    }

    private static void testS() {

        double[] out = new double[12500000];

        long start = System.currentTimeMillis();

        for (int i = 0; i < out.length; i++) {
            out[i] = getRelativeFixationProbability(Math.random());
        }


        long end = System.currentTimeMillis();


        System.out.printf("%s\n", (end - start) / 1000.0);

        System.out.printf("%s\n", CoreUtils.sum(out));

    }


    private static double getRelativeFixationProbability(double S) {
        // TODO: What's wrong with the simpler:
        return (S == 0) ? 1 : S / (1 - Math.exp(-S));

        //if (S == 0) return 1;
        // return S / (1 - Math.exp(-S));

        /*if (S == 0) return 1;
        if (S < -1e3) return 0;
        if (S > 1e3) return S;
        return S / -Math.expm1(-S);*/
    }

/*

    private static void testfmutsel() {

        GeneticCode.setCode(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);

        System.out.printf("%s\n", Ints.join(",", GeneticCode.getInstance().getSenseCodons()));
*/
/*
        BaseFrequencies d = new BaseFrequencies(new double[]{1, 1, 1, 1});
        System.out.printf("d = %s\n", Doubles.join(" ", d.get()));
        System.out.printf("out = %s\n", Doubles.join(" ", d.getOptimisable()));
        double[] out_ = d.getOptimisable();
        d.setOptimisable(d.getOptimisable());
        System.out.printf("d = %s\n", Doubles.join(" ", d.get()));
        System.exit(0);*//*



        Tree t = PhyloUtils.readTree("/Users/atamuri/Documents/2013/mutselbranch/data/small4s/codeml/tree.txt");
        //Tree t = PhyloUtils.readTree("/Users/atamuri/Documents/2013/tdg12/etc/all.but.ND6.AUT.FMutSel0.tree");

        MatrixArrayPool.treeSize = t.getInternalNodeCount();
        Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2013/mutselbranch/data/small4s/small4s_mtCDNA_500.txt");
        //Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2013/tdg12/etc/all.but.ND6.AUT.phys");



        ArgumentsProcessor p = new ArgumentsProcessor();
        Table<String, Integer, Byte> table = p.getCleanedSitesTable(a);

        table = GuavaUtils.getColumnSubset(table, Ints.asList(CoreUtils.seqi(1, 311)));

        int[] residueCount = new int[20];
        double[] observedCodons = new double[64];
        int observedSum = 0;




        for (int i : table.columnKeySet()) {
            Map<String, Byte> col = table.column(i);

            for (byte b : col.values()) {

                if (GeneticCode.getInstance().isSenseCodon(b)) {
                    observedCodons[b]++;
                    observedSum++;
                }

                if (!GeneticCode.getInstance().isUnknownCodonState(b) &&
                        GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(b) >= 0) {
                    residueCount[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(b)]++;
                }

            }

        }


        for (int i = 0; i < 64; i++) {
            observedCodons[i] /= observedSum;
        }

        System.out.printf("%s\n", Ints.join(" ", residueCount));
        System.out.printf("%s\n", Doubles.join(" ", observedCodons));



        double[] observedBases = new double[]{0.2786667, 0.2895000, 0.3008333, 0.1310000};


        double[] estimatedFitness = new double[20];


        for (int i = 0; i < 20; i++) {

            int[] codons = GeneticCode.getInstance().getCodonIndexFromAminoAcidIndex(i);

            double codonFreqSum = 0;
            double nucFreqSum = 0;

            for (int c : codons) {
                codonFreqSum += observedCodons[c];


                char[] nucs = GeneticCode.getInstance().getNucleotidesFromCodonIndex(c);

                double prod = 1;
                for (char n : nucs) {
                    prod *= observedBases[GeneticCode.getInstance().getNucleotideIndexByChar(n)];
                }

                nucFreqSum += prod;
            }


            estimatedFitness[i] = Math.log(codonFreqSum / nucFreqSum);

            if (Double.NEGATIVE_INFINITY == estimatedFitness[i]) estimatedFitness[i] = -25;

        }





        int fixedResidue = 7;

        double adjust = 0;
        if (estimatedFitness[fixedResidue] > 0) {
            adjust = -estimatedFitness[fixedResidue];
        } else if (estimatedFitness[fixedResidue] < 0) {
            adjust = estimatedFitness[fixedResidue];
        }

        for (int i = 0; i < 20; i++) {
            estimatedFitness[i] += adjust;
        }


        System.out.printf("%s\n", Doubles.join(" ", estimatedFitness));




        //table = GuavaUtils.getColumnSubset(table, Lists.newArrayList(1));

        MultiThreadedRunner r = new MultiThreadedRunner(1);
        r.setSites(table);

        TsTvRatio k = new TsTvRatio(5.0);

        SubstitutionModel model = new FMutSel0(
                k,
                new Omega(0.001),
                new BaseFrequencies(observedBases), // initial pi are observed frequencies
                new Fitness(
                        new int[]{7, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
                        //new double[]{-0.978738, -2.741638, -3.007054, -2.449005, -25.90481, -3.087357, -2.614215, 0.0, -3.078092, -1.279761, -1.310963, -3.866942, -1.893621, -0.793491, -2.537767, -2.17783, -2.577635, -1.070589, -3.024342, -0.95129}));
                        estimatedFitness));
        model.build();

        long start1, end1;
        for (int i = 0; i < 10000; i++) {
            start1 = System.currentTimeMillis();
            model.build();
            r.getLogLikelihood(t, model);
            // double xxx = CoreUtils.sum(r.getLogLikelihoodForTree(t, model).values());
            end1 =  System.currentTimeMillis();
            if (i % 100 == 0) System.out.printf("%s\n", (double) end1 - start1 );
        }


        CodeTimer.printAll();
        System.exit(0);
*/
/*
        TsTvRatio k = new TsTvRatio(7.96527);

        SubstitutionModel model = new FMutSel0(
                k,
                new Omega(0.05327),
                new BaseFrequencies(new double[]{0.17572, 0.28418, 0.45954, 0.08056}),
                new Fitness(
                        new int[]{7, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
                        new double[]{-0.978738, -2.741638, -3.007054, -2.449005, -25.90481, -3.087357, -2.614215, 0.0, -3.078092, -1.279761, -1.310963, -3.866942, -1.893621, -0.793491, -2.537767, -2.17783, -2.577635, -1.070589, -3.024342, -0.95129}));
        // model.build();*//*


        // -3544.1868293400503 - original


*/
/*

zy's pi parameters - one fixed to 1


e.g.

> y <- c(0.17571676, 0.28418062, 0.45954575, 0.08055687)
> y * (1 / y[4])
[1] 2.181276 3.527702 5.704613 1.000000


and backwards:

> x <- c(2.181276, 3.527702, 5.704613)
> x <- c(x, 1)
> x / sum(x)
[1] 0.17571676 0.28418062 0.45954575 0.08055687
ARNDCQEGHILKMFPSTWYV
                                                A          R           N           D          C            Q           E            G          H           I           L           K            M           F            P           S           T           W           Y           omega
7.965275    2.181286    3.527704    5.704601   -0.027447   -1.790347   -2.055764   -1.497714  -24.953519   -2.136067   -1.662924    0.951290   -2.126801   -0.328471   -0.359673   -2.915651   -0.942330    0.157800   -1.586476   -1.226539   -1.626344   -0.119299   -2.073051    0.053266





 *//*



        // 6441, -3559 - one at a time (10 iterations)
        // 5229 -3553 = neldermead
        //



        model.build();

        long start = CodeTimer.start();
        //Map<Integer, Double> out = r.getLogLikelihoodForTree(t, model);

        Pair<Map<Integer, Double>, SubstitutionModel> result = r.optimiseModel(t, model);

        CodeTimer.store("FMutSel0 ", start);


        double sum = 0;

        System.out.printf("sum = %s\n", CoreUtils.sum(result.first.values()));
        CodeTimer.printAll();


        System.out.printf("total evals = %s\n", r.totalevals);

        FMutSel0 outmodel = (FMutSel0) result.second;

        List<swmutsel.model.parameters.Parameter> x = outmodel.getParameters();

        for (swmutsel.model.parameters.Parameter para : x) {
            if (para.getClass() == TsTvRatio.class) {
                System.out.printf("kappa = %s\n", ((TsTvRatio) para).get());
            } else if (para.getClass() == BaseFrequencies.class) {
                System.out.printf("basefreqs = %s\n", Doubles.join(",", ((BaseFrequencies) para).get()));

            } else if (para.getClass() == Omega.class) {
                System.out.printf("omega = %s\n", ((Omega) para).get());

            } else if (para.getClass() == Fitness.class) {
                System.out.printf("fitness = %s\n", Doubles.join(",", ((Fitness) para).get()));

            }
        }

        r.shutdown();


    }
*/


    private static void alignmentTranslate() {
        Alignment a = PhyloUtils.readAlignment("/Users/atamuri/Documents/2014/rg_hiv/hM_cz_vpu/hM_cz_vpu_codon.phy");
        System.out.printf("%s\n", a.getSiteCount());

        DataTranslator dt = new DataTranslator(a);

        Codons codons = new Codons();
        SpecificAminoAcids aminoAcids = new SpecificAminoAcids();

        Alignment b = dt.toAlignment(aminoAcids, 0);


/*
        System.out.printf("%s\n", b.getSiteCount());

        System.out.printf("%s\n", Codons.getTLA(codons.getStateImpl(b.getData(0,0))));
*/


        for (int i = 0; i < b.getSequenceCount(); i++) {
            System.out.printf(">%s\n", b.getIdentifier(i).getName());
            for (int j = 0; j < b.getSiteCount(); j++) {
                System.out.printf("%s", b.getData(i, j));
            }
            System.out.println();
        }

        //dta.toAlignment(Molec, 0);



    }

    private static void singleSubstitutionCodons() {
        int count = 1;
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                if (i == j) continue;

                char[] ni = GeneticCode.getInstance().getNucleotidesFromCodonIndex(i);
                char[] nj = GeneticCode.getInstance().getNucleotidesFromCodonIndex(j);

                int changes = 0;

                for (int k = 0; k < 3; k++) {
                    if (ni[k] != nj[k]) {
                        changes++;
                    }
                }

                if (changes == 1) System.out.printf("%s: %s, %s\n", count++, i, j);
            }
        }
    }


    private static void testJcommander() {
        Options o = new Options();
        JCommander j = new JCommander(o);

        j.setAcceptUnknownOptions(true);
        j.parse("-five", "gah", "-one", "one", "-fitness", "1,2,3,4,5", "6,7,8,9,10", "-four", "four", "-twelve", "twelve", "eleven");


        // first parse only -m and -n -> (i) check name and (ii) run model


        System.out.printf("%s\n", o.toString());
    }

    @Parameters(separators = "=")
    private static class Options {
        @Parameter(names = "-test")
        public boolean test = false;

        @Parameter(names = "-one")
        public String one;

        @Parameter(names = "-two")
        public String two = "two";

        @Parameter(names = "-three")
        public String three;

        @Parameter(names = "-four")
        public String four;

        @Parameter(names = "-five")
        public String five;

        @Parameter(names = "-fitness", variableArity = true)
        public List<String> fitness = Lists.newArrayList();

        @Parameter(description = "other")
        public List<String> other = Lists.newArrayList();

        @Override
        public String toString() {
            return "Options{" +
                    "test = " + test +
                    ", one='" + one + '\'' +
                    ", two='" + two + '\'' +
                    ", three='" + three + '\'' +
                    ", four='" + four + '\'' +
                    ", five='" + five + '\'' +
                    ", \nfitness=" + fitness +
                    ", \nother=" + other +
                    '}';
        }
    }

    public static void testTreeRerooting() {
        // String path = "/Users/atamuri/Documents/2013/mutselbranch/data/small4s/small4s_mtCDNA_500.tree";

        String s1 = "(A:0.1, B:0.2, (C:0.3, D:0.4):0.5);";
        String s2 = "((A:0.1, B:0.2):0.24,(C:0.3, D:0.4):0.26):0.0;";


        PushbackReader pr = new PushbackReader(new StringReader(s2));

        SimpleTree t;
        try {
            t = new ReadTree(pr);
        } catch (Exception e) {
            throw new RuntimeException(e);
        }



        print(t);
        System.out.printf("\n\n\n");

        Tree tt = TreeManipulator.getUnrooted(t);

        print(tt);
        System.out.printf("\n\n\n");

        RerootedTreeIterator rti = new RerootedTreeIterator(tt);

        for (Tree t1 : rti) {
            print(t1);
            System.out.printf("\n\n\n");
        }

    }

    public static void print(Tree t) {
        for (Node n : PhyloUtils.internalNodes(t)) {
            System.out.printf("%s, %s, %s\n", n, n.getChildCount(), n.getNumber());

            for (int i = 0; i < n.getChildCount(); i++) {
                Node child = n.getChild(i);
                System.out.printf("\t%s, %s, %s, %s\n", child.getNumber(), myNumber(t, child), child.getIdentifier().getName(), child.getBranchLength());
            }
        }
    }

    private static int myNumber(Tree t, Node child) {
        if (child.isLeaf()) {
            return child.getNumber();
        } else {
            return t.getExternalNodeCount() + child.getNumber();
        }
    }
}


