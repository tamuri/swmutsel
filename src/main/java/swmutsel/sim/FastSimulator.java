package swmutsel.sim;

import com.google.common.primitives.Doubles;
import pal.tree.Node;
import pal.tree.SimpleTree;
import pal.tree.Tree;
import swmutsel.model.SwMut;
import swmutsel.model.SwMutSel;
import swmutsel.model.parameters.Fitness;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.PhyloUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.Random;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 12/11/2014 11:08
 */
public class FastSimulator {
    public static void main(String[] args) throws Exception {
        FastSimulator fs = new FastSimulator();
        fs.run(args);
    }

    private void run(String[] args) throws Exception {
        GeneticCode.setCode(GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE);
        final String fitnessFile = args[0];
        final String outFile = args[1];
        final String treeFile = args[2];
        final Tree tree = PhyloUtils.readTree(treeFile);

        final int sites = 100000;
        final double kappa = 2.0;
        final double[] pi = new double[]{0.25, 0.25, 0.25, 0.25};
        final double scaling = 1.0;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            tree.getExternalNode(i).setSequence(new byte[sites]);
        }

        SwMut mutation = new SwMut(kappa, scaling, pi);
        mutation.build();

        BufferedReader reader = Files.newBufferedReader(new File(fitnessFile).toPath());

        for (int i = 0; i < sites; i++) {

            String line = reader.readLine();
            String[] parts = line.split(" ");
            double[] fitness = new double[20];
            for (int j = 0; j < parts.length; j++) {
                fitness[j] = Double.parseDouble(parts[j]);
            }

            Fitness f = new Fitness(fitness);
            SwMutSel model = new SwMutSel(mutation, f);
            model.build();

            double[] freq = model.getCodonFrequencies();
            byte state = randomState(freq);

            simulate(model, tree.getRoot(), state, i);

            if (i % 500 == 0) System.out.printf("%d ", i);
        }


        // PrintWriter out = new PrintWriter(System.out);
        PrintWriter out = new PrintWriter(new FileWriter(outFile));
        printAlignment(out, tree);
        out.flush();
        out.close();

    }

    private void printAlignment(PrintWriter out, Tree tree) {
        // PrintWriter out = new PrintWriter(System.out);

        out.write("  " + tree.getExternalNodeCount() + " " + tree.getExternalNode(0).getSequence().length * 3 + "\n");

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            Node n = tree.getExternalNode(i);
            out.write(n.getIdentifier().getName() + "   ");
            byte[] sequence = n.getSequence();
            for (int j = 0; j < sequence.length; j++) {
                out.write(GeneticCode.getInstance().getCodonTLA(sequence[j]));
            }
            out.write("\n");
        }
        out.write("\n");
        out.flush();
    }

    private void simulate(SwMutSel model, Node node, byte state, int site) {

        for (int i = 0; i < node.getChildCount(); i++) {

            Node child = node.getChild(i);

            model.getPtCalculator().getTransitionProbabilities(Pt, child.getBranchLength());
            double[] Pij = Pt[state];
            byte newState = randomState(Pij);

            if (child.isLeaf()) {

                child.getSequence()[site] = newState;

            } else {

                simulate(model, child, newState, site);

            }

        }

    }

    private final Random random = new Random(1234567890l);
    private final double[][] Pt = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];

    // Given a probability vector, randomly select a state
    private byte randomState(double[] p) {
        double rnd = random.nextDouble();

        for (byte i = 0; i < p.length; i++) {
            if (rnd < p[i]) return i;
            rnd = rnd - p[i];
        }

        throw new RuntimeException("Couldn't pick random state from [ " + Doubles.join(" ", p) + " ]");
    }
}
