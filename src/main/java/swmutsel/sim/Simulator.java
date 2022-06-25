package swmutsel.sim;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.Constants;
import swmutsel.model.SwMut;
import swmutsel.model.SwMutSel;
import swmutsel.model.parameters.Fitness;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.PhyloUtils;

import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class Simulator {

    private Tree parsedTree;
    private int sites;
    private Map<Identifier, int[]> seqout;
    private SwMut globals;
    private List<String> heteroClades;
    private List<Integer> aminoAcids = Lists.newArrayList();
    private Map<String, SwMutSel> cladeModels = Maps.newHashMap();
    private double[][] Pt = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
    private double shiftFraction;

    public void initialise(int sites) {
        this.sites = sites;

        // Initialise sequence store
        seqout = Maps.newHashMap();
        seqout.put(parsedTree.getRoot().getIdentifier(), new int[sites]);
        for (int i = 0; i < parsedTree.getExternalNodeCount(); i++) {
            int[] store = new int[sites];
            seqout.put(parsedTree.getExternalNode(i).getIdentifier(), store);
        }
        for (int i = 0; i < parsedTree.getInternalNodeCount(); i++) {
            int[] store = new int[sites];
            seqout.put(parsedTree.getInternalNode(i).getIdentifier(), store);
        }
    }

    public void setTree(Tree tree) {
        parsedTree = tree;
    }

    public void setMutation(double kappa, double[] pi, double mu) {
        globals = new SwMut(kappa, mu, pi);
        globals.build();
    }

    public void simulate() {
        Node root = parsedTree.getRoot();

        // Root all uses the first model in the models list
        for (int i = 0; i < sites; i++) {
            seqout.get(root.getIdentifier())[i] = selectRandomCharacter(this.cladeModels.get(this.heteroClades.get(0)).getCodonFrequencies());
        }

        downTree(root);
    }

    private void downTree(Node parent) {
        for (int i = 0; i < parent.getChildCount(); i++) {
            Node child = parent.getChild(i);

            String model;

            // if we're using a homogeneous model
            if (this.heteroClades.size() == 1) {
                model = "ALL";

                this.cladeModels.get(model).getPtCalculator().getTransitionProbabilities(Pt, child.getBranchLength());

                for (int j = 0; j < this.sites; j++) {
                    int row = seqout.get(parent.getIdentifier())[j];
                    seqout.get(child.getIdentifier())[j] = selectRandomCharacter(Pt[row]);
                }

            } else {
                // this is a heterogeneous model - split the branch in two
                String modelA = (parent.getIdentifier().getName().length() == 0) ? this.heteroClades.get(0) : parent.getIdentifier().getName().substring(0, 2);
                String modelB = (child.getIdentifier().getName().length() == 0) ? this.heteroClades.get(0) : child.getIdentifier().getName().substring(0, 2);

                int[] switchPoint = new int[this.sites];

                // modelA to midpoint
                if (shiftFraction == 0) {
                    for (int j = 0; j < this.sites; j++) {
                        switchPoint[j] = seqout.get(parent.getIdentifier())[j];
                    }
                } else {
                    this.cladeModels.get(modelA).getPtCalculator().getTransitionProbabilities(Pt, child.getBranchLength() * shiftFraction);

                    for (int j = 0; j < this.sites; j++) {
                        int row = seqout.get(parent.getIdentifier())[j];
                        switchPoint[j] = selectRandomCharacter(Pt[row]);
                    }
                }

                // midpoint to modelB
                if (shiftFraction == 1) {
                    for (int j = 0; j < this.sites; j++) {
                        seqout.get(child.getIdentifier())[j] = switchPoint[j];
                    }
                } else {
                    this.cladeModels.get(modelB).getPtCalculator().getTransitionProbabilities(Pt, child.getBranchLength() * (1 - shiftFraction));

                    for (int j = 0; j < this.sites; j++) {
                        int row = switchPoint[j];
                        seqout.get(child.getIdentifier())[j] = selectRandomCharacter(Pt[row]);
                    }
                }
            }

            downTree(child);
        }

    }

    private int selectRandomCharacter(double[] array) {
        double rnd = Math.random();
        //Random random = new Random(123456789);
        //random.nextDouble();
        int character = -1;
        double sum = 0.0;
        for (int i = 0; i < array.length; i++) {
            sum += array[i];
            if (sum > rnd) {
                character = i;
                break;
            }
        }
        return character;
    }

    public void setClades(List<String> heteroClades) {
        this.heteroClades = heteroClades;

        for (String model : heteroClades) {
            this.cladeModels.put(model, null);
        }
    }

    public void setAminoAcids(char[] residues) {
        for (char c : residues) {
            aminoAcids.add(GeneticCode.getInstance().getAminoAcidIndexByChar(c));
        }
    }

    public void setCladeModel(String model, List<Double> fitness) {
        System.out.printf("%s has fitness %s for residues %s.\n", model, fitness, aminoAcids);


        double[] f = new double[GeneticCode.AMINO_ACID_STATES];
        Arrays.fill(f, -Constants.FITNESS_BOUND);


        for (int i = 0; i < aminoAcids.size(); i++) {
            int residue = aminoAcids.get(i);
            f[residue] = fitness.get(i);
        }




        SwMutSel codonModel = new SwMutSel(this.globals, new Fitness(f));
        codonModel.build();

        System.out.printf("Amino acid frequencies are:\n%s\n", CoreUtils.join("%.7f", " ", PhyloUtils.getAminoAcidFrequencies(codonModel.getCodonFrequencies())));
        this.cladeModels.put(model, codonModel);
    }

    public Map<String, int[]> getSimulatedData() {
        Map<String, int[]> data = Maps.newHashMap();
        for (int i = 0; i < parsedTree.getExternalNodeCount(); i++) {
            Identifier n = parsedTree.getExternalNode(i).getIdentifier();
            data.put(n.getName(), seqout.get(n));
        }
        return data;
    }

    public void printSimulatedData() {
        Map<String, int[]> data = getSimulatedData();

        for (Map.Entry<String, int[]> e : data.entrySet()) {
            System.out.printf("%s      ", e.getKey());
            for (int i = 0; i < sites; i++) {
                System.out.printf("%s", GeneticCode.getInstance().getCodonTLA(data.get(e.getKey())[i]));
            }
            System.out.println();
        }
    }

    public void setShiftFraction(double shiftFraction) {
        if (shiftFraction < 0 || shiftFraction > 1) {
            throw new RuntimeException("-shiftfrac should be between 0 and 1");
        }
        this.shiftFraction = shiftFraction;
    }
}
