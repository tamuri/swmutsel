package swmutsel.sim;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.google.common.base.Charsets;
import com.google.common.collect.*;
import com.google.common.io.Files;
import com.google.common.primitives.Ints;
import swmutsel.Constants;
import swmutsel.options.CharArrayConverter;
import swmutsel.options.DoubleArrayConverter;
import swmutsel.options.GeneticCodeConverter;
import swmutsel.utils.Functions;
import swmutsel.utils.GeneticCode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class AlignmentSimulator {

    public static void main(String[] args) throws Exception {
        AlignmentSimulator as = new AlignmentSimulator();
        JCommander jc = new JCommander(as);
        jc.setProgramName("java -cp " + Constants.PROGRAM_JAR + "swmutsel.sim.AlignmentSimulator ");

        try {
            //jc.parse("-tree", "/Users/atamuri/Documents/work/tdg12/etc/PB2_FMutSel0.tree.out", "-sites", "10", "-heteroclades", "Av,Hu", "-fitness", "1,2,3,4,5,6,7,8", "-fitness", "1,2,3,4,5,6,7,8", "-characters", "A,R,N,D,C,Q,E,H", "-tau", "1e-2", "-kappa", "7.5", "-pi", "0.25,0.25,0.25,0.25", "-mu", "2.0", "-gc", "standard");
            /*jc.parse("-tree", "/Users/atamuri/Documents/work/tdg12/etc/PB2_FMutSel0.tree.out", "-heteroclades", "Av,Hu",
                    "-fitnessfile", "/Users/atamuri/Documents/work/mitochondria/paper/response/pb2.no.pen/results.nonhomog/fitness.av.sorted.txt",
                    "-fitnessfile", "/Users/atamuri/Documents/work/mitochondria/paper/response/pb2.no.pen/results.nonhomog/fitness.hu.sorted.txt",
                    "-tau", "1.25010000000000e-02", "-kappa", "7.77498", "-pi", "0.22988,0.18954,0.37371,0.20687", "-mu", "2.902626667", "-gc", "standard",
                    "-output", "test.phylip");*/
            /*jc.parse("-tree", "/Users/atamuri/Documents/work/tdg12/etc/PB2_FMutSel0.tree.out", "-heteroclades", "Av,Hu",
                     "-fitness", "0.0, -15.319489016871987, -16.801001904423703, -16.526000568844385, -16.286280170254813, -19.729902687996834, -20.2977063058265, -20.974519424520118, -20.840313189353573, -18.284226997419616, -4.8233322459177455, -18.49536764263169, -19.943683230808453, -17.187213042879524, -18.63223866369517, -1.9348574797996214, -19.541364497586194, -17.540749030863793, -18.16927682552228, -18.698346806595154",
                     "-fitness", "0.0, 2.021057149653928, 0.19484651960663274, -3.064854162224571, -3.3316698123374264, -5.0690529846974295, -4.072399564874824, -6.814501440320555, -8.004932174701963, -9.135355785393472, 20.723108568839045, -15.931545405979719, -14.82133926454409, -16.616690992557334, -17.010476746771904, 20.44479597180679, -19.08063735172478, -15.120089935889876, -16.432255062972438, -17.699794989982887",
                     "-tau", "1e-2", "-kappa", "7.5", "-pi", "0.25,0.25,0.25,0.25", "-mu", "3.2", "-gc", "standard",
                     "-output", "test.phylip",
                     "-sites", "100");*/
            //jc.parse("-tree", "/Users/atamuri/Documents/work/tdg12/etc/PB2_FMutSel0.tree", "-sites", "10", "-fitness", "1,2,3,4,5,6,7,8", "-characters", "A,R,N,D,C,Q,E,H", "-tau", "1e-2", "-kappa", "7.5", "-pi", "0.25,0.25,0.25,0.25", "-mu", "2.0", "-gc", "standard");
            jc.parse(args);
        } catch (ParameterException pe) {
            System.out.printf("Error: %s\n\n", pe.getMessage());
            jc.usage();
            System.exit(0);
        }

        // fitnesses should have same size as characters (homogenous model) or size of characters * heteroclades (heterogeneous model)
        as.validate();
        as.run();
    }

    private void run() throws Exception {
        Simulator s = new Simulator();
        s.setTree(tree);
        s.setMutation(kappa, pi, mu);
        s.setClades(heteroClades);
        s.setAminoAcids(residues);
        s.setShiftFraction(shiftFrac);

        // If we're simulating a single set of fitnesses, specified using the -fitness option
        if (this.fitness.size() > 0) {
            s.initialise(sites);

            for (int i = 0; i < heteroClades.size(); i++) {
                s.setCladeModel(heteroClades.get(i), fitness.subList(i * residues.length, (i + 1) * residues.length));
            }

            s.simulate();

            Map<String, Collection<Integer>> simSites = Maps.newHashMap();
            for (Map.Entry<String, int[]> e : s.getSimulatedData().entrySet()) {
                simSites.put(e.getKey(), Ints.asList(e.getValue()));
            }

            writeOutput(simSites);
        } else {
            // We're reading fitnesses for each site from a file. We simulate each site once.
            s.initialise(1);

            // Fitnesses will be read from file(s). Hold each file reader by the clade key.
            Map<String, BufferedReader> cladeFitnessReaders = Maps.newHashMap();
            for (int i = 0; i < heteroClades.size(); i++) {
                cladeFitnessReaders.put(heteroClades.get(i), Files.newReader(new File(this.fitnessFiles.get(i)), Charsets.US_ASCII));
            }

            // Collect each simulated site here
            ListMultimap<String, Integer> simSites = ArrayListMultimap.create();

            // Start reading the fitnesses from the fitness file(s)
            int site = 1;
            String line;
            while ((line = cladeFitnessReaders.get(heteroClades.get(0)).readLine()) != null) {
                System.out.printf("Site %s:\n", site++);
                // "line" now holds the fitnesses for the first clade - set the clade model for the simulator
                s.setCladeModel(heteroClades.get(0), Lists.transform(Arrays.asList(line.split(" ")), Functions.stringToDouble()));

                // there may be multiple clades - read a line from the corresponding fitness file and set the clade model
                for (int i = 1; i < heteroClades.size(); i++) {
                    s.setCladeModel(heteroClades.get(i), Lists.transform(Arrays.asList(cladeFitnessReaders.get(heteroClades.get(i)).readLine().split(" ")), Functions.stringToDouble()));
                }

                // finally, simulate the site
                s.simulate();

                // we want to save this simulated site in our own set to have an alignment of multiple sites
                for (Map.Entry<String, int[]> e : s.getSimulatedData().entrySet()) {
                    simSites.put(e.getKey(), e.getValue()[0]); // remember, only one site
                }

            }

            // close all readers
            for (BufferedReader br : cladeFitnessReaders.values()) {
                br.close();
            }

            writeOutput(simSites.asMap());
        }
    }

    public void writeOutput(Map<String, Collection<Integer>> sequences) throws Exception {
        BufferedWriter out = Files.newWriter(new File(outputFile), Charsets.US_ASCII);

        // The PHYLIP header: <taxa>   <seq_length>
        out.write(sequences.keySet().size() + "   " + sequences.values().iterator().next().size() * 3);
        out.newLine();

        for (Map.Entry<String, Collection<Integer>> e : sequences.entrySet()) {
            Collection<String> codons = Collections2.transform(e.getValue(), Functions.codonIndexToTLA());
            out.write(e.getKey() + "     ");
            for (String c : codons) out.write(c);
            out.newLine();
        }

        out.close();
    }

    public void validate() {
        // Should have specified either -fitness or a file or fitness values using -fitnessfile
        if (this.fitness.size() == 0 && this.fitnessFiles.size() == 0) {
            throw new RuntimeException("You must specify either -fitness or -fitnessfile.");
        }

        // We can't have both -fitness and -fitnessfile specified
        if (this.fitness.size() > 0 && this.fitnessFiles.size() > 0) {
            throw new RuntimeException("You can only use -fitness or -fitnessfile, not both together.");
        }

        // If user has specified fitness using the -fitness parameter
        if (this.fitness.size() > 0) {
            // If we're running the homogeneous model
            if (this.heteroClades.size() == 1) {
                // We expect the same number of fitnesses to characters.
                if (this.fitness.size() != this.residues.length) {
                    throw new RuntimeException(
                            String.format("You have %s fitnesses but %s characters. They should be the same number for each.\n",
                                    this.fitness.size(), this.residues.length));
                }
                // Otherwise, we have a heterogeneous model
            } else {
                // We expect the number of fitnesses to be (number of characters * number of clades)
                if (this.fitness.size() != (this.residues.length * this.heteroClades.size())) {
                    throw new RuntimeException(
                            String.format("You have %s fitnesses for %s clades and %s characters. You should have %s fitnesses (%s for each clade).\n",
                                    this.fitness.size(), this.heteroClades.size(), this.residues.length, (this.residues.length * this.heteroClades.size()), this.residues.length));
                }
            }
        } else {
            // we're using fitness files via -fitnessfile option. Check that we have the right number of files
            if (this.heteroClades.size() != this.fitnessFiles.size()) {
                throw new RuntimeException("You have defined " + this.heteroClades.size() + " clades but have " + this.fitnessFiles.size() + " fitness file(s).");
            }
        }
    }

    // COMMAND-LINE ARGUMENTS (for JCommander)
    @Parameter(names = "-tree", description = "Tree file in Newick format", required = true)
    public String tree;

    @Parameter(names = "-sites", description = "Number of times to simulate each set of fitnesses.", required = false)
    public int sites = 1;

    // -fitness option can be specified multiple times for heterogeneous models
    @Parameter(names = "-fitness", description = "Comma-separated fitness coefficients.", required = false)
    public List<Double> fitness = Lists.newArrayList();

    // -fitnessfile option can be specified multiple times for heterogenous models
    @Parameter(names = "-fitnessfile", description = "A file containing space-separated fitness coefficients, one row per site.", required = false)
    public List<String> fitnessFiles = Lists.newArrayList();

    // The fitnesses in -fitness can be specified in any order. Default is canonical amino acid order
    @Parameter(names = "-characters", description = "Comma-separated amino acids (matching -fitness order).", converter = CharArrayConverter.class, required = false)
    public char[] residues = GeneticCode.AMINO_ACIDS;

    // You must specify this for heterogeneous models
    @Parameter(names = "-heteroclades", description = "If simulating heterogeneous model, specify the clade labels.", required = false)
    public List<String> heteroClades = Lists.newArrayList("ALL");


    @Parameter(names = "-output", description = "The name of the file to save the simulated alignment.", required = true)
    public String outputFile;

    @Parameter(names = "-shiftfrac", description = "How far down the branch to switch non-homogeneous models. 0 = start of branch, 1 = end of branch, x = percentage of branch")
    public double shiftFrac = 0.5;


    @Parameter(names = "-kappa", description = "Transition/transversion bias.", required = true)
    public double kappa;

    @Parameter(names = "-pi", description = "Comma-separated base nucleotide frequencies (T,C,A,G).", converter = DoubleArrayConverter.class, required = true)
    public double[] pi;

    @Parameter(names = "-mu", description = "Branch/rate scaling factor.", required = true)
    public double mu;

    @Parameter(names = "-gc", description = "The genetic code translation to use (standard or vertebrate_mit).", required = true, converter = GeneticCodeConverter.class)
    public GeneticCode geneticCode;
}
