package swmutsel.sim;

import com.google.common.base.Charsets;
import com.google.common.collect.*;
import com.google.common.io.Files;
import com.google.common.primitives.Ints;
import swmutsel.Handler;
import swmutsel.cli.SimulatorArguments;
import swmutsel.utils.Functions;

import java.io.*;
import java.util.Arrays;
import java.util.Collection;
import java.util.Map;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class SimulateHandler implements Handler {
    private final SimulatorArguments args;

    public SimulateHandler(SimulatorArguments args) {
        this.args = args;
        validate();
    }

    public void invoke() {
        Simulator s = new Simulator();
        s.setTree(args.tree);
        s.setMutation(args.kappa, args.pi, args.scaling);
        s.setClades(args.cladeModel);
        s.setAminoAcids(args.residues);
        s.setShiftFraction(args.shiftFrac);
        s.setQuiet(args.quiet);
        s.setCachePt(args.cachept);

        if (args.Zn == 0) {
            args.Zn = Double.POSITIVE_INFINITY;
        }

        // If we're simulating a single set of fitnesses, specified using the -fitness option
        if (args.fitness.size() > 0) {
            s.setFdsFitness(args.Z);
            s.initialise(args.sites);

            for (int i = 0; i < args.cladeModel.size(); i++) {
                s.setCladeModel(args.cladeModel.get(i), args.fitness.subList(i * args.residues.length, (i + 1) * args.residues.length));
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
            for (int i = 0; i < args.cladeModel.size(); i++) {
                try {
                    BufferedReader reader = Files.newReader(new File(args.fitnessFiles.get(i)), Charsets.US_ASCII);
                    cladeFitnessReaders.put(args.cladeModel.get(i), reader);
                } catch (FileNotFoundException e) {
                    throw new RuntimeException(e);
                }
            }

            // Collect each simulated site here
            ListMultimap<String, Integer> simSites = ArrayListMultimap.create();

            // Start reading the fitnesses from the fitness file(s)
            int site = 1;
            String line;
            int sitesWithFDS = 0;

            try {
                while ((line = cladeFitnessReaders.get(args.cladeModel.get(0)).readLine()) != null) {
                    if (args.Z > 0 && site <= args.Zn) {
                        s.setFdsFitness(args.Z);
                        sitesWithFDS++;
                    } else {
                        s.setFdsFitness(0);
                    }

                    if (!args.quiet) {
                        System.out.printf("Site %s:\n", site);
                    }

                    // "line" now holds the fitnesses for the first clade - set the clade model for the simulator
                    s.setCladeModel(args.cladeModel.get(0), Lists.transform(Arrays.asList(line.split(" ")), Functions.stringToDouble()));

                    // there may be multiple clades - read a line from the corresponding fitness file and set the clade model
                    for (int i = 1; i < args.cladeModel.size(); i++) {
                        s.setCladeModel(args.cladeModel.get(i), Lists.transform(Arrays.asList(cladeFitnessReaders.get(args.cladeModel.get(i)).readLine().split(" ")), Functions.stringToDouble()));
                    }

                    // finally, simulate the site
                    s.simulate();

                    // we want to save this simulated site in our own set to have an alignment of multiple sites
                    for (Map.Entry<String, int[]> e : s.getSimulatedData().entrySet()) {
                        simSites.put(e.getKey(), e.getValue()[0]); // remember, only one site
                    }

                    site++;
                }
            } catch (IOException e) {
                throw new RuntimeException(e);
            }

            System.out.println("Sites with FDS: " + sitesWithFDS);

            // close all readers
            for (BufferedReader br : cladeFitnessReaders.values()) {
                try {
                    br.close();
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }

            writeOutput(simSites.asMap());
        }
    }

    public void writeOutput(Map<String, Collection<Integer>> sequences) {
        try {
            BufferedWriter out = Files.newWriter(new File(args.identifier), Charsets.US_ASCII);

            // The PHYLIP header: <taxa>   <seq_length>
            out.write(sequences.keySet().size() + "   " + sequences.values().iterator().next().size() * 3);
            out.newLine();

            // sort labels
            String[] labels = sequences.keySet().toArray(new String[0]);
            Arrays.sort(labels);

            // for (Map.Entry<String, Collection<Integer>> e : sequences.entrySet()) {
            for (String label : labels) {
                Collection<Integer> value = sequences.get(label);
                Collection<String> codons = Collections2.transform(value, Functions.codonIndexToTLA());
                out.write(label + "     ");
                for (String c : codons) out.write(c);
                out.newLine();
            }

            out.close();
        } catch (IOException ioe) {
            throw new RuntimeException(ioe);
        }
    }

    public void validate() {
        // Should have specified either -fitness or a file or fitness values using -fitnessfile
        if (args.fitness.size() == 0 && args.fitnessFiles.size() == 0) {
            throw new RuntimeException("You must specify either -fitness or -fitnessfile.");
        }

        // We can't have both -fitness and -fitnessfile specified
        if (args.fitness.size() > 0 && args.fitnessFiles.size() > 0) {
            throw new RuntimeException("You can only use -fitness or -fitnessfile, not both together.");
        }

        // If user has specified fitness using the -fitness parameter
        if (args.fitness.size() > 0) {
            // If we're running the homogeneous model
            if (args.cladeModel.size() == 1) {
                // We expect the same number of fitnesses to characters.
                if (args.fitness.size() != args.residues.length) {
                    throw new RuntimeException(
                            String.format("You have %s fitnesses but %s characters. They should be the same number for each.\n",
                                    args.fitness.size(), args.residues.length));
                }
                // Otherwise, we have a heterogeneous model
            } else {
                // We expect the number of fitnesses to be (number of characters * number of clades)
                if (args.fitness.size() != (args.residues.length * args.cladeModel.size())) {
                    throw new RuntimeException(
                            String.format("You have %s fitnesses for %s clades and %s characters. You should have %s fitnesses (%s for each clade).\n",
                                    args.fitness.size(), args.cladeModel.size(), args.residues.length, (args.residues.length * args.cladeModel.size()), args.residues.length));
                }
            }
        } else {
            // we're using fitness files via -fitnessfile option. Check that we have the right number of files
            if (args.cladeModel.size() != args.fitnessFiles.size()) {
                throw new RuntimeException("You have defined " + args.cladeModel.size() + " clades but have " + args.fitnessFiles.size() + " fitness file(s).");
            }
        }
    }

}
