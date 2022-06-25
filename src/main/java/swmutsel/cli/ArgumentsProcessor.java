package swmutsel.cli;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import com.google.common.collect.*;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Ints;
import pal.alignment.Alignment;
import pal.alignment.AlignmentBuilder;
import pal.alignment.DataTranslator;
import pal.alignment.SitePattern;
import pal.datatype.Codons;
import pal.tree.Node;
import pal.tree.ReadTree;
import pal.tree.Tree;
import swmutsel.Constants;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.GuavaUtils;
import swmutsel.utils.PhyloUtils;

import java.io.*;
import java.nio.charset.Charset;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 18/10/2013 20:32
 */
public class ArgumentsProcessor {
    private ArgumentsProcessor() { }

    public static Arguments parse(final String[] args) {
        String[] merged = mergeFileArguments(args);

        SharedArguments shared = new SharedArguments();
        runJCommander(merged, shared, true);

        if (shared.showHelp) showHelpAndExit();

        Arguments opt;

        if (shared.fdr) {
            opt = new FDSArguments();
        } else if (shared.fmutsel0) {
            opt = new FMutSel0Arguments();
        } else if (shared.simulate) {
            // simulate swmutsel options
            opt = new SimulatorArguments();
        } else {
            // swmutsel options
            opt = new SwMutSelArguments();
        }

        runJCommander(merged, opt, false);

        try {
            opt.initialise();
        } catch (ParameterException e) {
            error(e, opt);
        }

        return opt;
    }

    private static void runJCommander(String[] arguments, Arguments object, boolean acceptUnknown) {
        try {
            JCommander jc = new JCommander(object);
            jc.setAcceptUnknownOptions(acceptUnknown);
            jc.parse(arguments);
        } catch (ParameterException e) {
            if (object instanceof SharedArguments && ((SharedArguments) object).showHelp) {
                showHelpAndExit();
            } else {
                error(e, object);
            }
        }
    }

    private static String getAllCommandsSynopsis() {
        InputStream is;
        String s = "";
        for (String helpFile : new String[]{Constants.HELP_COMMAND_SWMUTSEL, Constants.HELP_COMMAND_SIMULATE, Constants.HELP_COMMAND_FMUTSEL0}) {
            is = ArgumentsProcessor.class.getResourceAsStream(Constants.RESOURCE_PATH + helpFile);
            s += CoreUtils.convertStreamToString(is);
            s += "\n";
        }
        return s;
    }

    private static void error(ParameterException e, Arguments object) {
        String s = "";

        if (object instanceof SimulatorArguments || object instanceof FMutSel0Arguments) {
            InputStream is = ArgumentsProcessor.class.getResourceAsStream(Constants.RESOURCE_PATH + object.getCommandSynopsisPath());
            s += CoreUtils.convertStreamToString(is);
        } else {
            InputStream is = ArgumentsProcessor.class.getResourceAsStream(Constants.RESOURCE_PATH + Constants.HELP_COMMAND_SWMUTSEL);
            s += CoreUtils.convertStreamToString(is);
        }

        // All non-JCommander errors begin with 'ERROR:'
        if (!e.getMessage().startsWith("ERROR:")) {
            // show synopsis of command
            System.out.println("SYNOPSIS");
            System.out.println();
            System.out.println(s);
            System.out.println();
        }

        // show the specific error
        s = e.getMessage();
        if (s.contains("The following args are required")) {
            s = s.replaceAll(", ", "/");
            s = s.replaceAll("([A-z])(\\s)([-])", "$1  $3");
        }
        System.out.printf("%s\n", s);
        System.out.println();

        System.exit(0);
    }

    private static String[] mergeFileArguments(final String[] rawArguments) {

        List<String> arguments = Lists.newArrayList();
        
        for (String arg : rawArguments) {
            // If we've been given an input file
            if (arg.startsWith("@")) {
                String filename = arg.substring(1);
                try {
                    List<String> nonCommentLines = Files.readLines(new File(filename), Charset.defaultCharset(), new LineProcessor<List<String>>() {
                        List<String> lines = Lists.newArrayList();

                        @Override
                        public boolean processLine(String line) throws IOException {
                            // This is the end of the file
                            if (line.startsWith("//")) return false;
                            // Ignore comments and empty lines, ignore everything after a comment
                            if (!line.startsWith("#") && line.length() > 0) lines.add(line.split("#")[0]);
                            return true;
                        }

                        @Override
                        public List<String> getResult() {
                            return lines;
                        }
                    });

                    arguments.addAll(nonCommentLines);

                } catch (IOException e) {
                    CoreUtils.msg("ERROR: Could not read file '%s'\n", arg.substring(1));
                }

            } else {
                arguments.add(arg);
            }
        }

        return arguments.toArray(new String[arguments.size()]);
    }

    private static void showHelpAndExit() {
        InputStream is = ArgumentsProcessor.class.getClass().getResourceAsStream(Constants.RESOURCE_PATH + Constants.HELP_DESCRIPTION);
        String s = CoreUtils.convertStreamToString(is);
        System.out.println();
        System.out.println(s);
        System.out.println();

        s = getAllCommandsSynopsis();
        System.out.println("SYNOPSIS");
        System.out.println();
        System.out.println(s);
        System.out.println();

        is = ArgumentsProcessor.class.getClass().getResourceAsStream(Constants.RESOURCE_PATH + Constants.HELP_FULL);
        s = CoreUtils.convertStreamToString(is);
        System.out.println(s);
        System.out.println();

        System.exit(0);
    }

    public static FitnessStore loadFitnessStore(final List<String> list, final int required, boolean withIndex, final FitnessOrderingCallback orderingCallback) {
        int vectorSize = GeneticCode.AMINO_ACID_STATES;
        if (withIndex) vectorSize++; // 1 'index' and 20 fitnesses

        FitnessStore fitnesses = new FitnessStore();
        if (list == null) {
            for (int i = 1; i <= required; i++) {
                fitnesses.put(i, new Fitness(orderingCallback.getOrderingForOptimiser(i)));
            }
        } else {

            if (list.size() / vectorSize != required) {
                String msg = String.format("ERROR: %s required but %s sets of fitness values given.\n", required, list.size() / vectorSize);
                CoreUtils.msg(msg);
                throw new ParameterException(msg);
            }

            Iterator<String> items = list.iterator();

            int site = 0;

            while (items.hasNext()) {
                if (withIndex)
                    site = Integer.parseInt(items.next());
                else site++;

                double[] f = new double[GeneticCode.AMINO_ACID_STATES];
                for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                    f[i] = Double.parseDouble(items.next());
                }

                Fitness F = new Fitness(orderingCallback.getOrderingForOptimiser(site), f);

                fitnesses.put(site, F);
            }

            CoreUtils.msg("Parsed %s sets of fitness values.\n", fitnesses.size());

            list.clear();
        }

        return fitnesses;
    }

    public static Tree loadTree(final String stringOrFile) {
        Tree tree;
        if (new File(stringOrFile).exists()) {
            CoreUtils.msg("Loading tree file: %s\n", stringOrFile);
            tree = PhyloUtils.readTree(stringOrFile);
        } else {
            try {
                CoreUtils.msg("Loading tree string: %s...\n", stringOrFile.substring(0, stringOrFile.length() > 70 ? 70 : stringOrFile.length()));
                tree = new ReadTree(new PushbackReader(new StringReader(stringOrFile)));
            } catch (Exception e) {
                String msg = "ERROR: Tree '" + stringOrFile + "' is not a file or tree string (" + CoreUtils.getRootCause(e).getMessage()  + ").\n";
                CoreUtils.msg(msg);
                throw new ParameterException(msg);
            }
        }
        CoreUtils.msg("Tree has %s tips.\n", tree.getExternalNodeCount());

        // if all branch lengths are 0, set to initial branch length size
        boolean allZero = true;
        for (Node n : PhyloUtils.externalNodes(tree)) {
            if (n.getParent() != null) if (n.getBranchLength() > 0) allZero = false;
        }
        for (Node n : PhyloUtils.internalNodes(tree)) {
            if (n.getParent() != null) if (n.getBranchLength() > 0) allZero = false;
        }

        if (allZero) {
            for (Node n : PhyloUtils.externalNodes(tree)) {
                if (n.getParent() != null) n.setBranchLength(Constants.INITIAL_BRANCH_LENGTH);
            }
            for (Node n : PhyloUtils.internalNodes(tree)) {
                if (n.getParent() != null) n.setBranchLength(Constants.INITIAL_BRANCH_LENGTH);
            }
        } else {
            // some branch lengths might be 0 which is < minimum branch length
            for (Node n : PhyloUtils.externalNodes(tree)) {
                if (n.getParent() != null && n.getBranchLength() == 0) n.setBranchLength(Constants.MIN_BRANCH_LENGTH);
            }
            for (Node n : PhyloUtils.internalNodes(tree)) {
                if (n.getParent() != null && n.getBranchLength() == 0) n.setBranchLength(Constants.MIN_BRANCH_LENGTH);
            }
        }

        return tree;
    }

    public static Alignment loadPalAlignment(final String path) {
        if (path == null) {
            String msg = "ERROR: -sequences option required.\n";
            CoreUtils.msg(msg);
            throw new ParameterException(msg);
        }

        CoreUtils.msg("Loading alignment file: %s\n", path);

        Alignment nucleotideAlignment;
        try {
            nucleotideAlignment = PhyloUtils.readAlignment(path);
        } catch (Exception e) {
            String msg = "ERROR: Could not read alignment: " + CoreUtils.getRootCause(e).getMessage() + "\n";
            CoreUtils.msg(msg);
            throw new ParameterException(msg);
        }

        // Translate nucleotide alignment to codon alignment
        DataTranslator dataTranslator = new DataTranslator(nucleotideAlignment);
        Alignment codonAlignment = dataTranslator.toAlignment(new Codons(), 0);

        // Replace any unknown triplets and stop codons with gaps, build a clean alignment
        AlignmentBuilder builder = new AlignmentBuilder(codonAlignment.getSequenceCount());

        // For each sequence
        for (int i = 0; i < codonAlignment.getSequenceCount(); i++) {

            // Collect the sequence states
            int[] states = new int[codonAlignment.getSiteCount()];

            // For each site
            for (int j = 0; j < codonAlignment.getSiteCount(); j++) {

                // Get the PAL codon character and codon index
                char palCodonChar = codonAlignment.getData(i, j);
                int palCodonIndex = Codons.DEFAULT_INSTANCE.getState(palCodonChar);

                // Get the three-letter nucleotide
                char[] nucleotides = Codons.getNucleotidesFromCodonIndex(palCodonIndex);

                // Get swMutSel codon index and check whether its a sense codon (handles the genetic code table)
                int codonIndex = GeneticCode.getInstance().getCodonIndexFromNucleotides(nucleotides);

                if (GeneticCode.getInstance().isSenseCodon(codonIndex)) {
                    states[j] = palCodonIndex;
                } else {
                    states[j] = Codons.SUGGESTED_GAP_STATE;
                }
            }
            builder.addSequence(states, codonAlignment.getIdentifier(i).getName());
        }

        Alignment cleanedCodonAlignment = builder.generateAlignment(new Codons());

        return cleanedCodonAlignment;
    }

    public static Table<String, Integer, Byte> getAlignmentTable(Alignment alignment) {

        // Prepare PAL nucleotide alignment as codon sites, replace unknown and STOP codons with gaps
        List<String> sequences = Lists.newArrayList();
        for (int i = 0; i < alignment.getSequenceCount(); i++)
            sequences.add(alignment.getIdentifier(i).getName());

        List<Integer> sites = Ints.asList(CoreUtils.seq(1, alignment.getSiteCount()));

        ArrayTable<String, Integer, Byte> cleanedSites = ArrayTable.create(sequences, sites);

        for (int i = 0; i < alignment.getSiteCount(); i++) {

            for (int j = 0; j < alignment.getSequenceCount(); j++) {

                // Get the PAL codon character and codon index
                char palCodonChar = alignment.getData(j, i);
                int palCodonIndex = Codons.DEFAULT_INSTANCE.getState(palCodonChar);

                // Get the three-letter nucleotide
                char[] nucleotides = Codons.getNucleotidesFromCodonIndex(palCodonIndex);

                int codonIndex = GeneticCode.getInstance().getCodonIndexFromNucleotides(nucleotides);

                cleanedSites.put(alignment.getIdentifier(j).getName(), i + 1, (byte) codonIndex);

            }

        }



        CoreUtils.msg("Alignment has %s sequences and %s codon sites.\n",
                cleanedSites.rowKeySet().size(),
                cleanedSites.columnKeySet().size());
        return cleanedSites;
    }

    public static Table<String,Integer,Byte> findPatterns(
            final boolean noPatterns,
            final Table<String, Integer, Byte> sites,
            Alignment alignment,
            final Map<Integer, Integer> patternSiteMap,
            final Map<Integer, Integer> patternWeights) {

        patternWeights.clear();
        patternSiteMap.clear();

        if (noPatterns) {
            CoreUtils.msg("Not using site patterns.\n");

            // by default, each site has weight = 1
            for (int i : sites.columnKeySet()) patternWeights.put(i, 1);

            // by default each site's pattern is defined by itself
            for (int i : sites.columnKeySet()) patternSiteMap.put(i, i);

            return sites;
        }

        SitePattern patterns = new SitePattern(alignment);

        for (int i = 0; i < patterns.alias.length; i++) {

            int site = i + 1;
            int pattern = patterns.alias[i];

            int firstSiteWithPattern = Ints.indexOf(patterns.alias, pattern) + 1;

            if (!patternWeights.containsKey(firstSiteWithPattern)) {

                patternWeights.put(firstSiteWithPattern, 1);
                patternSiteMap.put(site, site);

            } else {

                patternWeights.put(firstSiteWithPattern, patternWeights.get(firstSiteWithPattern) + 1);
                patternSiteMap.put(site, firstSiteWithPattern);

            }
        }

        Collection<Integer> distinctSites = Lists.newArrayList();
        for (Map.Entry<Integer, Integer> e : patternWeights.entrySet()) {
            if (e.getValue() > 0) distinctSites.add(e.getKey());
        }

        CoreUtils.msg("%s site patterns.\n", distinctSites.size());
        return GuavaUtils.getColumnSubset(sites, distinctSites);
    }

    public static Table<String, Integer, Byte> loadSiteRange(final RangeSet<Integer> siteRange, final Table<String, Integer, Byte> sites, final FitnessStore fitnesses) {
        List<Integer> reduced = Lists.newArrayList();

        for (Range<Integer> r : siteRange.asRanges()) {
            for (int s : ContiguousSet.create(r, DiscreteDomain.integers()).asList()) {
                reduced.add(s);
            }
        }

        CoreUtils.msg("Site range specified %s -> %s site(s) to analyse.\n", siteRange.toString(), reduced.size());

        Iterator<Map.Entry<Integer, Fitness>> iterator = fitnesses.entrySet().iterator();
        while (iterator.hasNext()) {
            Map.Entry<Integer, Fitness> entry = iterator.next();
            if (!reduced.contains(entry.getKey())) {
                iterator.remove();
            }
        }

        return GuavaUtils.getColumnSubset(sites, reduced);
    }

    public static interface FitnessOrderingCallback {
        public int[] getOrderingForOptimiser(int set);
    }
}
