package swmutsel;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.ParameterException;
import com.google.common.collect.*;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Ints;
import pal.alignment.Alignment;
import pal.tree.ReadTree;
import pal.tree.Tree;
import pal.tree.TreeManipulator;
import swmutsel.model.parameters.Fitness;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.GuavaUtils;
import swmutsel.utils.PhyloUtils;

import java.io.*;
import java.nio.charset.Charset;
import java.util.*;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 18/10/2013 20:32
 */
public class ArgumentsProcessor {
    private String identifier;

    public static Arguments parse(String[] args) {
        ArgumentsProcessor ap = new ArgumentsProcessor();

        String[] preProcessedArgs = ap.preProcess(args);

        // writeInputArguments()

        Arguments arguments = ap.parseCommandLine(preProcessedArgs);

        ap.postProcess(arguments);

        if (!ap.valid(arguments)) {
            CoreUtils.msg("ERROR: Arguments not valid. See preceding messages.");
            System.exit(0);
        }

        return arguments;
    }

    public ArgumentsProcessor() {
        // this.identifier = new SimpleDateFormat("yyMMddHHmmss").format(new Date());
    }

    private String[] preProcess(String[] originalArguments) {

        List<String> arguments = Lists.newArrayList();

        // If we've been given an input file
        for (String arg : originalArguments) {
            if (arg.startsWith("@")) {
                String filename = arg.substring(1);

                try {
                    List<String> lines = Files.readLines(new File(filename), Charset.defaultCharset(), new LineProcessor<List<String>>() {
                        List<String> lines = Lists.newArrayList();

                        @Override
                        public boolean processLine(String line) throws IOException {

                            if (line.startsWith("//")) return false; // This is the end of the file
                            if (!line.startsWith("#") && line.length() > 0) lines.add(line); // Ignore comments and empty lines
                            return true;
                        }

                        @Override
                        public List<String> getResult() {
                            return lines;
                        }
                    });

                    arguments.addAll(lines);

                } catch (IOException e) {
                    // TODO: handle file not found!
                    throw new RuntimeException(e);
                }


            } else {
                arguments.add(arg);
            }
        }

        return arguments.toArray(new String[arguments.size()]);
    }

    private Arguments parseCommandLine(String[] args) {
        Arguments arguments = new Arguments();
        JCommander jc = new JCommander(arguments);
        try {
            jc.parse(args);
            if (arguments.showHelp) showHelpAndExit();
        } catch (Exception e) {
            if (arguments.showHelp) {
                showHelpAndExit();
            } else {
                InputStream is = this.getClass().getResourceAsStream("/resources/command.txt");
                String s = CoreUtils.convertStreamToString(is);

                System.out.println();
                System.out.println(s);
                System.out.println();

                s = e.getMessage();
                s = s.replaceAll(", ", "/");
                s = s.replaceAll("([A-z])(\\s)([-])", "$1  $3");

                System.out.printf("%s\n", s);
                System.out.println();
            }
            System.exit(0);
        }

        CoreUtils.msg("Working directory: %s\n", System.getProperty("user.dir"));
        return arguments;
    }
    
    private void showHelpAndExit() {
        InputStream is = this.getClass().getResourceAsStream("/resources/command.txt");
        String s = CoreUtils.convertStreamToString(is);

        System.out.println();
        System.out.println(s);
        System.out.println();

        is = this.getClass().getResourceAsStream("/resources/help.txt");
        s = CoreUtils.convertStreamToString(is);

        System.out.println(s);
        System.out.println();
    }

    private void postProcess(Arguments args) {
        // Identifier
        // args.identifier = this.identifier;

        // Read tree
        Tree tree;
        if (new File(args.treePath).exists()) {
            CoreUtils.msg("Loading tree file: %s\n", args.treePath);
            tree = PhyloUtils.readTree(args.treePath);
        } else {
            try {
                CoreUtils.msg("Loading tree string: %s...\n", args.treePath.substring(0, 70));
                tree = new ReadTree(new PushbackReader(new StringReader(args.treePath)));
            } catch (Exception e) {
                throw new ParameterException("Value " + args.treePath + " is not a path or tree string");
            }
        }

        args.initialTree = tree;
        args.tree = TreeManipulator.getUnrooted(tree);

        // Read alignment
        Alignment alignment = PhyloUtils.readAlignment(args._alignmentPath);
        CoreUtils.msg("Alignment has %s sequences and %s codon sites.\n",
                alignment.getSequenceCount(),
                alignment.getSiteCount() / 3);


        args.sites = getCleanedSitesTable(alignment);


        // Raw -_fitness as single List<String> should be parsed and put into FitnessStore
        if (args._fitness != null && args._fitness.size() > 0) {

            Iterator<String> items = args._fitness.iterator();

            while (items.hasNext()) {

                int site = Integer.parseInt(items.next());

                double[] f = new double[GeneticCode.AMINO_ACID_STATES];
                for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                    f[i] = Double.parseDouble(items.next());
                }

                args.fitnesses.set(site, new Fitness(Ints.toArray(PhyloUtils.getOrderedAminoAcids(args.sites.column(site).values())), f));
            }

            CoreUtils.msg("Parsed fitness vectors for %s sites.\n", args.fitnesses.size());

            args._fitness.clear();
            args._fitness = null;
        } else {
            for (int site : args.sites.columnKeySet()) {
                args.fitnesses.set(site, new Fitness(Ints.toArray(PhyloUtils.getOrderedAminoAcids(args.sites.column(site).values())), new double[GeneticCode.AMINO_ACID_STATES]));
            }
        }


        // Trim the sites and fitnesses to only those we need
        if (args.siteRange != null) {

            List<Integer> sites = Lists.newArrayList();

            for (Range<Integer> r : args.siteRange.asRanges()) {
                for (int s : ContiguousSet.create(r, DiscreteDomain.integers()).asList()) {
                    sites.add(s);
                }
            }

            CoreUtils.msg("Site range specified %s -> %s site(s) to analyse.\n", args.siteRange.toString(), sites.size());

            args.sites = GuavaUtils.getColumnSubset(args.sites, sites);

            Iterator<Map.Entry<Integer,Fitness>> iterator = args.fitnesses.entrySet().iterator();
            while (iterator.hasNext()) {
                Map.Entry<Integer,Fitness> entry = iterator.next();
                if(!sites.contains(entry.getKey())){
                    iterator.remove();
                }
            }

            // If we're doing site range, we cannot optimise mutational or branch length parameters
            if (args.fix.size() == 0) {
                args.fix.add("branches");
                args.fix.add("mutation");
            }

            // We also change the robust estimation interval (because we're not estimating everything)
            args.robustFitnessEstimationInterval = 1;
        }

        // by default, each site has weight = 1
        Map<Integer, Integer> patternWeights = Maps.newHashMap();
        for (int i : args.sites.columnKeySet()) patternWeights.put(i, 1);

        // by default each site's pattern is defined by itself
        Map<Integer, Integer> patternSiteMap = Maps.newHashMap();
        for (int i : args.sites.columnKeySet()) patternSiteMap.put(i, i);

        if (args.dontUseSitePatterns) {
            CoreUtils.msg("Not using site patterns.\n");
        } else {
            // Reduce the sites we analyse by using site patterns

            List<Integer> allSites = Lists.newArrayList(args.sites.columnKeySet());
            // step through sites, from end to beginning
            for (int i = allSites.size() - 1; i >= 0; i--) {
                int siteI = allSites.get(i);
                // compare this site with all other sites (with which is hasn't been compared)
                for (int j = i - 1; j >= 0; j--) {
                    int siteJ = allSites.get(j);
                    if (siteI == siteJ) continue;
                    if (areColumnsEqual(args.sites.column(siteI), args.sites.column(siteJ))) {
                        // the two sites have identical patterns
                        patternWeights.put(siteJ, patternWeights.get(siteJ) + patternWeights.get(siteI));
                        patternWeights.put(siteI, 0); // this site now has weight 0
                        patternWeights.remove(siteI); // we only need the weight of the sites we'll actually use!
                        patternSiteMap.put(siteI, siteJ); // siteI is the same as siteJ

                        // update the site -> pattern mapping for any sites that were pointing to siteI
                        for (int k = allSites.size() - 1; k > j; k--) {
                            int siteK = allSites.get(k);
                            if (patternSiteMap.get(siteK) == siteI) {
                                patternSiteMap.put(siteK, siteJ);
                            }
                        }
                        break;
                    }
                }
            }

            Collection<Integer> distinctSites = Lists.newArrayList();
            for (Map.Entry<Integer, Integer> e : patternWeights.entrySet()) {
                if (e.getValue() > 0) distinctSites.add(e.getKey());
            }

            CoreUtils.msg("%s observed site patterns.\n", distinctSites.size());
            args.sites = GuavaUtils.getColumnSubset(args.sites, distinctSites);
        }

        args.patternSiteMap = patternSiteMap;
        args.patternWeight = patternWeights;

    }

    public boolean areColumnsEqual(Map<String, Byte> column1, Map<String, Byte> column2) {
        // assume that the keys are identical, check every key has same value
        for (String s : column1.keySet()) {
            if (!column1.get(s).equals(column2.get(s))) return false;
        }
        return true;
    }

    public Table<String, Integer, Byte> getCleanedSitesTable(Alignment alignment) {
        // Prepare alignment as codon sites, replace unknown and STOP codons with gaps
        List<String> sequences = Lists.newArrayList();
        for (int i = 0; i < alignment.getSequenceCount(); i++)
            sequences.add(alignment.getIdentifier(i).getName());

        List<Integer> sites = Ints.asList(CoreUtils.seqi(1, alignment.getSiteCount() / 3));

        ArrayTable<String, Integer, Byte> cleanedSites = ArrayTable.create(sequences, sites);

        for (int i : sites) {
            Map<String, Integer> site = PhyloUtils.getCleanedCodons(alignment, i);
            for (Map.Entry<String, Integer> codon : site.entrySet()) {
                cleanedSites.put(codon.getKey(), i, codon.getValue().byteValue());
            }
        }
        return cleanedSites;
    }

    private boolean valid(Arguments args) {
        CoreUtils.msg("Summary: %s\n", args.summary());

        CoreUtils.msg("Tree has %s tips.\n", args.tree.getExternalNodeCount());

        // Check alignment and tree are congruent
        if (!PhyloUtils.isComplementary(args.tree, args.sites.rowKeySet())) {
            CoreUtils.msg("ERROR: tree and alignment do not have the same taxa.");
            return false;
        }

        // If any _fitness parameters were provided
        if (args._fitness != null && args._fitness.size() > 0) {
            if (args._fitness.size() / 21 != args.sites.columnKeySet().size()) {
                CoreUtils.msg("ERROR: alignment has %s sites but %s fitness vectors given.", args.sites.columnKeySet().size(), args._fitness.size() / 21);
                return false;
            }
        }

        /*
        // If we're running clade model, all parameters should have been provided
        if (args.cladeModel.size() > 0) {
            if (args._fitness == null || args.kappa == 0 || args.scaling == 0 || args.pi == null
                    || )
        }*/

        return true;
    }

}
