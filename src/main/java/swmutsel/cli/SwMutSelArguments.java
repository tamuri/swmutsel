package swmutsel.cli;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.converters.CommaParameterSplitter;
import com.google.common.collect.*;
import com.google.common.io.Files;
import com.google.common.primitives.Ints;
import pal.alignment.Alignment;
import pal.tree.Tree;
import swmutsel.Constants;
import swmutsel.Handler;
import swmutsel.SwMutSelHandler;
import swmutsel.cli.jc.DoubleArrayConverter;
import swmutsel.cli.jc.PenaltyConverter;
import swmutsel.cli.jc.SiteRangeConverter;
import swmutsel.model.Penalty;
import swmutsel.model.parameters.*;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.PhyloUtils;

import java.io.BufferedReader;
import java.io.File;
import java.nio.charset.Charset;
import java.util.*;

import static swmutsel.cli.ArgumentWriter.stringField;
import static swmutsel.cli.ArgumentWriter.stringFitnesses;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 23/03/2014 14:33
 */
public class SwMutSelArguments extends SharedArguments {
    // NOTE: inherits arguments from SharedArguments!

    @Parameter(names = {"-s", "-sequences"}, required=true)
    public String alignmentPath;

    @Parameter(names = {"-c", "-scaling"})
    public double scaling = BranchScaling.getDefault();

    @Parameter(names = {"-p", "-penalty"}, converter = PenaltyConverter.class)
    public Penalty penalty = null;

    @Parameter(names = "-clademodel", splitter = CommaParameterSplitter.class)
    public List<String> cladeModel = Lists.newArrayList("ALL");

    @Parameter(names = {"-k", "-kappa"})
    public double kappa = TsTvRatio.getDefault();

    @Parameter(names = "-pi", converter = DoubleArrayConverter.class)
    public double[] pi = BaseFrequencies.getDefault();

    @Parameter(names = "-tau", hidden = true)
    public double tau = Tau.getDefault();

    // Variable arity can be "all string up to next valid option". See https://github.com/cbeust/jcommander/pull/151
    // Not to be used in main code - only ArgumentsProcessor
    @Parameter(names = {"-f", "-fitness"}, variableArity = true)
    public List<String> fitnessList;

    // Optimisation args
    @Parameter(names = "-fix")
    public Set<String> fix = Sets.newHashSet();

    @Parameter(names = {"-T", "-threads"})
    public int threads = 1;

    @Parameter(names = { "-D", "-distributed"})
    public boolean distributed = false;

    @Parameter(names = {"-H","-hosts"}, splitter = CommaParameterSplitter.class)
    public List<String> hosts = Lists.newArrayList();

    @Parameter(names = "-sites", converter = SiteRangeConverter.class)
    public RangeSet<Integer> siteRange;

    @Parameter(names = "-restart-opt")
    public int optRestart = 1;

    @Parameter(names = "-restart-int")
    public int optRestartInterval = 5;

    @Parameter(names = "-nopatterns")
    public boolean noPatterns = false;

    @Parameter(names = "-hessian")
    public boolean hessian = false;

    // To be filled in post-processing of arguments - see ArgumentsProcessor
    public Tree tree;
    public FitnessStore fitnesses = new FitnessStore();
    public Table<String, Integer, Byte> sites;
    public Map<Integer, Integer> patternSiteMap = Maps.newHashMap();
    public Map<Integer, Integer> patternWeight = Maps.newHashMap();

    @Override
    public void initialise() {
        this.tree = ArgumentsProcessor.loadTree(this.treePath);

        final Alignment palAlignment = ArgumentsProcessor.loadPalAlignment(this.alignmentPath);
        this.sites = ArgumentsProcessor.getAlignmentTable(palAlignment);

        // Older versions of the program allowed comma-separated fix lists
        Set<String> split = Sets.newHashSet();
        Iterator<String> fixedIterator = this.fix.iterator();
        while (fixedIterator.hasNext()) {
            String s = fixedIterator.next();
            if (s.contains(",")) {
                Collections.addAll(split, s.split(","));
                fixedIterator.remove();
            }
        }
        if (split.size() > 0) this.fix.addAll(split);

        if (!PhyloUtils.isComplementary(tree, this.sites.rowKeySet())) {
            throw new ParameterException("ERROR: tree and alignment do not have the same taxa.\n");
        }

        // fitnesses are in a separate file to be loaded (useful for very large alignments)
        if (fitnessList != null && fitnessList.size() == 1) {
            final String fitnessFile = this.fitnessList.get(0);
            try {
                BufferedReader br = Files.newReader(new File(fitnessFile), Charset.defaultCharset());

                fitnessList.clear();

                String line;
                while ((line = br.readLine()) != null) {
                    String[] S = line.split(",");

                    if (S.length != 21)
                        throw new RuntimeException("Can't parse fitness for '" + line + "'");

                    Collections.addAll(fitnessList, S);
                }

            } catch (Exception e) {
                throw new RuntimeException("Can't load fitnesses from file '" + this.fitnessList.get(0) + "'\n" + e.getMessage());
            }
        }

        this.fitnesses = ArgumentsProcessor.loadFitnessStore(this.fitnessList, this.sites.columnKeySet().size(), true, new ArgumentsProcessor.FitnessOrderingCallback() {
            @Override
            public int[] getOrderingForOptimiser(int set) {
                // The SwMutSel0 model orders amino acids by their observed frequency
                return Ints.toArray(PhyloUtils.getOrderedAminoAcids(sites.column(set).values()));
            }
        });


        this.sites = ArgumentsProcessor.findPatterns(this.noPatterns, this.sites, palAlignment, this.patternSiteMap, this.patternWeight);

        // Trim the sites and fitnesses to only those we need
        if (this.siteRange != null) {
            this.sites = ArgumentsProcessor.loadSiteRange(this.siteRange, this.sites, this.fitnesses);

            // If we're doing site range, we cannot optimise mutational or branch length parameters
            this.fix.add("branches");
            this.fix.add("mutation");

            // We also change the robust estimation interval (because we're not estimating everything)
            this.optRestartInterval = 1;
        }

        treeRootedWarning(tree, this.fix);
    }

    public String getCheckpointString() {
        StringBuilder s = new StringBuilder();

        stringField(s, "-sequences", this.alignmentPath);
        stringField(s, "-geneticcode", GeneticCode.getInstance().getCurrentCodeName());

        if (this.fix.contains("branches") || this.fix.contains("all")) {
            stringField(s, "-tree", this.treePath);
        } else {
            stringField(s, "-tree", this.tree.toString());
        }

        stringField(s, "-kappa", String.format("%.7f", this.kappa));
        stringField(s, "-pi", String.format("%.7f,%.7f,%.7f,%.7f", this.pi[0], this.pi[1], this.pi[2], this.pi[3]));
        stringField(s, "-scaling", String.format("%.7f", this.scaling));

        if (this.penalty != null) stringField(s, "-penalty", this.penalty.toString());

        for (String fixed : this.fix) stringField(s, "-fix", fixed);

        stringFitnesses(s, this.fitnesses, true);
        s.append("\n\n");
        s.append("//\n\n");
        return s.toString();
    }

    @Override
    public String getCommandSynopsisPath() {
        return Constants.HELP_COMMAND_SWMUTSEL;
    }

    @Override
    public Handler getHandler() {
        return new SwMutSelHandler(this);
    }
}
