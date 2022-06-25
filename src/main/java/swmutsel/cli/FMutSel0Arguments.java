package swmutsel.cli;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.converters.CommaParameterSplitter;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.collect.Table;
import pal.alignment.Alignment;
import pal.tree.Tree;
import swmutsel.Constants;
import swmutsel.FMutSelHandler;
import swmutsel.Handler;
import swmutsel.cli.jc.DoubleArrayConverter;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.PhyloUtils;

import java.util.List;
import java.util.Map;
import java.util.Set;

import static swmutsel.cli.ArgumentWriter.stringField;
import static swmutsel.cli.ArgumentWriter.stringFitnesses;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 23/03/2014 14:33
 */
public class FMutSel0Arguments extends SharedArguments {

    @Parameter(names = {"-s", "-sequences"}, required=true)
    public String alignmentPath;

    @Parameter(names = {"-k", "-kappa"})
    public double kappa = TsTvRatio.getDefault();

    @Parameter(names = "-pi", converter = DoubleArrayConverter.class)
    public double[] pi;

    @Parameter(names = "-omega")
    public double omega = Omega.getDefault();

    @Parameter(names = "-classes")
    public int classes = 1;

    @Parameter(names = "-weights", converter = DoubleArrayConverter.class)
    public double[] weights;

    @Parameter(names = {"-T", "-threads"})
    public int threads = 1;

    @Parameter(names = "-nopatterns")
    public boolean noPatterns = false;

    @Parameter(names = "-fix", hidden = true)
    public Set<String> fix = Sets.newHashSet();

    @Parameter(names = "-hessian")
    public boolean hessian = false;

    // Variable arity can be "all string up to next valid option". See https://github.com/cbeust/jcommander/pull/151
    // Not to be used in main code - only ArgumentsProcessor
    @Parameter(names = {"-f", "-fitness"}, variableArity = true)
    public List<String> fitnessList;

    @Parameter(names = { "-D", "-distributed"})
    public boolean distributed = false;

    @Parameter(names = {"-H","-hosts"}, splitter = CommaParameterSplitter.class)
    public List<String> hosts = Lists.newArrayList();

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

        if (!PhyloUtils.isComplementary(tree, this.sites.rowKeySet())) {
            String msg = "ERROR: tree and alignment do not have the same taxa.\n";
            throw new ParameterException(msg);
        }

        // if user didn't provide fitness vector for the [single class] FMutSel0 model
        if (this.fitnessList == null && this.classes == 1) {
            // populate with fitness derived from observed amino acids
            this.fitnesses.put(1, new Fitness(PhyloUtils.aminoAcidToFitness(PhyloUtils.getObservedAminoAcidFrequencies(this.sites), CoreUtils.seq(0, 19))));
        } else {
            this.fitnesses = ArgumentsProcessor.loadFitnessStore(this.fitnessList, this.classes, false, new ArgumentsProcessor.FitnessOrderingCallback() {
                @Override
                public int[] getOrderingForOptimiser(int set) {
                    return CoreUtils.seq(0, 19);
                }
            });
        }

        if (this.pi == null) {
            this.pi = PhyloUtils.getObservedNucleotideFrequencies(this.sites);
        }

        this.sites = ArgumentsProcessor.findPatterns(this.noPatterns, this.sites, palAlignment, this.patternSiteMap, this.patternWeight);

        treeRootedWarning(tree, Sets.<String>newHashSet());
    }

    public String getCheckpointString() {
        StringBuilder s = new StringBuilder();

        stringField(s, "-fmutsel0", null);

        stringField(s, "-sequences", this.alignmentPath);
        stringField(s, "-geneticcode", GeneticCode.getInstance().getCurrentCodeName());

        if (this.fix.contains("branches") || this.fix.contains("all")) {
            stringField(s, "-tree", this.treePath);
        } else {
            stringField(s, "-tree", this.tree.toString());
        }

        for (String fixed : this.fix) stringField(s, "-fix", fixed);

        stringField(s, "-kappa", String.format("%.7f", this.kappa));
        stringField(s, "-pi", String.format("%.7f,%.7f,%.7f,%.7f", this.pi[0], this.pi[1], this.pi[2], this.pi[3]));
        stringField(s, "-omega", String.format("%.7f", this.omega));

        if (this.classes > 1) {
            stringField(s, "-classes", Integer.toString(this.classes));
            StringBuilder w = new StringBuilder();
            for (int i = 0; i < this.weights.length; i++) {
                if (i > 0) w.append(",");
                w.append(String.format("%.7f", this.weights[i]));
            }
            stringField(s, "-weights", w.toString());
        }

        stringFitnesses(s, this.fitnesses, false);
        s.append("\n\n");
        s.append("//\n\n");
        return s.toString();
    }

    @Override
    public String getCommandSynopsisPath() {
        return Constants.HELP_COMMAND_FMUTSEL0;
    }

    @Override
    public Handler getHandler() {
        return new FMutSelHandler(this);
    }
}
