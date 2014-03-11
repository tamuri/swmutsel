package swmutsel;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.CommaParameterSplitter;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.RangeSet;
import com.google.common.collect.Table;
import com.google.common.primitives.Doubles;
import pal.tree.Tree;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.model.Penalty;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.options.DoubleArrayConverter;
import swmutsel.options.GeneticCodeConverter;
import swmutsel.options.PenaltyConverter;
import swmutsel.options.SiteRangeConverter;
import swmutsel.utils.GeneticCode;

import java.util.List;
import java.util.Map;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 11:35
 */
@Parameters(separators = "=")
public class Arguments {

    @Parameter(names = "-help", required = false, hidden = true)
    public boolean showHelp = false;

    // for all models - required parameters

    @Parameter(names = {"-s", "-sequences"}, required = true)
    public String _alignmentPath;

    @Parameter(names = {"-t", "-tree"}, required = true)
    public String treePath;

    @Parameter(names = {"-gc", "-geneticcode"}, required = true, converter = GeneticCodeConverter.class)
    public GeneticCode geneticCode;

/*    @Parameter(names = {"-m", "-model"}, required = true)
    public String model;*/

    @Parameter(names = {"-n", "-name"}, required = true)
    public String identifier;


    // swMutSel parameters (-m swmutsel)

    @Parameter(names = {"-k", "-kappa"})
    public double kappa = TsTvRatio.getDefault();

    @Parameter(names = {"-c", "-scaling"})
    public double scaling = BranchScaling.getDefault();

    @Parameter(names = "-pi", converter = DoubleArrayConverter.class)
    public double[] pi = BaseFrequencies.getDefault();

    @Parameter(names = {"-p", "-penalty"}, converter = PenaltyConverter.class)
    public Penalty penalty = null;

    // Variable arity can be "all string up to next valid option". See https://github.com/cbeust/jcommander/pull/151
    @Parameter(names = {"-f", "-fitness"}, variableArity = true)
    public List<String> _fitness;

    @Parameter(names = "-clademodel", splitter = CommaParameterSplitter.class)
    public List<String> cladeModel = Lists.newArrayList("ALL");

    @Parameter(names = "-nopatterns")
    public boolean dontUseSitePatterns = false;


    // optimisation run options

    @Parameter(names = "-fix", splitter = CommaParameterSplitter.class)
    public List<String> fix = Lists.newArrayList();

    @Parameter(names = {"-T", "-threads"})
    public int threads = 1;

    @Parameter(names = { "-D", "-distributed"})
    public boolean distributed = false;

    @Parameter(names = {"-H","-hosts"}, splitter = CommaParameterSplitter.class)
    public List<String> hosts = Lists.newArrayList();

    @Parameter(names = "-sites", converter = SiteRangeConverter.class)
    public RangeSet<Integer> siteRange;

    @Parameter(names = "-restart-opt")
    public int numberOfOptimRestarts = 1;

    @Parameter(names = "-restart-int")
    public int robustFitnessEstimationInterval = 5;

    @Parameter(names = "-hessian")
    public boolean hessian = false;


    // To be filled in post-processing of arguments - see ArgumentsProcessor
    public Tree initialTree; // may be rooted
    public Tree tree;
    public FitnessStore fitnesses = new FitnessStore();
    public Table<String, Integer, Byte> sites;
    public Map<Integer, Integer> patternSiteMap;
    public Map<Integer, Integer> patternWeight;

    public String summary() {

        return String.format(
                "-kappa %s -scaling %s -pi %s -fix %s%s%s%s",
                kappa,
                scaling,
                Doubles.join(",", pi),
                Joiner.on(',').join(fix),
                distributed ? " -distributed" : "",
                (threads > 1) ? " -threads " + threads : "",
                (hosts.size() > 0) ? " -hosts " + Joiner.on(" -hosts ").join(hosts) : ""
        );

    }
}


