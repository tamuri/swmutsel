package swmutsel.results;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.collect.Lists;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.options.DoubleArrayConverter;
import swmutsel.options.GeneticCodeConverter;
import swmutsel.utils.GeneticCode;

import java.util.List;

/**
 * Runs all the various procedures for creating the distribution of selection coefficients from raw results files
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class All {
    private All() {
    }

    public static void main(String[] args) throws Exception {
        // Parse command-line options
        Options options = new Options();
        JCommander jc = new JCommander(options);
        jc.parse(args);


        // Write fitnesses to F.txt given the raw results file
        FitnessExtractor fe = new FitnessExtractor();
        fe.extract(options._fitness);

        // Write the model parameter values given the fitnesses and global parameters
        ModelWriter mw = new ModelWriter(options);
        mw.run();

        // Use the model parameter values to create the distribution of selection coefficients
        DistributionWriter dw = new DistributionWriter();
        dw.run();
    }

    /**
     * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
     */
    static class Options {

        @Parameter(description = "Misc")
        public List<String> _other = Lists.newArrayList();

        @Parameter(names = {"-gc", "-geneticcode"}, required = true, converter = GeneticCodeConverter.class)
        public GeneticCode geneticCode;

        @Parameter(names = {"-k", "-kappa"})
        public double kappa = TsTvRatio.getDefault();

        @Parameter(names = {"-scaling", "-c"})
        public double scaling = BranchScaling.getDefault();

        @Parameter(names = "-pi", converter = DoubleArrayConverter.class)
        public double[] pi = BaseFrequencies.getDefault();

        // Variable arity can be "all string up to next valid option". See https://github.com/cbeust/jcommander/pull/151
        @Parameter(names = "-fitness", variableArity = true)
        public List<String> _fitness;

    }
}
