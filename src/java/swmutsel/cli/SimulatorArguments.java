package swmutsel.cli;

import com.beust.jcommander.Parameter;
import com.google.common.collect.Lists;
import pal.tree.Tree;
import swmutsel.Constants;
import swmutsel.Handler;
import swmutsel.cli.jc.CharArrayConverter;
import swmutsel.cli.jc.DoubleArrayConverter;
import swmutsel.sim.SimulateHandler;
import swmutsel.utils.GeneticCode;

import java.util.List;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 23/03/2014 14:35
 */
public class SimulatorArguments extends SharedArguments {

    @Parameter(names = "-sites", required = false)
    public int sites = 1;

    // -fitness option can be specified multiple times for heterogeneous models
    @Parameter(names = "-fitness", required = false)
    public List<Double> fitness = Lists.newArrayList();

    // -fitnessfile option can be specified multiple times for heterogenous models
    @Parameter(names = "-fitnessfile", required = false)
    public List<String> fitnessFiles = Lists.newArrayList();

    // The fitnesses in -fitness can be specified in any order. Default is canonical amino acid order
    @Parameter(names = "-characters", converter = CharArrayConverter.class, required = false)
    public char[] residues = GeneticCode.AMINO_ACIDS;

    // You must specify this for heterogeneous models
    @Parameter(names = "-clademodel", required = false)
    public List<String> cladeModel = Lists.newArrayList("ALL");

    @Parameter(names = "-shiftfrac")
    public double shiftFrac = 0.5;

    @Parameter(names = "-kappa", required = true)
    public double kappa;

    @Parameter(names = "-pi", converter = DoubleArrayConverter.class, required = true)
    public double[] pi;

    @Parameter(names = "-scaling", required = true)
    public double scaling;

    public Tree tree;

    @Override
    public void initialise() {
        this.tree = ArgumentsProcessor.loadTree(this.treePath);
    }

    @Override
    public String getCommandSynopsisPath() {
        return Constants.HELP_COMMAND_SIMULATE;
    }

    @Override
    public Handler getHandler() {
        return new SimulateHandler(this);
    }
}
