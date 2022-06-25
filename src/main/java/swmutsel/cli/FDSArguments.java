package swmutsel.cli;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import pal.alignment.Alignment;
import pal.tree.Tree;
import swmutsel.ArrayPool;
import swmutsel.Constants;
import swmutsel.FDSHandler;
import swmutsel.Handler;
import swmutsel.cli.jc.DoubleArrayConverter;
import swmutsel.cli.jc.PenaltyConverter;
import swmutsel.model.Penalty;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.utils.PhyloUtils;

import java.util.List;
import java.util.Map;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 23/03/2014 14:33
 */
public class FDSArguments extends SharedArguments {
    // NOTE: inherits arguments from SharedArguments!

    @Parameter(names = {"-s", "-sequences"}, required=true)
    public String alignmentPath;

    @Parameter(names = {"-c", "-scaling"}, required=true)
    public double scaling = BranchScaling.getDefault();

    @Parameter(names = {"-p", "-penalty"}, converter = PenaltyConverter.class)
    public Penalty penalty = null;

    @Parameter(names = {"-k", "-kappa"}, required=true)
    public double kappa = TsTvRatio.getDefault();

    @Parameter(names = "-pi", converter = DoubleArrayConverter.class, required=true)
    public double[] pi = BaseFrequencies.getDefault();

    @Parameter(names = "-fitness", required = false)
    public List<Double> fitness = Lists.newArrayList();

    @Parameter(names = "-site", required = true)
    public int site;

    @Parameter(names = "-Z", required = false)
    public double Z;

    @Parameter(names = "-Zpenalty", required = false)
    public double Zpenalty;

    @Parameter(names = "-noOpt")
    public boolean noOpt = false;

    // To be filled in post-processing of arguments - see ArgumentsProcessor
    public Tree tree;
    public Table<String, Integer, Byte> sites;
    public Map<Integer, Integer> patternSiteMap = Maps.newHashMap();
    public Map<Integer, Integer> patternWeight = Maps.newHashMap();

    @Override
    public void initialise() {
        this.tree = ArgumentsProcessor.loadTree(this.treePath);

        final Alignment palAlignment = ArgumentsProcessor.loadPalAlignment(this.alignmentPath);
        this.sites = ArgumentsProcessor.getAlignmentTable(palAlignment);

        if (!PhyloUtils.isComplementary(tree, this.sites.rowKeySet())) {
            throw new ParameterException("ERROR: tree and alignment do not have the same taxa.\n");
        }

        ArrayPool.setTreeSize(tree.getInternalNodeCount() + 2);
    }

    @Override
    public String getCommandSynopsisPath() {
        return Constants.HELP_COMMAND_SWMUTSEL;
    }

    @Override
    public Handler getHandler() {
        return new FDSHandler(this);
    }
}
