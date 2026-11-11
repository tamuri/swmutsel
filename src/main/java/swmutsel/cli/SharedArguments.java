package swmutsel.cli;

import com.beust.jcommander.Parameter;
import swmutsel.Handler;
import swmutsel.cli.jc.GeneticCodeConverter;
import swmutsel.utils.GeneticCode;

/**
 * Arguments that are either:
 * 1. Required
 * 2. For determining what action we're running
 *
 * <p/>
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 23/03/2014 14:32
 */
public class SharedArguments extends Arguments {
    @Parameter(names = "-help", required = false, hidden = true)
    public boolean showHelp = false;

    @Parameter(names = {"-t", "-tree"}, required = true)
    public String treePath;

    @Parameter(names = {"-gc", "-geneticcode"}, converter = GeneticCodeConverter.class, required = true)
    public GeneticCode geneticCode;

    @Parameter(names = {"-n", "-name"}, required = true)
    public String identifier;

    @Parameter(names = "-simulate", required = false)
    public boolean simulate = false;

    @Parameter(names = "-fmutsel0", required = false)
    public boolean fmutsel0 = false;

    @Parameter(names = "-fds", required = false)
    public boolean fdr = false;

    @Override
    public void initialise() {
        throw new RuntimeException("Not implemented.");
    }

    @Override
    public Handler getHandler() {
        throw new RuntimeException("Not implemented.");
    }

    @Override
    public String getCommandSynopsisPath() {
        return null; // show everything
    }
}


