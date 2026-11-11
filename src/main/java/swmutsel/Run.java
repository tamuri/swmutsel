package swmutsel;

import swmutsel.cli.*;
import swmutsel.utils.CoreUtils;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 13:58
 */
public class Run {
    public static void main(String[] args) {
        new Run().dispatch(args);
    }

    private void dispatch(String... arguments) {
        Arguments args = ArgumentsProcessor.parse(arguments);
        CoreUtils.msg("Run ID %s\n", ((SharedArguments) args).identifier);
        args.getHandler().invoke();
    }
}
