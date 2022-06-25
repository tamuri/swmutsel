package swmutsel.cli;

import pal.tree.Tree;
import swmutsel.Handler;
import swmutsel.utils.CoreUtils;

import java.util.Set;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 23/03/2014 15:01
 */
public abstract class Arguments {
    public abstract void initialise();

    public abstract Handler getHandler();

    public abstract String getCommandSynopsisPath();

    protected static void treeRootedWarning(Tree tree, Set<String> fix) {
        if (tree.getRoot().getChildCount() == 2 &&
                !(fix.contains("branches") || fix.contains("all"))) {
            CoreUtils.msg("WARNING: This is a rooted tree and you are optimising branch lengths. The tree :WARNING%n");
            CoreUtils.msg("WARNING: will be unrooted during branch length optimisation and the program    :WARNING%n");
            CoreUtils.msg("WARNING: will try to put back the root afterwards. Please check the output     :WARNING%n");
            CoreUtils.msg("WARNING: carefully.                                                            :WARNING%n");
        }

    }
}
