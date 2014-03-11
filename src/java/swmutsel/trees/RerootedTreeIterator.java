package swmutsel.trees;

import com.beust.jcommander.internal.Lists;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeManipulator;
import pal.tree.TreeUtils;
import swmutsel.utils.PhyloUtils;

import java.util.Iterator;
import java.util.List;

/**
 * Takes a standard rooted/unrooted tree and handles the iteration over all the possible internal-node-rooted trees.
 * Importantly, it saves the original rooted (by saving an outgroup) and restores it after iteration. This keeps
 * any new branch lengths if they've been updated outside this class.
 */
public class RerootedTreeIterator implements Iterable<Tree> {
    private Tree tree;
    private List<Node> internalNodes = Lists.newArrayList();
    private List<String> outgroups = Lists.newArrayList();

    public RerootedTreeIterator(Tree t) {
        // Save the outgroup taxa from one child off the root node, so we can restore it once we're done
        if (t.getRoot().getChildCount() == 2) {
            findOutgroups(t.getRoot().getChild(0));
            //System.out.printf("Outgroup(s): %s\n\n", outgroups.toString());
        } else {
            //System.out.println("Root is trifurcation - no defined outgroup.");
        }

        // Unroot the tree - this is what we'll work with
        this.tree = TreeManipulator.getUnrooted(t);

        // Keep a list of all the internal nodes of the unrooted tree
        for (Node node : PhyloUtils.internalNodes(tree)) this.internalNodes.add(node);
    }

    /**
     * Descends the tree looking for a leaf nodes, and saves those to the outgroups list
     */
    public void findOutgroups(Node node) {
        if (node.isLeaf()) {
            outgroups.add(node.getIdentifier().getName());
        }

        for (int i = 0; i < node.getChildCount(); i++) {
            findOutgroups(node.getChild(i));
        }
    }

    /**
     * If we have an outgroup, returns the tree with the original rooting.
     */
    public Tree getOriginalRooting() {
        if (outgroups.size() > 0) {
            String[] outs = outgroups.toArray(new String[outgroups.size()]);
            return TreeManipulator.getRootedBy(tree, outs);
        } else {
            return tree;
        }
    }

    @Override
    public Iterator<Tree> iterator() {
        return new Iterator<Tree>() {
            int position  = 0;

            @Override
            public boolean hasNext() {
                return position < internalNodes.size();
            }

            @Override
            public Tree next() {
                TreeUtils.reroot(tree, internalNodes.get(position));
                position++;
                return tree;
            }

            @Override
            public void remove() {
                throw new RuntimeException("Not implemented!");
            }
        };
    }
}
