package swmutsel.trees;

import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeManipulator;
import pal.tree.TreeUtils;
import swmutsel.utils.PhyloUtils;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Takes a standard rooted/unrooted tree and handles the iteration over all the possible internal-node-rooted trees.
 * Importantly, it saves the original rooted (by saving an out-group) and restores it after iteration. This keeps
 * any new branch lengths if they've been updated outside this class.
 */
public class RerootedTreeIterator implements Iterable<Tree> {
    private static final int ROOTED_DEGREE = 2;

    private final Tree tree;
    private final Node[] internalNodes;
    private final TreeRootSaver treeRootSaver;

    public RerootedTreeIterator(Tree t) {
        treeRootSaver = new TreeRootSaver(t);

        if (t.getExternalNodeCount() == t.getRoot().getChildCount()) {
            this.tree = t;
            System.out.println(this.tree.toString());
        }

        else if (t.getRoot().getChildCount() == ROOTED_DEGREE) {
            this.tree = TreeManipulator.getUnrooted(t);

        } else {
            this.tree = t;
        }

        internalNodes = new Node[this.tree.getInternalNodeCount()];

        // Keep a list of all the internal nodes of the tree for iterator
        int pos = 0;
        for (Node node : PhyloUtils.internalNodes(this.tree)) {
            internalNodes[pos++] = node;
        }
    }

    public Tree getOriginalRooting() {
        return treeRootSaver.getOriginalRooting(this.tree);
    }

    @Override
    public Iterator<Tree> iterator() {
        return new Iterator<Tree>() {
            int position  = 0;

            @Override
            public boolean hasNext() {
                return position < internalNodes.length;
            }

            @Override
            public Tree next() {
                if(position > internalNodes.length - 1)
                    throw new NoSuchElementException();

                TreeUtils.reroot(tree, internalNodes[position]);
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
