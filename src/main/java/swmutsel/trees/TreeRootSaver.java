package swmutsel.trees;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeManipulator;
import pal.tree.TreeUtils;
import swmutsel.utils.Pair;
import swmutsel.utils.PhyloUtils;
import swmutsel.utils.Triple;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 13/03/2014 22:19
 */
public class TreeRootSaver {
    private final List<String> outGroups;
    private final Node originalRootNode;
    private final Map<Set<String>, Pair<Set<String>, Set<String>>> descendantsOfNodes;
    private Triple<Set<String>, Set<String>, Set<String>> orderAtRoot;

    private static final int CHILD_0 = 0;
    private static final int CHILD_1 = 1;
    private static final int CHILD_2 = 2;
    private static final int UNROOTED_DEGREE = 3;
    private static final int ROOTED_DEGREE = 2;


    public TreeRootSaver(Tree t) {

        descendantsOfNodes = Maps.newHashMap();

        // Make sure we give back the tree with exactly the same node ordering (i.e. children appear in same order)
        for (Node node : PhyloUtils.internalNodes(t)) {

            Set<String> child0 = Sets.newHashSet();
            getDescendants(node.getChild(CHILD_0), child0);
            Set<String> child1 = Sets.newHashSet();
            getDescendants(node.getChild(CHILD_1), child1);

            Set<String> nodeKey = Sets.newHashSet();
            nodeKey.addAll(child0);
            nodeKey.addAll(child1);

            descendantsOfNodes.put(nodeKey, Pair.of(child0, child1));

            if (node.isRoot()) {
                Set<String> child2;
                if (node.getChildCount() == UNROOTED_DEGREE) {
                    child2 = Sets.newHashSet();
                    getDescendants(node.getChild(CHILD_2), child2);
                    nodeKey.addAll(child2);
                } else {
                    child2 = null;
                }
                orderAtRoot = Triple.of(child0, child1, child2);
            }
        }

        outGroups = Lists.newArrayList();

        // If the root is bifurcating, tree is ROOTED
        if (t.getRoot().getChildCount() == ROOTED_DEGREE) {
            getDescendants(t.getRoot().getChild(CHILD_0), outGroups);
            originalRootNode = null;
        } else {
            originalRootNode = t.getRoot();
        }

    }

    public void getDescendants(Node node, Collection<String> names) {
        if (node.isLeaf()) names.add(node.getIdentifier().getName());

        for (int i = 0; i < node.getChildCount(); i++) {
            getDescendants(node.getChild(i), names);
        }
    }


    /**
     * If we have an out-group, returns the tree with the original rooting.
     */
    public Tree getOriginalRooting(Tree tree) {
        if (outGroups.size() > 0) {
            String[] outGroupNames = outGroups.toArray(new String[outGroups.size()]);
            Tree rootedTree = TreeManipulator.getRootedBy(tree, outGroupNames);
            orderNodesInTree(rootedTree);
            return rootedTree;
        } else {
            TreeUtils.reroot(tree, originalRootNode);
            orderNodesInTree(tree);
            return tree;
        }
    }

    private void orderNodesInTree(Tree t) {
        for (Node n : PhyloUtils.internalNodes(t)) {

            Node n0 = n.getChild(CHILD_0);
            Set<String> child0 = Sets.newHashSet();
            getDescendants(n0, child0);

            Node n1 = n.getChild(CHILD_1);
            Set<String> child1 = Sets.newHashSet();
            getDescendants(n1, child1);

            Set<String> nodeKey = Sets.newHashSet();
            nodeKey.addAll(child0);
            nodeKey.addAll(child1);

            if (n.getChildCount() == UNROOTED_DEGREE) {
                Set<String> child2 = Sets.newHashSet();
                getDescendants(n.getChild(CHILD_2), child2);
                nodeKey.addAll(child2);
            }

            for (Set<String> s : descendantsOfNodes.keySet()) {
                if (s.equals(nodeKey)) {
                    if (n.getChildCount() == UNROOTED_DEGREE) {
                        if (!n.isRoot())
                            throw new RuntimeException("ERROR: Trifurcating node isn't root. I don't know what to do!");


                        Node n2 = n.getChild(CHILD_2);

                        n.removeChild(0);
                        n.removeChild(0);
                        n.removeChild(0);

                        Set<String> originalChild0 = orderAtRoot.first;
                        Set<String> originalChild1 = orderAtRoot.second;
                        Set<String> originalChild2 = orderAtRoot.third;

                        if (originalChild0.equals(child0)) {
                            n.insertChild(n0, CHILD_0);
                        } else if (originalChild0.equals(child1)) {
                            n.insertChild(n1, CHILD_0);
                        } else {
                            n.insertChild(n2, CHILD_0);
                        }

                        if (originalChild1.equals(child0)) {
                            n.insertChild(n0, CHILD_1);
                        } else if (originalChild1.equals(child1)) {
                            n.insertChild(n1, CHILD_1);
                        } else {
                            n.insertChild(n2, CHILD_1);
                        }

                        if (originalChild2.equals(child0)) {
                            n.insertChild(n0, CHILD_2);
                        } else if (originalChild2.equals(child1)) {
                            n.insertChild(n1, CHILD_2);
                        } else {
                            n.insertChild(n2, CHILD_2);
                        }

                    } else {
                        if (!child0.equals(descendantsOfNodes.get(s).first)) {
                            n.removeChild(0);
                            n.removeChild(0);

                            n.insertChild(n1, CHILD_0);
                            n.insertChild(n0, CHILD_1);
                        }
                    }
                }
            }
        }
    }



}
