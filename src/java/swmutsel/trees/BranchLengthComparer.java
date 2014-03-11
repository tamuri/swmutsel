package swmutsel.trees;

import com.beust.jcommander.internal.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeManipulator;
import swmutsel.utils.PhyloUtils;

import java.util.*;

/**
 * Takes two trees and compares their branch lengths. They can be rooted differently but must have identical taxa,
 * which is identified by its string label.
 */
public class BranchLengthComparer {
    public static void main(String[] args) {
        BranchLengthComparer blc = new BranchLengthComparer();
        blc.run(args[0], args[1]);
    }

    private void run(String tree1Path, String tree2Path) {
        Tree tree1 = TreeManipulator.getUnrooted(PhyloUtils.readTree(tree1Path));
        Tree tree2 = TreeManipulator.getUnrooted(PhyloUtils.readTree(tree2Path));

        // The key for every node in the tree is the group of external nodes which are its children
        Map<Set<String>, Node> tree1Nodes = getNodeLeaves(tree1);
        Map<Set<String>, Node> tree2Nodes = getNodeLeaves(tree2);

        // Sanity check - the trees to have the same topology
        System.out.printf("Tree1Nodes: %s\n", tree1Nodes.size());
        System.out.printf("Tree2Nodes: %s\n", tree2Nodes.size());
        System.out.printf("Intersection: %s\n", Sets.intersection(tree1Nodes.keySet(), tree2Nodes.keySet()).size());
        System.out.printf("Difference: %s\n", Sets.difference(tree1Nodes.keySet(), tree2Nodes.keySet()).size());
        System.out.println();

        // Just for pretty printing
        List<Set<String>> keys = Lists.newArrayList(tree1Nodes.keySet());
        Collections.sort(keys, new Comparator<Set<String>>() {
            @Override
            public int compare(Set<String> s1, Set<String> s2) {
                return s1.iterator().next().compareTo(s2.iterator().next());
            }
        });

        // Output the branch lengths of each matching node
        double totalLength1 = 0.0;
        double totalLength2 = 0.0;

        System.out.println("BranchLengths:");
        for (Set<String> group : keys) {
            Node node1 = tree1Nodes.get(group);
            Node node2 = tree2Nodes.get(group);

            System.out.printf("  - [ %s, %s, \"%s\" ]%n", node1.getBranchLength(), node2.getBranchLength(), group);

            totalLength1 += node1.getBranchLength();
            totalLength2 += node2.getBranchLength();
        }

        System.out.println();

        // Normalise both trees so total tree length is the same, so we can compare relative branch lengths
        double newTotal = 100.0;

        for (Set<String> group : keys) {
            tree1Nodes.get(group).setBranchLength(tree1Nodes.get(group).getBranchLength() / totalLength1 * newTotal);
            tree2Nodes.get(group).setBranchLength(tree2Nodes.get(group).getBranchLength() / totalLength2 * newTotal);
        }

        System.out.printf("FirstTree: >%n  %s%n", PhyloUtils.prettyTreeString(tree1).replaceAll("\n", "\n  "));
        System.out.printf("SecondTree: >%n  %s%n", PhyloUtils.prettyTreeString(tree2).replaceAll("\n", "\n  "));
    }

    private Map<Set<String>, Node> getNodeLeaves(Tree tree) {
        Map<Set<String>, Node> treeNodes = Maps.newHashMap();

        // Get taxa connect to internal nodes
        for (Node node : PhyloUtils.internalNodes(tree)) {
            Set<String> leaves = Sets.newTreeSet();
            getLeavesFromNode(node, leaves);
            treeNodes.put(leaves, node);
        }

        // Get the external nodes themselves
        for (Node node : PhyloUtils.externalNodes(tree)) {
            Set<String> leaves = Sets.newTreeSet();
            getLeavesFromNode(node, leaves);
            treeNodes.put(leaves, node);
        }

        return treeNodes;
    }

    private void getLeavesFromNode(Node node, Set<String> leaves) {
        if (node.isLeaf()) {
            leaves.add(node.getIdentifier().getName());
        } else {
            for (int i = 0; i < node.getChildCount(); i++) {
                Node child = node.getChild(i);
                getLeavesFromNode(child, leaves);
            }
        }
    }

}

/*

In R:
x <- read.csv(pipe('pbpaste'), head=F, sep='\t')
x$V4 <- x$V1 / sum(x$V1) * 100
x$V5 <- x$V2 / sum(x$V2) * 100
x$V6 <- abs(x$V4 - x$V5)
par(pty='s')
plot(x$V4, x$V5)
abline(0,1)


 */
