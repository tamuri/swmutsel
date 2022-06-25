package swmutsel.trees;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeUtils;
import swmutsel.utils.Pair;
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
        // blc.branchLengthByNodeHeight(args[0], args[1]);
    }

    private void branchLengthByNodeHeight(String tree1Path, String tree2Path) {
        Tree tree1 = PhyloUtils.readTree(tree1Path);
        Tree tree2 = PhyloUtils.readTree(tree2Path);

        Map<Set<String>, Node> tree1Nodes = getNodeLeaves(tree1);
        Map<Set<String>, Node> tree2Nodes = getNodeLeaves(tree2);


        for (Map.Entry<Set<String>, Node> entry : tree1Nodes.entrySet()) {
            Node n1 = entry.getValue();
            Set<String> key1 = entry.getKey();

            Node n2 = tree2Nodes.get(key1);

            if (!n1.isRoot())
                System.out.printf("%s,%s\n", n1.getNodeHeight(), n1.getBranchLength());

        }

        System.out.println();
        System.out.println();

        for (Map.Entry<Set<String>, Node> entry : tree1Nodes.entrySet()) {
            Node n1 = entry.getValue();
            Set<String> key1 = entry.getKey();

            Node n2 = tree2Nodes.get(key1);

            if (!n1.isRoot())
            System.out.printf("%s,%s\n", n1.getNodeHeight(), n2.getBranchLength());

        }


    }


    private void run(String tree1Path, String tree2Path) {
        //Tree tree1 = TreeManipulator.getUnrooted(PhyloUtils.readTree(tree1Path));
        //Tree tree2 = TreeManipulator.getUnrooted(PhyloUtils.readTree(tree2Path));

        Tree tree1 = PhyloUtils.readTree(tree1Path);
        Tree tree2 = PhyloUtils.readTree(tree2Path);

        // The key for every node in the tree is the group of external nodes which are its children
        Map<Set<String>, Node> tree1Nodes = getNodeLeaves(tree1);
        Map<Set<String>, Node> tree2Nodes = getNodeLeaves(tree2);

        // Sanity check - the trees to have the same topology
        // System.out.printf("Tree1Nodes: %s\n", tree1Nodes.size());
        // System.out.printf("Tree2Nodes: %s\n", tree2Nodes.size());
        // System.out.printf("Intersection: %s\n", Sets.intersection(tree1Nodes.keySet(), tree2Nodes.keySet()).size());
        // System.out.printf("Difference: %s\n", Sets.difference(tree1Nodes.keySet(), tree2Nodes.keySet()).size());
        // System.out.println();



        // External node distances
        Map<String, Integer> tree1NodeNumber = Maps.newHashMap();
        Map<String, Integer> tree2NodeNumber = Maps.newHashMap();
        for (Node n : PhyloUtils.externalNodes(tree1)) {
            tree1NodeNumber.put(n.getIdentifier().getName(), n.getNumber());
        }
        for (Node n : PhyloUtils.externalNodes(tree2)) {
            tree2NodeNumber.put(n.getIdentifier().getName(), n.getNumber());
        }

        //System.out.println("ExternalDistance:");
        for (int i = 0; i < tree1.getExternalNodeCount(); i++) {
            Node nFrom = tree1.getExternalNode(i);
            String nameFrom = nFrom.getIdentifier().getName();
            int tree1From = tree1NodeNumber.get(nameFrom);
            int tree2From = tree2NodeNumber.get(nameFrom);

            for (int j = i + 1; j < tree1.getExternalNodeCount(); j++) {

                Node nTo = tree1.getExternalNode(j);
                String nameTo = nTo.getIdentifier().getName();
                int tree1To = tree1NodeNumber.get(nameTo);
                int tree2To = tree2NodeNumber.get(nameTo);

                double tree1Distance = TreeUtils.computeDistance(tree1, tree1From, tree1To);
                double tree2Distance = TreeUtils.computeDistance(tree2, tree2From, tree2To);

                // System.out.printf("  - [%s, %s, %.7f, %.7f]%n",
                System.out.printf("%s, %s, %.7f, %.7f%n",
                        nameFrom,
                        nameTo,
                        tree1Distance,
                        tree2Distance
                        );

            }

        }

        System.out.println( );
         System.exit(0);






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


        // balance the branches off the root node
        double totalRootBL = tree2.getRoot().getChild(0).getBranchLength();
        totalRootBL += tree2.getRoot().getChild(1).getBranchLength();
        /*System.out.println(NodeUtils.getLeafCount(tree2.getRoot().getChild(0)));
        System.out.println(NodeUtils.getLeafCount(tree2.getRoot().getChild(1)));
        System.out.println(totalRootBL);
*/
        double half = totalRootBL / 2;
        tree2.getRoot().getChild(0).setBranchLength(totalRootBL);
        tree2.getRoot().getChild(1).setBranchLength(0);
  //      System.out.println(half);





        System.out.println("BranchLengths:");
        for (Set<String> group : keys) {
            Node node1 = tree1Nodes.get(group);
            Node node2 = tree2Nodes.get(group);
            Pair<Integer,Double> node1Step = nodeFromRoot(node1);
            Pair<Integer,Double> node2Step = nodeFromRoot(node2);

            System.out.printf("  - [ %s, %s, %.7f, %.7f, %s, %s, \"%s\", %s ]%n",
                    node1.getBranchLength(), node2.getBranchLength(),
                    node1Step.second, node2Step.second,
                    node1Step.first, node2Step.first,
                    group,
                    group.size());

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

    private Pair<Integer,Double> nodeFromRoot(Node node) {
        Node parent = node;

        int step = 0;
        double height = 0;
        while (parent != null && !parent.isRoot()) {
            step++;
            height += parent.getBranchLength();
            parent = parent.getParent();
        }

        return Pair.of(step,height);
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
y <- yaml.load_file(pipe('pbpaste'))
yy <- data.frame(matrix(unlist(y$BranchLengths), nrow=length(y$BranchLengths), byrow=T))
yy$X3 <- NULL
plot(as.numeric(as.vector(yy$X2)))

x <- read.csv(pipe('pbpaste'), head=F, sep='\t')
x$V4 <- x$V1 / sum(x$V1) * 100
x$V5 <- x$V2 / sum(x$V2) * 100
x$V6 <- abs(x$V4 - x$V5)
par(pty='s')
plot(x$V4, x$V5)
abline(0,1)


 */
