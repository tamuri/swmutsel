package swmutsel.trees;

import com.google.common.base.Charsets;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.io.Files;
import pal.tree.Node;
import pal.tree.Tree;
import pal.tree.TreeTool;
import pal.tree.TreeUtils;
import swmutsel.utils.PhyloUtils;

import java.io.File;
import java.io.PrintWriter;
import java.util.List;
import java.util.Set;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class TreeNodeLabeler {
    public static void main(String[] args) throws Exception {
        TreeNodeLabeler m = new TreeNodeLabeler();
        m.run(args[0]);
    }

    private void run(String fileIn) throws Exception {
        System.out.printf("Reading file %s.\n", fileIn);
        Tree t = TreeTool.readTree(Files.newReader(new File(fileIn), Charsets.US_ASCII));

        List<Node> unknownNodes = Lists.newArrayList();
        List<Node> knownNodes = Lists.newArrayList();

        // to begin with all internal nodes are unnamed
        for (Node n : PhyloUtils.internalNodes(t)) unknownNodes.add(n);

        int preUnknownNodeCount = unknownNodes.size();
        int postUnknownNodeCount = -1;

        while (!unknownNodes.isEmpty() && preUnknownNodeCount != postUnknownNodeCount) {
            System.out.printf("Unknown nodes: %s.\n", unknownNodes.size());

            preUnknownNodeCount = unknownNodes.size();

            // loop through each of the unknown nodes
            for (Node n : unknownNodes) {
                // get a distinct list of all hosts for the node
                Set<String> s = Sets.newHashSet();

                // loop over every child node
                for (int j = 0; j < n.getChildCount(); j++) {
                    Node c = n.getChild(j);
                    String childNodeName = c.getIdentifier().getName();
                    if (c.isLeaf()) {
                        s.add(childNodeName.substring(0, 2));
                    } else {
                        // one of the children is an internal node with a name
                        if (childNodeName != null && childNodeName.length() > 0) {
                            s.add(childNodeName.substring(0, 2));
                        }
                    }
                }

                // list all the possible names for this node
                System.out.printf("Node %s possible clade(s): %s\n", n.getNumber(), Joiner.on(" ").join(s));

                if (s.size() == 1) {
                    n.getIdentifier().setName(s.iterator().next());
                    knownNodes.add(n);
                }
            }

            // remove all known nodes
            for (Node n : knownNodes) {
                unknownNodes.remove(n);
            }

            postUnknownNodeCount = unknownNodes.size();
        }

        // now the only remaining unknown nodes are those that neighbour two different clades (and the root node)
        for (Node n : unknownNodes) {
            // the name of this hostshift node = parent node name (remember, could be the root node = leave unlabelled)
            if (n.getParent() == null) {
                n.getIdentifier().setName("ROOT");
            } else {
                String parentNodeName = n.getParent().getIdentifier().getName();
                if (parentNodeName != null && parentNodeName.length() > 0) {
                    n.getIdentifier().setName(parentNodeName.substring(0, 2) + "_HS");
                    System.out.printf("Node %s has change of clade: %s\n", n.getNumber(), n.getIdentifier().getName());
                }
            }
        }

        PrintWriter pw = new PrintWriter(Files.newWriter(new File(fileIn + ".out"), Charsets.US_ASCII));
        TreeUtils.printNH(t, pw);
        pw.close();
    }

}

