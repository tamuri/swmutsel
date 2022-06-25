package swmutsel.trees;

import com.google.common.collect.Lists;
import com.google.common.primitives.Ints;
import pal.tree.*;
import swmutsel.utils.CoreUtils;

import java.io.PrintWriter;
import java.io.PushbackReader;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.Collections;
import java.util.List;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 04/06/2014 12:58
 */
public class TreeNodeRank {
    public static void main(String[] args) throws Exception{
        TreeNodeRank tnr = new TreeNodeRank();

        for (int i = 0; i < 10; i++) {
            tnr.run();

        }
    }

    private void run() throws Exception {
        String testTree = "( ( (a,b),((c,d),(j,k)) ), ( (e,f),(g,(h,i)) ) );";
        //String testTree = "( (a,b),(c,d) );";

        // TODO: test tree is rooted and bifurcating
        Tree tree = new ReadTree(new PushbackReader(new StringReader(testTree)));

        List<Integer> ranks = Lists.newArrayList(Ints.asList(CoreUtils.seq(0, tree.getExternalNodeCount() - 2)));

        assignRandomRanks(tree.getRoot(), ranks);

        NodeUtils.heights2Lengths(tree.getRoot(), true);

        StringWriter sw = new StringWriter();
        PrintWriter pw = new PrintWriter(sw);
        TreeUtils.printNH(tree, pw, true, true);
        pw.flush();
        pw.close();
        System.out.print(sw.toString());
        sw.close();
    }

    private void assignRandomRanks(Node node, List<Integer> availableRanks) {
        // Leaf nodes are not assigned a rank
        if (node.isLeaf()) return;

        Node left = node.getChild(0);
        Node right = node.getChild(1);

        // The rank of this node is the smallest possible rank in the list; set and remove
        int nodeRank = Collections.min(availableRanks);
        availableRanks.remove(availableRanks.indexOf(nodeRank));

        // Shuffle the remaining ranks
        Collections.shuffle(availableRanks);

        int leftRankSize;

        if (left.isLeaf())
            leftRankSize = 0;
        else
            leftRankSize = NodeUtils.getLeafCount(left) - 1;

        List<Integer> leftRanks = Lists.newArrayList(availableRanks.subList(0, leftRankSize));
        List<Integer> rightRanks = Lists.newArrayList(availableRanks.subList(leftRankSize, availableRanks.size()));

        assignRandomRanks(left, leftRanks);
        assignRandomRanks(right, rightRanks);

        node.getIdentifier().setName("" + nodeRank);

        if (!node.isRoot()) {
            double[] heights = new double[]{0.825017,0.7710149,0.6178068,0.5023164,0.4061649,0.3726295,0.2379134,0.08356747,0.04558889};
            node.setNodeHeight(heights[nodeRank - 1]);
        } else {
            node.setNodeHeight(1);

        }


    }

}
