package swmutsel.trees;

import com.google.common.collect.Lists;
import pal.misc.Identifier;
import pal.tree.Node;
import pal.tree.NodeFactory;
import pal.tree.SimpleTree;
import swmutsel.utils.PhyloUtils;

import java.util.List;

import static swmutsel.utils.CoreUtils.seq;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 25/04/2014 21:19
 */
public class TreeCreator {
    public static void main(String[] args) {
        TreeCreator tc = new TreeCreator();
        tc.run(args[0], Integer.parseInt(args[1]), Double.parseDouble(args[2]));
    }

    private void run(String command, int totalNodes, double length) {

        if (command.equals("unbalanced_fixed")) {
            unbalanced_fixed(totalNodes, length);
        } else if (command.equals("unbalanced_grow")) {
        } else if (command.equals("balanced_fixed")) {
            balanced_fixed(totalNodes, length);
        } else if (command.equals("balanced_halved")) {
            balanced_halved(totalNodes, length);
        }
        else {
            throw new RuntimeException("Unknown command '" + command + "'");
        }

    }


    private void balanced_halved(int nodes, double height) {
        if (!isPowerOfTwo(nodes))
            throw new RuntimeException("Does not implement balanced tree where number_of_taxa != 2^n");

        Node root = NodeFactory.createNode();

        List<Node> currentTips = Lists.newArrayList(root);

        int totalTips = currentTips.size();
        int depth = 0;

        while (totalTips < nodes) {
            depth++;

            List<Node> newTips = Lists.newArrayList();

            for (Node n : currentTips) {

                double l;
                l = height / Math.pow(2, depth);

                if (depth == 7) l += l;

                Node n1 = NodeFactory.createNode();
                n1.setBranchLength(l);
                newTips.add(n1);

                Node n2 = NodeFactory.createNode();
                n2.setBranchLength(l );
                newTips.add(n2);

                n.addChild(n1);
                n.addChild(n2);

            }

            currentTips.clear();
            currentTips.addAll(newTips);

            totalTips = currentTips.size();
        }


        SimpleTree st = new SimpleTree(root);

        int seqNumber = 1;
        for (Node tip : PhyloUtils.externalNodes(st)) {
            tip.setIdentifier(new Identifier("seq_" + seqNumber++));
        }

        System.out.println(st.toString());
    }



    private boolean isPowerOfTwo(int x) {
        return (x != 0) && ((x & (x - 1)) == 0);
    }

    private void balanced_fixed(int nodes, double length) {

        if (!isPowerOfTwo(nodes))
            throw new RuntimeException("Does not implement balanced tree where number_of_taxa != 2^n");

        Node root = NodeFactory.createNode();

        List<Node> currentTips = Lists.newArrayList(root);


        int totalTips = currentTips.size();

        while (totalTips < nodes) {
            List<Node> newTips = Lists.newArrayList();

            for (Node n : currentTips) {
                Node n1 = NodeFactory.createNode();
                n1.setBranchLength(length);
                newTips.add(n1);

                Node n2 = NodeFactory.createNode();
                n2.setBranchLength(length);
                newTips.add(n2);

                n.addChild(n1);
                n.addChild(n2);
            }

            currentTips.clear();
            currentTips.addAll(newTips);

            totalTips = currentTips.size();
        }


        SimpleTree st = new SimpleTree(root);

        int seqNumber = 1;
        for (Node tip : PhyloUtils.externalNodes(st)) {
            tip.setIdentifier(new Identifier("seq_" + seqNumber++));
        }

        System.out.println(st.toString());


    }

    private void unbalanced_fixed(int totalNodes, double stepLength) {

        int nodeNumber = totalNodes;

        Node n1 = NodeFactory.createNode();
        n1.setIdentifier(new Identifier("seq_" + nodeNumber--));
        n1.setBranchLength(stepLength);

        Node n2 = NodeFactory.createNode();
        n2.setIdentifier(new Identifier("seq_" + nodeNumber--));
        n2.setBranchLength(stepLength);

        Node n3 = NodeFactory.createNode();
        n3.addChild(n1);
        n3.addChild(n2);

        Node root = n3;

        for (int i : seq(1, totalNodes - 2)) {

            root.setBranchLength(stepLength);

            Node n4 = NodeFactory.createNode();

            n4.addChild(root);

            root = n4;


            double length = stepLength * (i + 1);
            Node n5 = NodeFactory.createNode();
            n5.setBranchLength(length);
            n5.setIdentifier(new Identifier("seq_" + nodeNumber--));


            n4.addChild(n5);

        }

        SimpleTree st = new SimpleTree(root);
        System.out.println(st.toString());

    }
}
