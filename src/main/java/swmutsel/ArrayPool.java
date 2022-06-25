package swmutsel;

import com.google.common.collect.Lists;
import swmutsel.utils.GeneticCode;

import java.util.LinkedList;

/**
 * User: atamuri
 * Date: 01/03/2013 22:54
 */
// TODO: Move these pools to the classes where they are used: LikelihoodCalculator and TDGCodonModel!
public final class ArrayPool {

    private final static LinkedList<double[][]> NODES_BY_STATES_MATRIX = Lists.newLinkedList();
    private final static LinkedList<double[][]> STATES_BY_STATES_MATRIX = Lists.newLinkedList();

    private static int treeSize = -1; // This should be set by caller (n > 0)
    public static synchronized void setTreeSize(int size) {
        treeSize = size;
    }

    public static synchronized double[][] popNodesByStatesMatrix() {
        if (NODES_BY_STATES_MATRIX.isEmpty())
            return new double[treeSize][GeneticCode.CODON_STATES];
        else
            return NODES_BY_STATES_MATRIX.pop();
    }

    public static synchronized void pushNodesByStatesMatrix(double[][] matrix) {
        NODES_BY_STATES_MATRIX.push(matrix);
    }

    public static synchronized void pushStatesByStatesMatrix(double[][] matrix) {
        STATES_BY_STATES_MATRIX.push(matrix);
    }

    public static synchronized double[][] popStatesByStatesMatrix() {
        if (STATES_BY_STATES_MATRIX.isEmpty())
            return new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
        else
            return STATES_BY_STATES_MATRIX.pop();
    }
}
