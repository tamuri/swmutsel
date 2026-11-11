package swmutsel.model;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 14/03/2014 04:50
 */
public interface TransProbCalculator {
    void getTransitionProbabilities(final double[][] Pt, final double branchLength);
    void calculatePartialLikelihood(final double[] Ldown, final double[] Lup, double t);
    void calculatePartialLikelihoodPair(final double[] parent, final double[] leftChild, final double leftLength, final double[] rightChild, double rightLength, double[] pTemp, double[] qTemp);
}
