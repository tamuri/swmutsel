package swmutsel.model;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 14/03/2014 04:50
 */
public interface TransProbCalculator {
    public void getTransitionProbabilities(final double[][] Pt, final double branchLength);
    public void calculatePartialLikelihood(final double[] Ldown, final double[] Lup, double t);
}
