package swmutsel.model;

import org.apache.commons.math3.analysis.MultivariateFunction;
import swmutsel.Constants;
import swmutsel.model.parameters.Mapper;
import swmutsel.runner.Runner;
import swmutsel.utils.CoreUtils;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 12/03/2014 10:40
 */
public class SwMutFunction implements MultivariateFunction {
    private int evaluations = 0;
    private double current = Double.NEGATIVE_INFINITY;
    private final SwMut mutation;
    private final Runner runner;

    public SwMutFunction(SwMut mutation, Runner runner) {
        this.mutation = mutation;
        this.runner = runner;
    }

    @Override
    public double value(double[] point) {
        Mapper.setOptimisable(mutation.getParameters(), point);

        mutation.build();

        double lnL = CoreUtils.sum(runner.getLogLikelihood(mutation, false).values());

        evaluations++;

        // Running output
        if (evaluations == 1) System.out.printf("%.7f ", lnL);
        if (evaluations % 25 == 0) {
            if (lnL < current) {
                System.out.printf("↓"); // getting worse!
            } else {
                if ((lnL - current) > Constants.SMALL_CHANGE ) {
                    System.out.printf("↑"); // getting better
                } else {
                    System.out.printf("→"); // small improvement
                }
            }
            current = lnL;
        }
        return lnL;
    }

    public int getEvaluations() {
        return evaluations;
    }
}
