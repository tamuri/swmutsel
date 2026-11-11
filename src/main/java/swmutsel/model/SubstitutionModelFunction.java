package swmutsel.model;

import org.apache.commons.math3.analysis.MultivariateFunction;
import swmutsel.Constants;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Parameter;
import swmutsel.runner.Runner;
import swmutsel.utils.CoreUtils;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 21/03/2014 17:25
 */
public class SubstitutionModelFunction  implements MultivariateFunction {
    private int evaluations = 0;
    private double current = Double.NEGATIVE_INFINITY;
    private final SubstitutionModel model;
    private final Runner runner;

    public SubstitutionModelFunction(SubstitutionModel model, Runner runner) {
        this.model = model;
        this.runner = runner;
    }

    @Override
    public double value(double[] point) {
        Mapper.setOptimisable(model.getParameters(), point);

        for (Parameter p : model.getParameters()) {
            if (p instanceof Fitness) {
                for (double d : ((Fitness)p).get()) {
                    if (d < -(Constants.FITNESS_BOUND + 1) || d > (Constants.FITNESS_BOUND + 1)) {
                        return Constants.VERY_BAD_LIKELIHOOD;
                    }
                }
            }
        }

        model.build();

        double lnL = CoreUtils.sum(runner.getLogLikelihood(model).values());

        evaluations++;

        // Running output
        if (evaluations == 1) {
            System.out.printf("%.7f ", lnL);
            // System.out.println(model.toString());
        }
        /*if (evaluations % 25 == 0) {
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
        }*/

        return lnL;
    }

    public int getEvaluations() {
        return evaluations;
    }
}
