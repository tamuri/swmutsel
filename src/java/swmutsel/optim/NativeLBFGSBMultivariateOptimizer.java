package swmutsel.optim;

import lbfgsb.*;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import swmutsel.Constants;
import swmutsel.utils.MathUtils;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 12/03/2014 10:08
 */
public class NativeLBFGSBMultivariateOptimizer extends MultivariateOptimizer {
    /**
     * @param checker Convergence checker.
     */
    public NativeLBFGSBMultivariateOptimizer(ConvergenceChecker<PointValuePair> checker) {
        super(checker);
    }

    private MultivariateFunction function;

    private MultivariateFunction getFunction() {
        return function;
    }

    private void setFunction(MultivariateFunction function) {
        this.function = function;
    }

    @Override
    public PointValuePair optimize(OptimizationData... optData) throws TooManyEvaluationsException {
        return super.optimize(optData);
    }

    @Override
    protected PointValuePair doOptimize() {

        try {
            FunctionWrapper fun = new FunctionWrapper(getFunction());
            Minimizer alg = new Minimizer();
            alg.setNoBounds(getStartPoint().length);
            // alg.getStopConditions().setFunctionReductionFactor(1e12);
            alg.setIterationFinishedListener(new DebugListener());
            Result result = alg.run(fun, getStartPoint());
            return new PointValuePair(result.point, -result.functionValue);
        } catch (LBFGSBException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }

    @Override
    public int getEvaluations() {
        // computeObjectValue() calls + gradient using finite differences
        return super.getEvaluations() + (2 * getStartPoint().length * super.getEvaluations());
    }

    @Override
    protected void parseOptimizationData(OptimizationData... optData) {
        super.parseOptimizationData(optData);
        for (OptimizationData data : optData) {
            if (data instanceof ObjectiveFunction) {
                setFunction(((ObjectiveFunction) data).getObjectiveFunction());
                continue;
            }
        }
    }

    private class FunctionWrapper implements DifferentiableFunction {
        private MultivariateFunction function;

        private FunctionWrapper(MultivariateFunction function) {
            this.function = function;
        }

        @Override
        public FunctionValues getValues(double[] point) {
            double v = -computeObjectiveValue(point);
            double[] grad = MathUtils.gradient(function, point);
            return new FunctionValues(v, grad);
        }
    }

    static private class DebugListener implements IterationFinishedListener {
        int i = 0;
        double current = Double.POSITIVE_INFINITY;
        @Override
        public boolean iterationFinished(double[] point,
                                         double functionValue, double[] gradient) {

            // Running output
            if (i == 0) {
                i++;
                return true;
            }

            if (functionValue > current) {
                System.out.printf("↓"); // getting worse!
            } else {
                if ((current - functionValue) > Constants.SMALL_CHANGE ) {
                    System.out.printf("↑"); // getting better
                } else {
                    System.out.printf("→"); // small improvement
                }
            }

            if (current - functionValue < Constants.CONVERGENCE_TOL) {
                return false;
            }

            current = functionValue;
            i++;

            return true;
        }
    }

}
