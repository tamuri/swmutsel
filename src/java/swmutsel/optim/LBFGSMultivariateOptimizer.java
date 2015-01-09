package swmutsel.optim;

import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.OptimizationData;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import pal.math.MachineAccuracy;
import riso.numerical.LBFGS;
import swmutsel.Constants;
import swmutsel.utils.MathUtils;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 12/03/2014 10:08
 */
public class LBFGSMultivariateOptimizer extends MultivariateOptimizer {
    /**
     * @param checker Convergence checker.
     */
    public LBFGSMultivariateOptimizer(ConvergenceChecker<PointValuePair> checker) {
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
        double[] x = getStartPoint();
        int n = x.length;
        int m = 5;

        double[] diag = new double[n];

        int iprint[] = new int[2];
        iprint [ 1 -1] = 1;
        iprint [ 2 -1] = 0;

        int iflag[] = new int[1];
        iflag[0]=0;

        int icall = 0;

        double eps= Constants.CONVERGENCE_TOL;
        double xtol = MachineAccuracy.EPSILON;

        double f;

        do {

            f = -function.value(x);
            double[] g = MathUtils.gradient(getFunction(), x);
            try {
                LBFGS.lbfgs(n, m, x, f, g, false, diag, iprint, eps, xtol, iflag);
            } catch (LBFGS.ExceptionWithIflag e) {
                System.out.printf("x = %s\n", Doubles.join(",", x));
                throw new RuntimeException(e);
            }

            icall++;

        } while (iflag[0] != 0 && icall <= 500);

        return new PointValuePair(x, -f);

    }

    @Override
    public int getEvaluations() {
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

}
