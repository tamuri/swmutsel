package swmutsel.model;


import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import swmutsel.Constants;

public class MVNPenalty implements Penalty {

    private static final long serialVersionUID = 2735656410360111302L;

    public static final String LABEL = "mvn";

    private final double coefficient;
    private final double sigma;

    public MVNPenalty(final double sigma) {
        this.sigma = sigma;
        this.coefficient = 1 / (2 * Math.pow(sigma, 2));
    }

    @Override
    public double calculate(final double[] parameters) {
        // 20 fitness values
        DoubleMatrix1D f = DoubleFactory1D.dense.make(20);

        // 19 parameters are estimated, passed in parameters[]
        for (int i = 0; i < parameters.length; i++) f.setQuick(i, parameters[i]);

        // One is fixed
        f.setQuick(19, Constants.FITNESS_FIXED_FOR_RELATIVE);

        double mean = f.zSum() / 20.0;
        f.assign(Functions.minus(mean));
        return -coefficient * f.aggregate(Functions.plus, Functions.pow(2.0));
    }

    @Override
    public String toString() {
        return LABEL + "," + sigma;
    }
}

