package swmutsel.model;


import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;

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
        DoubleMatrix1D f = DoubleFactory1D.dense.make(parameters);
        return -coefficient * f.aggregate(Functions.plus, Functions.pow(2.0));
    }

    @Override
    public String toString() {
        return LABEL + "," + sigma;
    }
}

