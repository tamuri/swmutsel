package swmutsel.model.parameters;

import org.apache.commons.math3.util.MathArrays;
import swmutsel.utils.CoreUtils;

/**
 * A general purpose parameter of n-probabilities summing to 1
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 02/04/2014 13:55
 */
public class Probabilities extends Parameter {
    private static final long serialVersionUID = 3323655769041648532L;
    private double[] probabilities;

    public Probabilities(int n) {
        this(CoreUtils.rep(1.0, n));
    }

    public Probabilities(double[] probabilities) {
        set(probabilities);
    }

    public double[] get() {
        return this.probabilities;
    }

    public void set(double[] x) {
        this.probabilities = MathArrays.normalizeArray(x, 1);
    }

    @Override
    public double[] getOptimisable() {
        double[] x = new double[probabilities.length - 1];
        System.arraycopy(this.probabilities, 0, x, 0, x.length);
        return CoreUtils.alr(x);
    }

    @Override
    public void setOptimisable(double[] params) {
        double[] in = CoreUtils.alr_inv(params);
        double sum = 0;
        for (int i = 0; i < in.length; i++) {
            this.probabilities[i] = in[i];
            sum += in[i];
        }
        this.probabilities[this.probabilities.length - 1] = 1 - sum;
    }

    @Override
    public int getOptimisableCount() {
        return this.probabilities.length - 1;
    }

    @Override
    public String toString() {
        return "Probabilities{" +
                "p=" +
                CoreUtils.join("%.7f", ", ", this.probabilities) +
                '}';
    }
}
