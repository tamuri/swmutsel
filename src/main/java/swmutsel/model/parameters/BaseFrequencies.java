package swmutsel.model.parameters;

import org.apache.commons.math3.util.MathArrays;
import swmutsel.utils.CoreUtils;

import java.util.Arrays;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 21:32
 */
public class BaseFrequencies extends Parameter {

    private static final long serialVersionUID = 2674842918979799977L;

    public BaseFrequencies() {
        setArgument("-pi");
    }

    private double[] frequencies;

    public BaseFrequencies(double[] pi) {
        this();
        double sum = CoreUtils.sum(pi);
        if (sum != 1)
            CoreUtils.msg("WARNING: Nucleotide frequencies [ %s ] do not sum to 1. Normalising.\n", CoreUtils.join("%.5f", ", ", pi));
        set(pi);
    }

    public static double[] getDefault() {
        return new double[]{0.25, 0.25, 0.25, 0.25};
    }

    public double[] get() {
        return this.frequencies;
    }

    public void set(double[] params) {
        // params.length == 4 && sum(params) == 1
        this.frequencies = MathArrays.normalizeArray(params, 1);
    }

    @Override
    public double[] getOptimisable() {
        return CoreUtils.alr(new double[]{
                this.frequencies[0],
                this.frequencies[1],
                this.frequencies[2]
        });

    }

    @Override
    public void setOptimisable(double[] params) {
        // params.length == 3
        double[] frequencies = CoreUtils.alr_inv(params);
        this.frequencies = new double[]{
                frequencies[0],
                frequencies[1],
                frequencies[2],
                1 - (frequencies[0] + frequencies[1] + frequencies[2])
        };

    }

    @Override
    public int getOptimisableCount() {
        return 3;
    }

    @Override
    public String toString() {
        return String.format("BaseFrequencies{frequencies=%s}", CoreUtils.join("%.7f", ", ", frequencies));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BaseFrequencies that = (BaseFrequencies) o;

        if (!Arrays.equals(frequencies, that.frequencies)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return Arrays.hashCode(frequencies);
    }
}
