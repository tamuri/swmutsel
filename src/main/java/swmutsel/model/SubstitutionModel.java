package swmutsel.model;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import com.google.common.collect.Lists;
import swmutsel.model.parameters.Parameter;
import swmutsel.utils.MathUtils;

import java.io.Serializable;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 15:23
 */
public abstract class SubstitutionModel implements Serializable {
    private static final long serialVersionUID = 6408557902068035305L;
    private final List<Parameter> parameters = Lists.newArrayList();

    protected static double getRelativeFixationProbability(final double S) {
        if (S == 0) return 1;
        if (S < -1e3) return 0;
        if (S > 1e3) return S;
        return S / -Math.expm1(-S);

//      return (S == 0) ? 1 : S / (1 - Math.exp(-S));
    }

    /**
     * B = Π^(½) * Q * Π^(-½) - Yang CME - pg. 69, eq. 2.26
     */
    protected static DoubleMatrix2D makeB(final double[] Q, final int numSenseCodons, final double[] senseCodonFrequencies) {
        DoubleMatrix2D B = DoubleFactory2D.dense.make(numSenseCodons, numSenseCodons);
        for (int i = 0; i < numSenseCodons; i++) {
            for (int j = 0; j < numSenseCodons; j++) {
                B.setQuick(i, j, Q[i * numSenseCodons + j] * Math.sqrt(senseCodonFrequencies[i]) / Math.sqrt(senseCodonFrequencies[j]));
            }
        }
        return B;
    }

    /**
     * Diagonalise B = R * Λ * R⁻¹ - Yang CME - pg. 69, eq. 2.27
     *
     * where B = Π^(½) * Q * Π^(-½)
     *
     * so
     *
     * Q = (Π^(-½) * R) * Λ * (R⁻¹ * Π^(½))       eq. 2.28
     *
     * and
     *
     * U   = Π^(-½) * R
     * U⁻¹ = R⁻¹ * Π^(½)
     *
     * Note, R⁻¹ == transpose(R)
     *
     * and so
     *
     * P(t) =  U * exp(Λ * t) * U⁻¹
     *
     */
    protected static void setEigenFactors(final DoubleMatrix2D B, final double[] lambda, final int states, final double[] freq, final double[][] U, final double[][] Uinv) {
        double[][] R = new double[states][states];

        // Symmetric models only
        MathUtils.setEigenFactors(B, lambda, R, states);
        double[] U_row, Uinv_row;

        for (int row = 0; row < states; row++) {
            double piInvSqrt = 1 / Math.sqrt(freq[row]);
            U_row = U[row];
            Uinv_row = Uinv[row];
            for (int column = 0; column < states; column++) {
                U_row[column] = piInvSqrt * R[row][column]; // U = Π^(-½) * R
                Uinv_row[column] = R[column][row] * Math.sqrt(freq[column]); // U⁻¹ = R⁻¹ * Π^(½)
            }
        }

        /*
        // For non-symmetric models but can't use spectral decomposition TODO: should be removed!
        EigenvalueDecomposition eigen = new EigenvalueDecomposition(B);
        DoubleMatrix2D RR = eigen.getV();
        DoubleMatrix1D LL = DoubleFactory2D.dense.diagonal(eigen.getD());
        double[] U_row, Uinv_row;
        DoubleMatrix2D diag = DoubleFactory2D.dense.make(states, states);
        for (int row = 0; row < states; row++) {
            double piInvSqrt = 1 / Math.sqrt(freq[row]);
            U_row = U[row];
            Uinv_row = Uinv[row];
            for (int column = 0; column < states; column++) {
                U_row[column] = piInvSqrt * RR.getQuick(row, column); // U = Π^(-½) * R
                Uinv_row[column] = RR.getQuick(column, row) * Math.sqrt(freq[column]); // U⁻¹ = R⁻¹ * Π^(½)
            }
            lambda[row] = LL.getQuick(row);
            diag.setQuick(row, row, lambda[row]);
        }
        */
    }

    public abstract TransProbCalculator getPtCalculator();

    public abstract double[] getCodonFrequencies();
    public abstract void build();
    public abstract boolean parametersValid();

    public void addParameters(Parameter... parameters) {
        Collections.addAll(this.parameters, parameters);
    }

    public void addParameters(Collection<? extends Parameter> parameters) {
        this.parameters.addAll(parameters);
    }

    public void clearParameters() {
        this.parameters.clear();
    }



    public List<Parameter> getParameters() {
        return this.parameters;
    }
}
