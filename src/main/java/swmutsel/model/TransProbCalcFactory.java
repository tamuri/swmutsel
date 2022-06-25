package swmutsel.model;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.LUDecompositionQuick;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import org.apache.commons.math3.util.FastMath;
import swmutsel.ArrayPool;
import swmutsel.colt.Algebra;
import swmutsel.utils.ColtMatrixUtils;
import swmutsel.utils.GeneticCode;

import java.util.Arrays;

/**
 * Most expensive function in the entire app - get rid of JVM array bounds
 * checking by having classes with number of sense codons hard-coded.
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 14/03/2014 04:41
 */
class TransProbCalcFactory {
// -------------------------- STATIC METHODS --------------------------

    static TransProbCalculator getPadeCalculator(double[] Q, int[] senseCodons, double scaling) {
        return new PadeApproximantCalculator(Q, senseCodons, scaling);
    }

    static TransProbCalculator getTaylorCalculator(double[] Q, int[] senseCodons, double scaling) {
        return new TaylorExpansionCalculator(Q, senseCodons, scaling);
    }

    static TransProbCalculator getDecompositionCalculator(double[] lambda, double[][] U, double[][] UInv, int[] senseCodons, boolean cache) {
        TransProbCalculator calculator = new DecompositionCalculator(lambda, U, UInv, senseCodons);
        if (cache) {
            return new Cached(calculator);
        }
        return calculator;
    }

    private static void matrixVectorMult(final double[][] leftM, final double[] leftV, final double[] rightV) {
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            final double[] rowLeftM = leftM[i];
            double b = 0;
            for (int j = 0; j < GeneticCode.CODON_STATES; j++){
                b += leftV[j] * rowLeftM[j];
            }
            rightV[i] *= b;
        }
    }
// -------------------------- INNER CLASSES --------------------------

    /**
     * You can cache calculated transition probabilities if the substitution model
     * does not change from site-to-site (i.e. the usual case). Approximate space
     * required a bit more than 64*64 matrix of 8-byte doubles * (number of taxa - 1 or 2)
     * e.g. 200 taxa ~ 6.25 mb, 512 taxa ~ 16 mb, 4096 ~ 128 mb
     *
     */
    private static final class Cached extends CacheLoader<Double, double[][]> implements TransProbCalculator {
        private final TransProbCalculator calculator;
        private final LoadingCache<Double, double[][]> cache;

        private Cached(TransProbCalculator calculator) {
            this.calculator = calculator;
            this.cache = CacheBuilder.newBuilder().maximumSize(1024).build(this);
        }

        @Override
        public void getTransitionProbabilities(final double[][] Pt, final double branchLength) {
            double[][] T = this.cache.getUnchecked(branchLength);

            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                System.arraycopy(T[i], 0, Pt[i], 0, GeneticCode.CODON_STATES);
            }
        }

        @Override
        public void calculatePartialLikelihood(final double[] Ldown, final double[] Lup, final double t) {
            final double[][] Pt = ArrayPool.popStatesByStatesMatrix();
            getTransitionProbabilities(Pt, t);
            matrixVectorMult(Pt, Ldown, Lup);
            ArrayPool.pushStatesByStatesMatrix(Pt);
        }

        @Override
        public void calculatePartialLikelihoodPair(double[] parent, double[] leftChild, double leftLength, double[] rightChild, double rightLength, double[] a, double[] b) {
            throw new RuntimeException("Not implemented");
        }

        @Override
        public double[][] load(final Double branchLength) throws Exception {
            double[][] Pt = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            this.calculator.getTransitionProbabilities(Pt, branchLength);
            return Pt;
        }
    }

    private static final class DecompositionCalculator implements TransProbCalculator {
        private final double[] lambda;
        private final double[][] U;
        private final double[][] UInv;
        private final double[][] UInv_cm; // UInv stored in column-major style
        private final double[][] U_cm;

        private DecompositionCalculator(final double[] lambda_in, final double[][] U_in, final double[][] UInv_in, final int[] senseCodons) {
            this.lambda = new double[GeneticCode.CODON_STATES];
            this.U = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            this.UInv = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            this.UInv_cm = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            this.U_cm = new double[64][64];

            // input arguments are in 'sense codon' vectors - store in 64-states
            for (int i = 0; i < senseCodons.length; i++) {
                this.lambda[senseCodons[i]] = lambda_in[i];
                double[] U_r = this.U[senseCodons[i]];
                double[] UInv_r = this.UInv[senseCodons[i]];
                double[] U_in_r = U_in[i];
                double[] UInv_in_r = UInv_in[i];
                for (int j = 0; j < senseCodons.length; j++) {
                    U_r[senseCodons[j]] = U_in_r[j];
                    UInv_r[senseCodons[j]] = UInv_in_r[j];
                    UInv_cm[senseCodons[j]][senseCodons[i]] = UInv_in_r[j];
                    U_cm[senseCodons[j]][senseCodons[i]] = U_in_r[j];
                }
            }
        }

/*        public final double[][] getTransitionProbabilities(final double branchLength) {
            final double[][] Pt = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            getTransitionProbabilities(Pt, branchLength);
            return Pt;
        }*/

        @Override
        public final void getTransitionProbabilities(final double[][] Pt, final double branchLength) {
            for (int i = 0; i < GeneticCode.CODON_STATES; i++)
                Arrays.fill(Pt[i], 0.0);

            for (int k = 0; k < GeneticCode.CODON_STATES; k++) {

                final double L = FastMath.exp(branchLength * lambda[k]);
                final double[] U_r = U_cm[k];
                final double[] UInv_r = UInv[k];

                for (int i = 0; i < GeneticCode.CODON_STATES; i++) {

                    final double[] Pt_r = Pt[i];
                    final double UL = U_r[i] * L;

                    for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                        Pt_r[j] += UL * UInv_r[j];
                    }
                }
            }
        }

        @Override
        public final void calculatePartialLikelihood(final double[] Ldown, final double[] Lup, final double t) {

            double[] T = new double[GeneticCode.CODON_STATES]; // TODO: only one needed per-thread

            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                double[] UInv_r = UInv[i];
                double sum = 0;
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    sum += UInv_r[j] * Ldown[j];
                }
                T[i] = sum * FastMath.exp(lambda[i] * t);
            }

            matrixVectorMult(U, T, Lup);

/*            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                double[] U_c = U[i];
                double sum = 0;
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    sum += U_c[j] * T[j];
                }
                Lup[i] *= sum;
            }*/

        }

/*
        public final void
        calculatePartialLikelihoodPair(final double[] parent,
                                       final double[] pCond, final double pLen,
                                       final double[] qCond, final double qLen) {
            final double[] pTemp = new double[GeneticCode.CODON_STATES];
            final double[] qTemp = new double[GeneticCode.CODON_STATES];
            calculatePartialLikelihoodPair(parent, pCond, pLen, qCond, qLen, pTemp, qTemp);
        }
*/

        @Override
        public final void
        calculatePartialLikelihoodPair(final double[] parent,
                                       final double[] pCond, final double pLen,
                                       final double[] qCond, final double qLen,
                                       final double[] pTemp, final double[] qTemp) {

            Arrays.fill(pTemp, 0);
            Arrays.fill(qTemp, 0);

            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {

                final double[] UInv_i = UInv[i];
                double p = 0, q = 0;
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    p += UInv_i[j] * pCond[j];
                    q += UInv_i[j] * qCond[j];
                }

                final double pC = p * FastMath.exp(lambda[i] * pLen);
                final double qC = q * FastMath.exp(lambda[i] * qLen);

                final double[] U_i = U_cm[i];
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    pTemp[j] += U_i[j] * pC;
                    qTemp[j] += U_i[j] * qC;
                }
            }

            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                parent[i] = pTemp[i] * qTemp[i];
            }
        }
    }

    private static final class TaylorExpansionCalculator implements TransProbCalculator {
        // Scaling & squaring parameters (using Taylor expansion)
        private final int j = 6;
        private final int m = (int) Math.pow(2, j);

        private final double scaling;
        private final DoubleMatrix2D Q1 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES, 0);
        private final DoubleMatrix2D Q2 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private final DoubleMatrix2D Q3 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private final DoubleMatrix2D PtOut = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private final DoubleMatrix2D tmp64x64 =  DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);

        TaylorExpansionCalculator(final double[] Q, final int[] senseCodons, final double scaling) {
            this.scaling = scaling;

            for (int i1 = 0; i1 < senseCodons.length; i1++) {
                for (int i2 = 0; i2 < senseCodons.length; i2++) {
                    this.Q1.setQuick(senseCodons[i1], senseCodons[i2], Q[i1 * senseCodons.length + i2]);
                }
            }

            tmp64x64.assign(this.Q1);
            ColtMatrixUtils.taylorExpTerm(tmp64x64, 2, Q2);
            tmp64x64.assign(this.Q1);
            ColtMatrixUtils.taylorExpTerm(tmp64x64, 3, Q3);
        }

        @Override
        public void getTransitionProbabilities(final double[][] Pt, final double branchLength) {
            final double f1 = branchLength * scaling / m;
            final double f2 = Math.pow(f1, 2);
            final double f3 = Math.pow(f1, 3);

            // Step 1, Taylor expansion
            // exp(Qt / m) ~ I + (Qt / m) + 1/2! (Qt / m)^2 + 1/3! (Qt / m)^3
            // exp(Qt / m) ~ I + Q1 * (t/m) + 1/2! Q1^2 * (t/m)^2 + 1/3! Q1^3 * (t/m)^3
            for (int i1 = 0; i1 < GeneticCode.CODON_STATES; i1++) {
                for (int i2 = 0; i2 < GeneticCode.CODON_STATES; i2++) {
                    PtOut.setQuick(i1, i2,
                            Q1.getQuick(i1, i2) * f1 +
                                    Q2.getQuick(i1, i2) * f2 +
                                    Q3.getQuick(i1, i2) * f3 +
                                    (i1 == i2 ? 1 : 0));
                }
            }

            // Step 2: exp(Qt) = exp(Qt/m)^m, where m = 2^j
            for (int i = 0; i < j; i++) {
                PtOut.zMult(PtOut, tmp64x64);
                PtOut.assign(tmp64x64);
            }

            for (int i1 = 0; i1 < GeneticCode.CODON_STATES; i1++) {
                for (int i2 = 0; i2 < GeneticCode.CODON_STATES; i2++) {
                    Pt[i1][i2] = PtOut.getQuick(i1, i2);
                }
            }
        }

        @Override
        public void calculatePartialLikelihood(final double[] Ldown, final double[] Lup, final double t) {
            final double[][] Pt = ArrayPool.popStatesByStatesMatrix();
            getTransitionProbabilities(Pt, t);
            matrixVectorMult(Pt, Ldown, Lup);
            ArrayPool.pushStatesByStatesMatrix(Pt);
        }

        @Override
        public void calculatePartialLikelihoodPair(double[] parent, double[] leftChild, double leftLength, double[] rightChild, double rightLength, double[] pTemp, double[] qTemp) {
            throw new RuntimeException("not implemented");
        }
    }

    private static final class PadeApproximantCalculator implements TransProbCalculator {

        // WARNING: This class is not thread-safe!!
        private static final DoubleMatrix2D Q1 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private static final DoubleMatrix2D Q2 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private static final DoubleMatrix2D Q3 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private static final DoubleMatrix2D Q4 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private static final DoubleMatrix2D Q5 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private static final DoubleMatrix2D PT = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private static final DoubleMatrix2D TMP =  DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private static final DoubleMatrix2D D = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);
        private static final LUDecompositionQuick LU = new LUDecompositionQuick(0);

        private final double scaling;
        private final double Qnorm;

        // private final double[] c = new double[]{1.0, 0.5, 0.1, 0.008333333333333333};
        // private final double[] c = new double[]{1.0, 0.5, 0.10714285714285714, 0.011904761904761904, 5.952380952380953E-4};
        private final double[] c = new double[]{1.0, 0.5, 0.1111111111111111, 0.013888888888888888, 9.92063492063492E-4, 3.306878306878307E-5};

        PadeApproximantCalculator(final double[] Qin, final int[] senseCodons, final double scaling) {
            this.scaling = scaling;

            for (int i1 = 0; i1 < senseCodons.length; i1++) {
                for (int i2 = 0; i2 < senseCodons.length; i2++) {
                    Q1.setQuick(senseCodons[i1], senseCodons[i2], Qin[i1 * senseCodons.length + i2]);
                }
            }

            // Pre-calculate powers of Q1
            Q1.zMult(Q1, Q2);
            Q2.zMult(Q1, Q3);
            Q3.zMult(Q1, Q4);
            Q4.zMult(Q1, Q5);

            Qnorm = Algebra.DEFAULT.norm1(Q1);
        }

        @Override
        public void getTransitionProbabilities(final double[][] Pt, final double branchLength) {
            getTransitionProbabilities(branchLength);

            for (int i1 = 0; i1 < GeneticCode.CODON_STATES; i1++) {
                for (int i2 = 0; i2 < GeneticCode.CODON_STATES; i2++) {
                    Pt[i1][i2] = PT.getQuick(i1, i2);
                }
            }
        }

        private void getTransitionProbabilities(final double branchLength) {

            // Based on Table 1 in Moler & Van Loan (2003), resulting in lnL values within 1e-3 of
            // those using spectral decomposition method on real data sets (ha & mtcyb)
            final double norm = Qnorm * branchLength * scaling; // use ||Qt|| to decide how many squarings
            int j;
            if (norm <= 1) {
                j = 0;
            } else if (norm <= 10) {
                j = 1;
            } else if (norm <= 100) {
                j = 5;
            } else {
                j = 8;
            }
            final double m = Math.pow(2, j);

            final double a1 = branchLength * scaling / m;
            final double a2 = a1 * a1;
            final double a3 = a1 * a2;
            final double a4 = a1 * a3;
            final double a5 = a1 * a4;

            for (int i1 = 0; i1 < GeneticCode.CODON_STATES; i1++) {
                for (int i2 = 0; i2 < GeneticCode.CODON_STATES; i2++) {
                    final double diag = i1 == i2 ? 1.0 : 0.0;
                    final double f1 = Q1.getQuick(i1, i2) * a1 * c[1];
                    final double f2 = Q2.getQuick(i1, i2) * a2 * c[2];
                    final double f3 = Q3.getQuick(i1, i2) * a3 * c[3];
                    final double f4 = Q4.getQuick(i1, i2) * a4 * c[4];
                    final double f5 = Q5.getQuick(i1, i2) * a5 * c[5];

                    D.setQuick(i1, i2, diag + -f1 + f2 + -f3 + f4 + -f5);

                    // building Nqq here
                    PT.setQuick(i1, i2, diag + f1 + f2 + f3 + f4 + f5);
                }
            }

            // Instead of doing x = A^(-1)B, we solve Ax=B
            LU.decompose(D); // in: Dqq out: L-U so L * U = D
            LU.solve(PT); // in: Nqq matrix; out: x = D^(-1)N (from solving Dx=N)

            // Step 2: exp(Qt) = exp(Qt/m)^m, where m = 2^j
            for (int i = 0; i < j; i++) {
                PT.zMult(PT, TMP);
                PT.assign(TMP);
            }
        }

        @Override
        public void calculatePartialLikelihood(final double[] Ldown, final double[] Lup, final double t) {
            getTransitionProbabilities(t);

            double b1, b2, b3, b4;
            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                {
                    b1 = b2 = b3 = b4 = 0;
                    b1 += Ldown[0]  * PT.getQuick(i, 0);
                    b1 += Ldown[1]  * PT.getQuick(i, 1);
                    b1 += Ldown[2]  * PT.getQuick(i, 2);
                    b1 += Ldown[3]  * PT.getQuick(i, 3);
                    b1 += Ldown[4]  * PT.getQuick(i, 4);
                    b1 += Ldown[5]  * PT.getQuick(i, 5);
                    b1 += Ldown[6]  * PT.getQuick(i, 6);
                    b1 += Ldown[7]  * PT.getQuick(i, 7);
                    b1 += Ldown[8]  * PT.getQuick(i, 8);
                    b1 += Ldown[9]  * PT.getQuick(i, 9);
                    b1 += Ldown[10] * PT.getQuick(i,10);
                    b1 += Ldown[11] * PT.getQuick(i,11);
                    b1 += Ldown[12] * PT.getQuick(i,12);
                    b1 += Ldown[13] * PT.getQuick(i,13);
                    b1 += Ldown[14] * PT.getQuick(i,14);
                    b1 += Ldown[15] * PT.getQuick(i,15);
                    b2 += Ldown[16] * PT.getQuick(i,16);
                    b2 += Ldown[17] * PT.getQuick(i,17);
                    b2 += Ldown[18] * PT.getQuick(i,18);
                    b2 += Ldown[19] * PT.getQuick(i,19);
                    b2 += Ldown[20] * PT.getQuick(i,20);
                    b2 += Ldown[21] * PT.getQuick(i,21);
                    b2 += Ldown[22] * PT.getQuick(i,22);
                    b2 += Ldown[23] * PT.getQuick(i,23);
                    b2 += Ldown[24] * PT.getQuick(i,24);
                    b2 += Ldown[25] * PT.getQuick(i,25);
                    b2 += Ldown[26] * PT.getQuick(i,26);
                    b2 += Ldown[27] * PT.getQuick(i,27);
                    b2 += Ldown[28] * PT.getQuick(i,28);
                    b2 += Ldown[29] * PT.getQuick(i,29);
                    b2 += Ldown[30] * PT.getQuick(i,30);
                    b2 += Ldown[31] * PT.getQuick(i,31);
                    b3 += Ldown[32] * PT.getQuick(i,32);
                    b3 += Ldown[33] * PT.getQuick(i,33);
                    b3 += Ldown[34] * PT.getQuick(i,34);
                    b3 += Ldown[35] * PT.getQuick(i,35);
                    b3 += Ldown[36] * PT.getQuick(i,36);
                    b3 += Ldown[37] * PT.getQuick(i,37);
                    b3 += Ldown[38] * PT.getQuick(i,38);
                    b3 += Ldown[39] * PT.getQuick(i,39);
                    b3 += Ldown[40] * PT.getQuick(i,40);
                    b3 += Ldown[41] * PT.getQuick(i,41);
                    b3 += Ldown[42] * PT.getQuick(i,42);
                    b3 += Ldown[43] * PT.getQuick(i,43);
                    b3 += Ldown[44] * PT.getQuick(i,44);
                    b3 += Ldown[45] * PT.getQuick(i,45);
                    b3 += Ldown[46] * PT.getQuick(i,46);
                    b3 += Ldown[47] * PT.getQuick(i,47);
                    b4 += Ldown[48] * PT.getQuick(i,48);
                    b4 += Ldown[49] * PT.getQuick(i,49);
                    b4 += Ldown[50] * PT.getQuick(i,50);
                    b4 += Ldown[51] * PT.getQuick(i,51);
                    b4 += Ldown[52] * PT.getQuick(i,52);
                    b4 += Ldown[53] * PT.getQuick(i,53);
                    b4 += Ldown[54] * PT.getQuick(i,54);
                    b4 += Ldown[55] * PT.getQuick(i,55);
                    b4 += Ldown[56] * PT.getQuick(i,56);
                    b4 += Ldown[57] * PT.getQuick(i,57);
                    b4 += Ldown[58] * PT.getQuick(i,58);
                    b4 += Ldown[59] * PT.getQuick(i,59);
                    b4 += Ldown[60] * PT.getQuick(i,60);
                    b4 += Ldown[61] * PT.getQuick(i,61);
                    b4 += Ldown[62] * PT.getQuick(i,62);
                    b4 += Ldown[63] * PT.getQuick(i,63);
                }
                Lup[i] *= b1 + b2 + b3 + b4;
            }
        }

        @Override
        public void calculatePartialLikelihoodPair(double[] parent, double[] leftChild, double leftLength, double[] rightChild, double rightLength, double[] pTemp, double[] qTemp) {
            throw new RuntimeException("not implemented");
        }
    }
}
