package swmutsel.model;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import org.apache.commons.math3.util.FastMath;
import swmutsel.ArrayPool;
import swmutsel.utils.GeneticCode;

/**
 * Most expensive function in the entire app - get rid of JVM array bounds
 * checking by having classes with number of sense codons hard-coded.
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 14/03/2014 04:41
 */
public class TransProbCalcFactory {
// -------------------------- STATIC METHODS --------------------------

    public static TransProbCalculator getPtCalculator(double[] lambda, double[][] U, double[][] UInv, int[] senseCodons, boolean cache) {
        TransProbCalculator calculator = new PtCalculator64(lambda, U, UInv, senseCodons);

        if (cache) {
            return new Cached(calculator);
        }

        return calculator;
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

            double b1, b2, b3, b4;
            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                final double[] probMatrix_r = Pt[i];
                {
                    b1 = b2 = b3 = b4 = 0;
                    b1 += Ldown[0] * probMatrix_r[0];
                    b1 += Ldown[1] * probMatrix_r[1];
                    b1 += Ldown[2] * probMatrix_r[2];
                    b1 += Ldown[3] * probMatrix_r[3];
                    b1 += Ldown[4] * probMatrix_r[4];
                    b1 += Ldown[5] * probMatrix_r[5];
                    b1 += Ldown[6] * probMatrix_r[6];
                    b1 += Ldown[7] * probMatrix_r[7];
                    b1 += Ldown[8] * probMatrix_r[8];
                    b1 += Ldown[9] * probMatrix_r[9];
                    b1 += Ldown[10] * probMatrix_r[10];
                    b1 += Ldown[11] * probMatrix_r[11];
                    b1 += Ldown[12] * probMatrix_r[12];
                    b1 += Ldown[13] * probMatrix_r[13];
                    b1 += Ldown[14] * probMatrix_r[14];
                    b1 += Ldown[15] * probMatrix_r[15];
                    b2 += Ldown[16] * probMatrix_r[16];
                    b2 += Ldown[17] * probMatrix_r[17];
                    b2 += Ldown[18] * probMatrix_r[18];
                    b2 += Ldown[19] * probMatrix_r[19];
                    b2 += Ldown[20] * probMatrix_r[20];
                    b2 += Ldown[21] * probMatrix_r[21];
                    b2 += Ldown[22] * probMatrix_r[22];
                    b2 += Ldown[23] * probMatrix_r[23];
                    b2 += Ldown[24] * probMatrix_r[24];
                    b2 += Ldown[25] * probMatrix_r[25];
                    b2 += Ldown[26] * probMatrix_r[26];
                    b2 += Ldown[27] * probMatrix_r[27];
                    b2 += Ldown[28] * probMatrix_r[28];
                    b2 += Ldown[29] * probMatrix_r[29];
                    b2 += Ldown[30] * probMatrix_r[30];
                    b2 += Ldown[31] * probMatrix_r[31];
                    b3 += Ldown[32] * probMatrix_r[32];
                    b3 += Ldown[33] * probMatrix_r[33];
                    b3 += Ldown[34] * probMatrix_r[34];
                    b3 += Ldown[35] * probMatrix_r[35];
                    b3 += Ldown[36] * probMatrix_r[36];
                    b3 += Ldown[37] * probMatrix_r[37];
                    b3 += Ldown[38] * probMatrix_r[38];
                    b3 += Ldown[39] * probMatrix_r[39];
                    b3 += Ldown[40] * probMatrix_r[40];
                    b3 += Ldown[41] * probMatrix_r[41];
                    b3 += Ldown[42] * probMatrix_r[42];
                    b3 += Ldown[43] * probMatrix_r[43];
                    b3 += Ldown[44] * probMatrix_r[44];
                    b3 += Ldown[45] * probMatrix_r[45];
                    b3 += Ldown[46] * probMatrix_r[46];
                    b3 += Ldown[47] * probMatrix_r[47];
                    b4 += Ldown[48] * probMatrix_r[48];
                    b4 += Ldown[49] * probMatrix_r[49];
                    b4 += Ldown[50] * probMatrix_r[50];
                    b4 += Ldown[51] * probMatrix_r[51];
                    b4 += Ldown[52] * probMatrix_r[52];
                    b4 += Ldown[53] * probMatrix_r[53];
                    b4 += Ldown[54] * probMatrix_r[54];
                    b4 += Ldown[55] * probMatrix_r[55];
                    b4 += Ldown[56] * probMatrix_r[56];
                    b4 += Ldown[57] * probMatrix_r[57];
                    b4 += Ldown[58] * probMatrix_r[58];
                    b4 += Ldown[59] * probMatrix_r[59];
                    b4 += Ldown[60] * probMatrix_r[60];
                    b4 += Ldown[61] * probMatrix_r[61];
                    b4 += Ldown[62] * probMatrix_r[62];
                    b4 += Ldown[63] * probMatrix_r[63];
                }
                Lup[i] *= b1 + b2 + b3 + b4;
            }

            ArrayPool.pushStatesByStatesMatrix(Pt);
        }

        @Override
        public double[][] load(final Double branchLength) throws Exception {
            double[][] Pt = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            this.calculator.getTransitionProbabilities(Pt, branchLength);
            return Pt;
        }
    }

    private static final class PtCalculator64 implements TransProbCalculator {
        private final double[] lambda;
        private final double[][] U;
        private final double[][] UInv;
        private final double[][] UInv_cm; // UInv stored in column-major style

        private PtCalculator64(final double[] lambda_in, final double[][] U_in, final double[][] UInv_in, final int[] senseCodons) {
            this.lambda = new double[GeneticCode.CODON_STATES];
            this.U = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            this.UInv = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            this.UInv_cm = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];

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
                }
            }
        }

        @Override
        public final void getTransitionProbabilities(final double[][] Pt, final double branchLength) {
            double p, p1, p2, p3 ,p4 ,p5;
            double[] U_r;
            double[] Pt_r;

            double[] L = new double[GeneticCode.CODON_STATES];
            double[] T = new double[GeneticCode.CODON_STATES];

            // LAMBDA = exp(Λ * branch)
            L[0] = FastMath.exp(branchLength * lambda[0]);
            L[1] = FastMath.exp(branchLength * lambda[1]);
            L[2] = FastMath.exp(branchLength * lambda[2]);
            L[3] = FastMath.exp(branchLength * lambda[3]);
            L[4] = FastMath.exp(branchLength * lambda[4]);
            L[5] = FastMath.exp(branchLength * lambda[5]);
            L[6] = FastMath.exp(branchLength * lambda[6]);
            L[7] = FastMath.exp(branchLength * lambda[7]);
            L[8] = FastMath.exp(branchLength * lambda[8]);
            L[9] = FastMath.exp(branchLength * lambda[9]);
            L[10] = FastMath.exp(branchLength * lambda[10]);
            L[11] = FastMath.exp(branchLength * lambda[11]);
            L[12] = FastMath.exp(branchLength * lambda[12]);
            L[13] = FastMath.exp(branchLength * lambda[13]);
            L[14] = FastMath.exp(branchLength * lambda[14]);
            L[15] = FastMath.exp(branchLength * lambda[15]);
            L[16] = FastMath.exp(branchLength * lambda[16]);
            L[17] = FastMath.exp(branchLength * lambda[17]);
            L[18] = FastMath.exp(branchLength * lambda[18]);
            L[19] = FastMath.exp(branchLength * lambda[19]);
            L[20] = FastMath.exp(branchLength * lambda[20]);
            L[21] = FastMath.exp(branchLength * lambda[21]);
            L[22] = FastMath.exp(branchLength * lambda[22]);
            L[23] = FastMath.exp(branchLength * lambda[23]);
            L[24] = FastMath.exp(branchLength * lambda[24]);
            L[25] = FastMath.exp(branchLength * lambda[25]);
            L[26] = FastMath.exp(branchLength * lambda[26]);
            L[27] = FastMath.exp(branchLength * lambda[27]);
            L[28] = FastMath.exp(branchLength * lambda[28]);
            L[29] = FastMath.exp(branchLength * lambda[29]);
            L[30] = FastMath.exp(branchLength * lambda[30]);
            L[31] = FastMath.exp(branchLength * lambda[31]);
            L[32] = FastMath.exp(branchLength * lambda[32]);
            L[33] = FastMath.exp(branchLength * lambda[33]);
            L[34] = FastMath.exp(branchLength * lambda[34]);
            L[35] = FastMath.exp(branchLength * lambda[35]);
            L[36] = FastMath.exp(branchLength * lambda[36]);
            L[37] = FastMath.exp(branchLength * lambda[37]);
            L[38] = FastMath.exp(branchLength * lambda[38]);
            L[39] = FastMath.exp(branchLength * lambda[39]);
            L[40] = FastMath.exp(branchLength * lambda[40]);
            L[41] = FastMath.exp(branchLength * lambda[41]);
            L[42] = FastMath.exp(branchLength * lambda[42]);
            L[43] = FastMath.exp(branchLength * lambda[43]);
            L[44] = FastMath.exp(branchLength * lambda[44]);
            L[45] = FastMath.exp(branchLength * lambda[45]);
            L[46] = FastMath.exp(branchLength * lambda[46]);
            L[47] = FastMath.exp(branchLength * lambda[47]);
            L[48] = FastMath.exp(branchLength * lambda[48]);
            L[49] = FastMath.exp(branchLength * lambda[49]);
            L[50] = FastMath.exp(branchLength * lambda[50]);
            L[51] = FastMath.exp(branchLength * lambda[51]);
            L[52] = FastMath.exp(branchLength * lambda[52]);
            L[53] = FastMath.exp(branchLength * lambda[53]);
            L[54] = FastMath.exp(branchLength * lambda[54]);
            L[55] = FastMath.exp(branchLength * lambda[55]);
            L[56] = FastMath.exp(branchLength * lambda[56]);
            L[57] = FastMath.exp(branchLength * lambda[57]);
            L[58] = FastMath.exp(branchLength * lambda[58]);
            L[59] = FastMath.exp(branchLength * lambda[59]);
            L[60] = FastMath.exp(branchLength * lambda[60]);
            L[61] = FastMath.exp(branchLength * lambda[61]);
            L[62] = FastMath.exp(branchLength * lambda[62]);
            L[63] = FastMath.exp(branchLength * lambda[63]);

            // U (MATRIX) * Λ (DIAGONAL) = T (MATRIX)
            // T (MATRIX) * U⁻¹ (MATRIX) = Pt
            // We only keep one row of T for each loop
            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                Pt_r = Pt[i];
                U_r = U[i];
                T[0] = U_r[0] * L[0];
                T[1] = U_r[1] * L[1];
                T[2] = U_r[2] * L[2];
                T[3] = U_r[3] * L[3];
                T[4] = U_r[4] * L[4];
                T[5] = U_r[5] * L[5];
                T[6] = U_r[6] * L[6];
                T[7] = U_r[7] * L[7];
                T[8] = U_r[8] * L[8];
                T[9] = U_r[9] * L[9];
                T[10] = U_r[10] * L[10];
                T[11] = U_r[11] * L[11];
                T[12] = U_r[12] * L[12];
                T[13] = U_r[13] * L[13];
                T[14] = U_r[14] * L[14];
                T[15] = U_r[15] * L[15];
                T[16] = U_r[16] * L[16];
                T[17] = U_r[17] * L[17];
                T[18] = U_r[18] * L[18];
                T[19] = U_r[19] * L[19];
                T[20] = U_r[20] * L[20];
                T[21] = U_r[21] * L[21];
                T[22] = U_r[22] * L[22];
                T[23] = U_r[23] * L[23];
                T[24] = U_r[24] * L[24];
                T[25] = U_r[25] * L[25];
                T[26] = U_r[26] * L[26];
                T[27] = U_r[27] * L[27];
                T[28] = U_r[28] * L[28];
                T[29] = U_r[29] * L[29];
                T[30] = U_r[30] * L[30];
                T[31] = U_r[31] * L[31];
                T[32] = U_r[32] * L[32];
                T[33] = U_r[33] * L[33];
                T[34] = U_r[34] * L[34];
                T[35] = U_r[35] * L[35];
                T[36] = U_r[36] * L[36];
                T[37] = U_r[37] * L[37];
                T[38] = U_r[38] * L[38];
                T[39] = U_r[39] * L[39];
                T[40] = U_r[40] * L[40];
                T[41] = U_r[41] * L[41];
                T[42] = U_r[42] * L[42];
                T[43] = U_r[43] * L[43];
                T[44] = U_r[44] * L[44];
                T[45] = U_r[45] * L[45];
                T[46] = U_r[46] * L[46];
                T[47] = U_r[47] * L[47];
                T[48] = U_r[48] * L[48];
                T[49] = U_r[49] * L[49];
                T[50] = U_r[50] * L[50];
                T[51] = U_r[51] * L[51];
                T[52] = U_r[52] * L[52];
                T[53] = U_r[53] * L[53];
                T[54] = U_r[54] * L[54];
                T[55] = U_r[55] * L[55];
                T[56] = U_r[56] * L[56];
                T[57] = U_r[57] * L[57];
                T[58] = U_r[58] * L[58];
                T[59] = U_r[59] * L[59];
                T[60] = U_r[60] * L[60];
                T[61] = U_r[61] * L[61];
                T[62] = U_r[62] * L[62];
                T[63] = U_r[63] * L[63];
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    double[] UInv_cm_r = UInv_cm[j];
                    p = p1 = p2 = p3 = p4 = p5 = 0;
                    p += T[0] * UInv_cm_r[0];
                    p += T[1] * UInv_cm_r[1];
                    p += T[2] * UInv_cm_r[2];
                    p += T[3] * UInv_cm_r[3];
                    p += T[4] * UInv_cm_r[4];
                    p += T[5] * UInv_cm_r[5];
                    p += T[6] * UInv_cm_r[6];
                    p += T[7] * UInv_cm_r[7];
                    p += T[8] * UInv_cm_r[8];
                    p += T[9] * UInv_cm_r[9];
                    p1 += T[10] * UInv_cm_r[10];
                    p1 += T[11] * UInv_cm_r[11];
                    p1 += T[12] * UInv_cm_r[12];
                    p1 += T[13] * UInv_cm_r[13];
                    p1 += T[14] * UInv_cm_r[14];
                    p1 += T[15] * UInv_cm_r[15];
                    p1 += T[16] * UInv_cm_r[16];
                    p1 += T[17] * UInv_cm_r[17];
                    p1 += T[18] * UInv_cm_r[18];
                    p1 += T[19] * UInv_cm_r[19];
                    p2 += T[20] * UInv_cm_r[20];
                    p2 += T[21] * UInv_cm_r[21];
                    p2 += T[22] * UInv_cm_r[22];
                    p2 += T[23] * UInv_cm_r[23];
                    p2 += T[24] * UInv_cm_r[24];
                    p2 += T[25] * UInv_cm_r[25];
                    p2 += T[26] * UInv_cm_r[26];
                    p2 += T[27] * UInv_cm_r[27];
                    p2 += T[28] * UInv_cm_r[28];
                    p2 += T[29] * UInv_cm_r[29];
                    p3 += T[30] * UInv_cm_r[30];
                    p3 += T[31] * UInv_cm_r[31];
                    p3 += T[32] * UInv_cm_r[32];
                    p3 += T[33] * UInv_cm_r[33];
                    p3 += T[34] * UInv_cm_r[34];
                    p3 += T[35] * UInv_cm_r[35];
                    p3 += T[36] * UInv_cm_r[36];
                    p3 += T[37] * UInv_cm_r[37];
                    p3 += T[38] * UInv_cm_r[38];
                    p3 += T[39] * UInv_cm_r[39];
                    p4 += T[40] * UInv_cm_r[40];
                    p4 += T[41] * UInv_cm_r[41];
                    p4 += T[42] * UInv_cm_r[42];
                    p4 += T[43] * UInv_cm_r[43];
                    p4 += T[44] * UInv_cm_r[44];
                    p4 += T[45] * UInv_cm_r[45];
                    p4 += T[46] * UInv_cm_r[46];
                    p4 += T[47] * UInv_cm_r[47];
                    p4 += T[48] * UInv_cm_r[48];
                    p4 += T[49] * UInv_cm_r[49];
                    p5 += T[50] * UInv_cm_r[50];
                    p5 += T[51] * UInv_cm_r[51];
                    p5 += T[52] * UInv_cm_r[52];
                    p5 += T[53] * UInv_cm_r[53];
                    p5 += T[54] * UInv_cm_r[54];
                    p5 += T[55] * UInv_cm_r[55];
                    p5 += T[56] * UInv_cm_r[56];
                    p5 += T[57] * UInv_cm_r[57];
                    p5 += T[58] * UInv_cm_r[58];
                    p5 += T[59] * UInv_cm_r[59];
                    p5 += T[60] * UInv_cm_r[60];
                    p5 += T[61] * UInv_cm_r[61];
                    p5 += T[62] * UInv_cm_r[62];
                    p5 += T[63] * UInv_cm_r[63];

                    Pt_r[j] = p + p1 + p2 + p3 + p4 + p5;
                }
            }
        }

        @Override
        public final void calculatePartialLikelihood(final double[] Ldown, final double[] Lup, final double t) {
            double[] T = new double[GeneticCode.CODON_STATES];

            double p1, p2, p3, p4;

            for (int i = 0; i < 64; i++) {
                double[] UInv_r = UInv[i];
                p1 = p2 = p3 = p4 = 0;

                p1 += UInv_r[0] * Ldown[0];
                p1 += UInv_r[1] * Ldown[1];
                p1 += UInv_r[2] * Ldown[2];
                p1 += UInv_r[3] * Ldown[3];
                p1 += UInv_r[4] * Ldown[4];
                p1 += UInv_r[5] * Ldown[5];
                p1 += UInv_r[6] * Ldown[6];
                p1 += UInv_r[7] * Ldown[7];
                p1 += UInv_r[8] * Ldown[8];
                p1 += UInv_r[9] * Ldown[9];
                p1 += UInv_r[10] * Ldown[10];
                p1 += UInv_r[11] * Ldown[11];
                p1 += UInv_r[12] * Ldown[12];
                p1 += UInv_r[13] * Ldown[13];
                p1 += UInv_r[14] * Ldown[14];
                p1 += UInv_r[15] * Ldown[15];
                p2 += UInv_r[16] * Ldown[16];
                p2 += UInv_r[17] * Ldown[17];
                p2 += UInv_r[18] * Ldown[18];
                p2 += UInv_r[19] * Ldown[19];
                p2 += UInv_r[20] * Ldown[20];
                p2 += UInv_r[21] * Ldown[21];
                p2 += UInv_r[22] * Ldown[22];
                p2 += UInv_r[23] * Ldown[23];
                p2 += UInv_r[24] * Ldown[24];
                p2 += UInv_r[25] * Ldown[25];
                p2 += UInv_r[26] * Ldown[26];
                p2 += UInv_r[27] * Ldown[27];
                p2 += UInv_r[28] * Ldown[28];
                p2 += UInv_r[29] * Ldown[29];
                p3 += UInv_r[30] * Ldown[30];
                p3 += UInv_r[31] * Ldown[31];
                p3 += UInv_r[32] * Ldown[32];
                p3 += UInv_r[33] * Ldown[33];
                p3 += UInv_r[34] * Ldown[34];
                p3 += UInv_r[35] * Ldown[35];
                p3 += UInv_r[36] * Ldown[36];
                p3 += UInv_r[37] * Ldown[37];
                p3 += UInv_r[38] * Ldown[38];
                p3 += UInv_r[39] * Ldown[39];
                p3 += UInv_r[40] * Ldown[40];
                p3 += UInv_r[41] * Ldown[41];
                p3 += UInv_r[42] * Ldown[42];
                p3 += UInv_r[43] * Ldown[43];
                p4 += UInv_r[44] * Ldown[44];
                p4 += UInv_r[45] * Ldown[45];
                p4 += UInv_r[46] * Ldown[46];
                p4 += UInv_r[47] * Ldown[47];
                p4 += UInv_r[48] * Ldown[48];
                p4 += UInv_r[49] * Ldown[49];
                p4 += UInv_r[50] * Ldown[50];
                p4 += UInv_r[51] * Ldown[51];
                p4 += UInv_r[52] * Ldown[52];
                p4 += UInv_r[53] * Ldown[53];
                p4 += UInv_r[54] * Ldown[54];
                p4 += UInv_r[55] * Ldown[55];
                p4 += UInv_r[56] * Ldown[56];
                p4 += UInv_r[57] * Ldown[57];
                p4 += UInv_r[58] * Ldown[58];
                p4 += UInv_r[59] * Ldown[59];
                p4 += UInv_r[60] * Ldown[60];
                p4 += UInv_r[61] * Ldown[61];
                p4 += UInv_r[62] * Ldown[62];
                p4 += UInv_r[63] * Ldown[63];

                T[i] = p1 + p2 + p3 + p4;
            }

            T[0] *= FastMath.exp(lambda[0] * t);
            T[1] *= FastMath.exp(lambda[1] * t);
            T[2] *= FastMath.exp(lambda[2] * t);
            T[3] *= FastMath.exp(lambda[3] * t);
            T[4] *= FastMath.exp(lambda[4] * t);
            T[5] *= FastMath.exp(lambda[5] * t);
            T[6] *= FastMath.exp(lambda[6] * t);
            T[7] *= FastMath.exp(lambda[7] * t);
            T[8] *= FastMath.exp(lambda[8] * t);
            T[9] *= FastMath.exp(lambda[9] * t);
            T[10] *= FastMath.exp(lambda[10] * t);
            T[11] *= FastMath.exp(lambda[11] * t);
            T[12] *= FastMath.exp(lambda[12] * t);
            T[13] *= FastMath.exp(lambda[13] * t);
            T[14] *= FastMath.exp(lambda[14] * t);
            T[15] *= FastMath.exp(lambda[15] * t);
            T[16] *= FastMath.exp(lambda[16] * t);
            T[17] *= FastMath.exp(lambda[17] * t);
            T[18] *= FastMath.exp(lambda[18] * t);
            T[19] *= FastMath.exp(lambda[19] * t);
            T[20] *= FastMath.exp(lambda[20] * t);
            T[21] *= FastMath.exp(lambda[21] * t);
            T[22] *= FastMath.exp(lambda[22] * t);
            T[23] *= FastMath.exp(lambda[23] * t);
            T[24] *= FastMath.exp(lambda[24] * t);
            T[25] *= FastMath.exp(lambda[25] * t);
            T[26] *= FastMath.exp(lambda[26] * t);
            T[27] *= FastMath.exp(lambda[27] * t);
            T[28] *= FastMath.exp(lambda[28] * t);
            T[29] *= FastMath.exp(lambda[29] * t);
            T[30] *= FastMath.exp(lambda[30] * t);
            T[31] *= FastMath.exp(lambda[31] * t);
            T[32] *= FastMath.exp(lambda[32] * t);
            T[33] *= FastMath.exp(lambda[33] * t);
            T[34] *= FastMath.exp(lambda[34] * t);
            T[35] *= FastMath.exp(lambda[35] * t);
            T[36] *= FastMath.exp(lambda[36] * t);
            T[37] *= FastMath.exp(lambda[37] * t);
            T[38] *= FastMath.exp(lambda[38] * t);
            T[39] *= FastMath.exp(lambda[39] * t);
            T[40] *= FastMath.exp(lambda[40] * t);
            T[41] *= FastMath.exp(lambda[41] * t);
            T[42] *= FastMath.exp(lambda[42] * t);
            T[43] *= FastMath.exp(lambda[43] * t);
            T[44] *= FastMath.exp(lambda[44] * t);
            T[45] *= FastMath.exp(lambda[45] * t);
            T[46] *= FastMath.exp(lambda[46] * t);
            T[47] *= FastMath.exp(lambda[47] * t);
            T[48] *= FastMath.exp(lambda[48] * t);
            T[49] *= FastMath.exp(lambda[49] * t);
            T[50] *= FastMath.exp(lambda[50] * t);
            T[51] *= FastMath.exp(lambda[51] * t);
            T[52] *= FastMath.exp(lambda[52] * t);
            T[53] *= FastMath.exp(lambda[53] * t);
            T[54] *= FastMath.exp(lambda[54] * t);
            T[55] *= FastMath.exp(lambda[55] * t);
            T[56] *= FastMath.exp(lambda[56] * t);
            T[57] *= FastMath.exp(lambda[57] * t);
            T[58] *= FastMath.exp(lambda[58] * t);
            T[59] *= FastMath.exp(lambda[59] * t);
            T[60] *= FastMath.exp(lambda[60] * t);
            T[61] *= FastMath.exp(lambda[61] * t);
            T[62] *= FastMath.exp(lambda[62] * t);
            T[63] *= FastMath.exp(lambda[63] * t);

            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                double[] U_c = U[i];
                p1 = p2 = p3 = p4 = 0;

                p1 += U_c[0] * T[0];
                p1 += U_c[1] * T[1];
                p1 += U_c[2] * T[2];
                p1 += U_c[3] * T[3];
                p1 += U_c[4] * T[4];
                p1 += U_c[5] * T[5];
                p1 += U_c[6] * T[6];
                p1 += U_c[7] * T[7];
                p1 += U_c[8] * T[8];
                p1 += U_c[9] * T[9];
                p1 += U_c[10] * T[10];
                p2 += U_c[11] * T[11];
                p2 += U_c[12] * T[12];
                p2 += U_c[13] * T[13];
                p2 += U_c[14] * T[14];
                p2 += U_c[15] * T[15];
                p2 += U_c[16] * T[16];
                p2 += U_c[17] * T[17];
                p2 += U_c[18] * T[18];
                p2 += U_c[19] * T[19];
                p2 += U_c[20] * T[20];
                p2 += U_c[21] * T[21];
                p2 += U_c[22] * T[22];
                p2 += U_c[23] * T[23];
                p2 += U_c[24] * T[24];
                p2 += U_c[25] * T[25];
                p2 += U_c[26] * T[26];
                p3 += U_c[27] * T[27];
                p3 += U_c[28] * T[28];
                p3 += U_c[29] * T[29];
                p3 += U_c[30] * T[30];
                p3 += U_c[31] * T[31];
                p3 += U_c[32] * T[32];
                p3 += U_c[33] * T[33];
                p3 += U_c[34] * T[34];
                p3 += U_c[35] * T[35];
                p3 += U_c[36] * T[36];
                p3 += U_c[37] * T[37];
                p3 += U_c[38] * T[38];
                p3 += U_c[39] * T[39];
                p3 += U_c[40] * T[40];
                p3 += U_c[41] * T[41];
                p3 += U_c[42] * T[42];
                p3 += U_c[43] * T[43];
                p3 += U_c[44] * T[44];
                p4 += U_c[45] * T[45];
                p4 += U_c[46] * T[46];
                p4 += U_c[47] * T[47];
                p4 += U_c[48] * T[48];
                p4 += U_c[49] * T[49];
                p4 += U_c[50] * T[50];
                p4 += U_c[51] * T[51];
                p4 += U_c[52] * T[52];
                p4 += U_c[53] * T[53];
                p4 += U_c[54] * T[54];
                p4 += U_c[55] * T[55];
                p4 += U_c[56] * T[56];
                p4 += U_c[57] * T[57];
                p4 += U_c[58] * T[58];
                p4 += U_c[59] * T[59];
                p4 += U_c[60] * T[60];
                p4 += U_c[61] * T[61];
                p4 += U_c[62] * T[62];
                p4 += U_c[63] * T[63];


                Lup[i] *= p1 + p2 + p3 + p4;
            }

/*            for (int i = 0; i < DIM; i++) {
                double[] UInv_r = UInv[i];
                double p = 0;
                for (int j = 0; j < DIM; j++) {
                    p += UInv_r[j] * Ldown[senseCodons[j]];
                }
                temp[i] = p;
            }

            for (int i = 0; i < DIM; i++) {
                temp[i] *= FastMath.exp(lambda[i] * t);
            }

            for (int i = 0; i < DIM; i++) {
                double p = 0;
                for (int j = 0; j < DIM; j++) {
                    p += U[i][j] * temp[j];
                }
                if (p < 0) p = 0;
                Lup[senseCodons[i]] *= p;
            }*/
        }
    }
}

