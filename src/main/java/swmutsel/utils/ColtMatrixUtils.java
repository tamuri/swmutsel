package swmutsel.utils;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.*;
import swmutsel.colt.Algebra;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class ColtMatrixUtils {

    public static DoubleMatrix1D make(double... doubles) {
        return DoubleFactory1D.dense.make(doubles);
    }

    public static DoubleMatrix2D make(int rows, double... doubles) {
        return DoubleFactory2D.dense.make(doubles, rows);
    }

    // From c++ source:
    // http://code.google.com/p/numerical-recipes-java/source/browse/trunk/numerical-recipes-j/core/src/main/cpp/eigen_unsym.h
    public static void balance(DoubleMatrix2D a) {
        double RADIX = 2;
        boolean done = false;
        double RADIX_SQ = RADIX * RADIX;

        int n = a.rows();

        while (!done) {
            done = true;
            for (int i = 0; i < n; i++) {
                double r = 0.0, c = 0.0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        c += Math.abs(a.getQuick(j, i));
                        r += Math.abs(a.getQuick(i, j));
                    }
                }
                if (c != 0.0 && r != 0.0) {
                    double g = r / RADIX;
                    double f = 1.0;
                    double s = c + r;
                    while (c < g) {
                        f *= RADIX;
                        c *= RADIX_SQ;
                    }
                    g = r * RADIX;
                    while (c > g) {
                        f /= RADIX;
                        c /= RADIX_SQ;
                    }
                    if ((c + r) / f < 0.95 * s) {
                        done = false;
                        g = 1.0 / f;
                        // scale[i] *= f;
                        for (int j = 0; j < n; j++)
                            a.setQuick(i, j, a.getQuick(i, j) * g);
                        for (int j = 0; j < n; j++)
                            a.setQuick(j, i, a.getQuick(j, i) * f);
                    }
                }
            }
        }

    }

    private static DoubleMatrix2D T = DoubleFactory2D.dense.make(64, 64);
    private static DoubleMatrix2D B = DoubleFactory2D.dense.make(64, 64);

    public static DoubleMatrix2D taylorExpTerm(DoubleMatrix2D Q, int k, DoubleMatrix2D out) {
        T.assign(0);
        B.assign(0);
        //Algebra.DEFAULT.pow(Q, k, out, T, B); <--- WHY DOESN'T THIS WORK??
        out.assign(Algebra.DEFAULT.pow(Q, k, out, T, B));
        out.assign(cern.jet.math.Functions.div(Arithmetic.factorial(k)));
        return out;
    }




}
