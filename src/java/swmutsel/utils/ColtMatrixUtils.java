package swmutsel.utils;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

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
}
