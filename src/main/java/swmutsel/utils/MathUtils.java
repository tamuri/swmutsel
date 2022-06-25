package swmutsel.utils;

import cern.colt.matrix.DoubleMatrix2D;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.analysis.MultivariateFunction;
import pal.math.MachineAccuracy;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 12/03/2014 10:57
 */
public class MathUtils {
    public static double[] gradient(MultivariateFunction f, double[] x) {
        double[] grad = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            double h = MachineAccuracy.SQRT_EPSILON * (Math.abs(x[i]) + 1.0);

            /* codeml step size: double h = 1e-6 * 2 * Math.abs(x[i] + 1.); */

            double oldx = x[i];
            x[i] = oldx + h;
            double fxplus = -f.value(x); // notice this is -ve! For L-BFGS-S [minimising] optimiser
            x[i] = oldx - h;
            double fxminus = -f.value(x); // notice this is -ve!
            x[i] = oldx;

            // Centered first derivative
            grad[i] = (fxplus - fxminus) / (2.0 * h);
        }
        return grad;
    }

    public static double sumLogs(double x, double y) {
        if (x > y)
            return x + Math.log(1 + Math.exp(y - x));
        else
            return y + Math.log(1 + Math.exp(x - y));
    }

    public static double sumLogs(double[] a) {
        double sum = a[0];
        for (int i = 1; i < a.length; i++) {
            sum = sumLogs(sum, a[i]);
        }
        return sum;
    }

    public static double[] smartRound(double[] in, int digits) {
        /* Based on the following R code, courtesy of Mikhail at http://stackoverflow.com/a/35930285
        smart.round <- function(x, digits = 0) {
          up <- 10 ^ digits
          x <- x * up
          y <- floor(x)
          indices <- tail(order(x-y), round(sum(x)) - sum(y))
          y[indices] <- y[indices] + 1
          y / up
        }
        */

        double up = Math.pow(10, digits);

        double[] x = new double[in.length];
        double[] y = new double[in.length];
        double[] decimals = new double[in.length];
        double xSum = 0;
        double ySum = 0;

        for (int i = 0; i < in.length; i++) {
            x[i] = in[i] * up;
            y[i] = Math.floor(x[i]);
            decimals[i] = x[i] - y[i];
            xSum += x[i];
            ySum += y[i];
        }

        int[] sortedIndices = IntStream.range(0, y.length).boxed().sorted(Comparator.comparingDouble(i -> decimals[i])).mapToInt(i -> i).toArray();
        int places = (int) Math.round(xSum - ySum);

        for (int j = y.length - places; j < y.length; j++) {
            int i = sortedIndices[j];
            y[i] = y[i] + 1;
        }

        for (int i = 0; i < y.length; i++) {
            y[i] = y[i] / up;
        }

        return y;
    }

    private static void test_smartRound() {
        // Sum using the exact values
        // double[] x = new double[]{0.13626332, 0.47989636, 0.09596008, 0.28788024};
        double[] x = new double[]{0.1913489,0.6738984,0.1347527};
        System.out.println(Doubles.join(", ", x));
        System.out.println(Arrays.stream(x).sum());
        System.out.println();

        // Standard rounds to 4 decimal places, and sum
        double[] xx = Arrays.stream(x).map(i -> round(i, 2)).toArray();
        System.out.println(Doubles.join(", ", xx));
        System.out.println(Arrays.stream(xx).sum());
        System.out.println();

        // Smart rounding procedure
        double[] y = smartRound(x, 2);
        System.out.println(Doubles.join(", ", y));
        System.out.println(Arrays.stream(y).sum());
    }

    public static void main(String[] args) {
        test_smartRound();
    }

    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();

        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_EVEN);
        return bd.doubleValue();
    }

    /**
     * COPIED FROM cern.colt.matrix.linalg.EigenvalueDecomposition
     *
     * This version:
     * 1) Assumes that the matrix to diagonalize is symmetric
     * 2) Populates values in arrays passed as arguments
     */
    public static void setEigenFactors(final DoubleMatrix2D matrix, final double[] diagEigenValue, final double[][] eigenVectorMatrix, final int size) {
        double[] e = new double[size];

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                eigenVectorMatrix[i][j] = matrix.getQuick(i,j);
            }
        }

        // Tridiagonalize.
        tred2(size, diagEigenValue, eigenVectorMatrix, e);

        // Diagonalize.
        tql2(size, diagEigenValue, eigenVectorMatrix, e);
    }

    /**
     * COPIED FROM cern.colt.matrix.linalg.EigenvalueDecomposition
     */
    private static void tql2 (int n, double[] d, double[][] V, double[] e) {

        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.

        for (int i = 1; i < n; i++) {
            e[i-1] = e[i];
        }
        e[n-1] = 0.0;

        double f = 0.0;
        double tst1 = 0.0;
        double eps = Math.pow(2.0,-52.0);
        for (int l = 0; l < n; l++) {

            // Find small subdiagonal element

            tst1 = Math.max(tst1,Math.abs(d[l]) + Math.abs(e[l]));
            int m = l;
            while (m < n) {
                if (Math.abs(e[m]) <= eps*tst1) {
                    break;
                }
                m++;
            }

            // If m == l, d[l] is an eigenvalue,
            // otherwise, iterate.

            if (m > l) {
                int iter = 0;
                do {
                    iter = iter + 1;  // (Could check iteration count here.)

                    // Compute implicit shift

                    double g = d[l];
                    double p = (d[l+1] - g) / (2.0 * e[l]);
                    double r = hypot(p, 1.0);
                    if (p < 0) {
                        r = -r;
                    }
                    d[l] = e[l] / (p + r);
                    d[l+1] = e[l] * (p + r);
                    double dl1 = d[l+1];
                    double h = g - d[l];
                    for (int i = l+2; i < n; i++) {
                        d[i] -= h;
                    }
                    f = f + h;

                    // Implicit QL transformation.

                    p = d[m];
                    double c = 1.0;
                    double c2 = c;
                    double c3 = c;
                    double el1 = e[l+1];
                    double s = 0.0;
                    double s2 = 0.0;
                    for (int i = m-1; i >= l; i--) {
                        c3 = c2;
                        c2 = c;
                        s2 = s;
                        g = c * e[i];
                        h = c * p;
                        r = hypot(p,e[i]);
                        e[i+1] = s * r;
                        s = e[i] / r;
                        c = p / r;
                        p = c * d[i] - s * g;
                        d[i+1] = h + s * (c * g + s * d[i]);

                        // Accumulate transformation.

                        for (int k = 0; k < n; k++) {
                            h = V[k][i+1];
                            V[k][i+1] = s * V[k][i] + c * h;
                            V[k][i] = c * V[k][i] - s * h;
                        }
                    }
                    p = -s * s2 * c3 * el1 * e[l] / dl1;
                    e[l] = s * p;
                    d[l] = c * p;

                    // Check for convergence.

                } while (Math.abs(e[l]) > eps*tst1);
            }
            d[l] = d[l] + f;
            e[l] = 0.0;
        }

        // Sort eigenvalues and corresponding vectors.

        for (int i = 0; i < n-1; i++) {
            int k = i;
            double p = d[i];
            for (int j = i+1; j < n; j++) {
                if (d[j] < p) {
                    k = j;
                    p = d[j];
                }
            }
            if (k != i) {
                d[k] = d[i];
                d[i] = p;
                for (int j = 0; j < n; j++) {
                    p = V[j][i];
                    V[j][i] = V[j][k];
                    V[j][k] = p;
                }
            }
        }
    }

    /**
     * COPIED FROM cern.colt.matrix.linalg.EigenvalueDecomposition
     */
    private static void tred2 (int n, double[] d, double[][] V, double[] e) {
        //  This is derived from the Algol procedures tred2 by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.


        for (int j = 0; j < n; j++) {
            d[j] = V[n-1][j];
        }


        // Householder reduction to tridiagonal form.

        for (int i = n-1; i > 0; i--) {

            // Scale to avoid under/overflow.

            double scale = 0.0;
            double h = 0.0;
            for (int k = 0; k < i; k++) {
                scale = scale + Math.abs(d[k]);
            }
            if (scale == 0.0) {
                e[i] = d[i-1];
                for (int j = 0; j < i; j++) {
                    d[j] = V[i-1][j];
                    V[i][j] = 0.0;
                    V[j][i] = 0.0;
                }
            } else {

                // Generate Householder vector.

                for (int k = 0; k < i; k++) {
                    d[k] /= scale;
                    h += d[k] * d[k];
                }
                double f = d[i-1];
                double g = Math.sqrt(h);
                if (f > 0) {
                    g = -g;
                }
                e[i] = scale * g;
                h = h - f * g;
                d[i-1] = f - g;
                for (int j = 0; j < i; j++) {
                    e[j] = 0.0;
                }

                // Apply similarity transformation to remaining columns.

                for (int j = 0; j < i; j++) {
                    f = d[j];
                    V[j][i] = f;
                    g = e[j] + V[j][j] * f;
                    for (int k = j+1; k <= i-1; k++) {
                        g += V[k][j] * d[k];
                        e[k] += V[k][j] * f;
                    }
                    e[j] = g;
                }
                f = 0.0;
                for (int j = 0; j < i; j++) {
                    e[j] /= h;
                    f += e[j] * d[j];
                }
                double hh = f / (h + h);
                for (int j = 0; j < i; j++) {
                    e[j] -= hh * d[j];
                }
                for (int j = 0; j < i; j++) {
                    f = d[j];
                    g = e[j];
                    for (int k = j; k <= i-1; k++) {
                        V[k][j] -= (f * e[k] + g * d[k]);
                    }
                    d[j] = V[i-1][j];
                    V[i][j] = 0.0;
                }
            }
            d[i] = h;
        }

        // Accumulate transformations.

        for (int i = 0; i < n-1; i++) {
            V[n-1][i] = V[i][i];
            V[i][i] = 1.0;
            double h = d[i+1];
            if (h != 0.0) {
                for (int k = 0; k <= i; k++) {
                    d[k] = V[k][i+1] / h;
                }
                for (int j = 0; j <= i; j++) {
                    double g = 0.0;
                    for (int k = 0; k <= i; k++) {
                        g += V[k][i+1] * V[k][j];
                    }
                    for (int k = 0; k <= i; k++) {
                        V[k][j] -= g * d[k];
                    }
                }
            }
            for (int k = 0; k <= i; k++) {
                V[k][i+1] = 0.0;
            }
        }
        for (int j = 0; j < n; j++) {
            d[j] = V[n-1][j];
            V[n-1][j] = 0.0;
        }
        V[n-1][n-1] = 1.0;
        e[0] = 0.0;
    }

    /**
     * COPIED FROM cern.colt.matrix.linalg.Algebra
     */
    private static double hypot(double a, double b) {
        double r;
        if (Math.abs(a) > Math.abs(b)) {
            r = b/a;
            r = Math.abs(a)*Math.sqrt(1+r*r);
        } else if (b != 0) {
            r = a/b;
            r = Math.abs(b)*Math.sqrt(1+r*r);
        } else {
            r = 0.0;
        }
        return r;
    }


}
