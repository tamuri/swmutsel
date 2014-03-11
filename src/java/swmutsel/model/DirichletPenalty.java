package swmutsel.model;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.math.Functions;

public class DirichletPenalty implements Penalty {

    private static final long serialVersionUID = -7423296435995772451L;

    public static final String LABEL = "dirichlet";

    private static final int DIM = 19;
    private final double alpha;

    public DirichletPenalty(double alpha) {
        this.alpha = alpha;
    }

    @Override
    public double calculate(final double[] fitness) {
        DoubleMatrix1D theta = DoubleFactory1D.dense.make(fitness);
        theta.assign(Functions.exp);
        theta.assign(Functions.div(1 + theta.zSum()));


        // original:
        return ((alpha - 1) * theta.aggregate(Functions.plus, Functions.log))
                + ((alpha - 1) * Math.log(1 - theta.zSum()))
                + Math.log(Math.abs(Algebra.DEFAULT.det(jacobian(theta))));


        // MdR's function in the manuscript:
        /*
        DoubleMatrix1D theta = DoubleFactory1D.dense.make(20, 0.0);
        for (int i = 0; i < fitness.length; i++) theta.setQuick(i, fitness[i]);

        theta.assign(Functions.exp);
        theta.assign(Functions.div(theta.zSum()));

        return alpha * theta.aggregate(Functions.plus, Functions.log);
        */


    }

    private DoubleMatrix2D jacobian(final DoubleMatrix1D theta) {
        DoubleMatrix2D jacobi = DoubleFactory2D.dense.make(DIM, DIM);
        for (int i = 0; i < DIM; i++) {
            for (int j = 0; j < DIM; j++) {
                if (i == j) {
                    jacobi.setQuick(i, j, theta.getQuick(i) * (1 - theta.getQuick(i)));
                } else {
                    jacobi.setQuick(i, j, -theta.getQuick(i) * theta.getQuick(j));
                }
            }
        }
        return jacobi;
    }

    @Override
    public String toString() {
        return LABEL + "," + alpha;
    }


  /*  public static void main(String[] args) {

        double[] f = new double[19];
        for (int i = 0; i < f.length; i++) {
            f[i] = Normal.staticNextDouble(0, 1);
        }

        f = new double[]{-1.0135940331511981, 0.9254958707503568, -0.29255612654565816, -0.22681279848729385, 0.47387836055700755, -0.15291969247995885, -0.8747788652108226, -1.7177774513544954, 1.3341023647790993, -1.264047307128678, -3.129557842977072, 0.10073966555682189, 1.3617232021873107, 1.1981880441978687, 0.44439692707114636, -0.5022015366261262, 1.676426029781541, 0.020532038066831744, 0.008691116712494517};

        System.out.printf("Fitness: %s\n", Doubles.join(", ", f));

        for (double a : new double[]{1.0, 0.5, 0.1, 0.01, 0.0001}) {
            DirichletPenalty dp = new DirichletPenalty(a);
            double cal1 = dp.calculate(f);

            System.out.printf("%s\t%s\n", a, cal1);
        }

        System.out.printf("length: %s\n", f.length);

        for (double a : new double[]{1.0, 0.1, 0.01, 0.0001}) {
            double sum1 = 0;
            for (double ff : f) {
                sum1 += Math.exp(ff);
            }
            sum1 += 1;
            sum1 = Math.log(sum1);
            System.out.printf("sum1 = %s\n", sum1);

            double sum2 = 0;
            for (double ff : f) {
                sum2 += ff - sum1;
            }
            System.out.printf("sum2 = %s\n", sum2);

            double p = sum2 * a;
            System.out.printf("%s\t%s\n", a, p);
        }
    }
*/

}
