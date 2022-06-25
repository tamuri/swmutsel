package riso.numerical;

import com.google.common.collect.Table;
import org.apache.commons.math3.analysis.MultivariateFunction;
import pal.alignment.Alignment;
import pal.tree.Tree;
import swmutsel.ArrayPool;
import swmutsel.Constants;
import swmutsel.cli.ArgumentsProcessor;
import swmutsel.model.FMutSel0;
import swmutsel.model.SubstitutionModel;
import swmutsel.model.parameters.*;
import swmutsel.runner.MultiThreadedRunner;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.MathUtils;
import swmutsel.utils.PhyloUtils;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 12/03/2014 12:36
 */
public class TestLBFGS {
    static int ndim = 2000 , msave = 7 ;
    static int nwork = ndim * ( 2 * msave + 1 ) + 2 * msave ;

    public static void main(String[] args) throws Exception {
        TestLBFGS t = new TestLBFGS();
        t.run2();
    }

    private void run2() throws Exception {
        final SubstitutionModel model;
        final MultiThreadedRunner r;
        final Tree t;

        t = PhyloUtils.readTree("/Users/atamuri/Documents/2014/mutselbranch/data/small4s/small4s_mtCDNA_500.tree");
        ArrayPool.setTreeSize(t.getExternalNodeCount());
        Alignment a = ArgumentsProcessor.loadPalAlignment("/Users/atamuri/Documents/2014/mutselbranch/data/small4s/small4s_mtCDNA_100.txt");

        Table<String, Integer, Byte> table = ArgumentsProcessor.getAlignmentTable(a);

        r = new MultiThreadedRunner(3);
        r.setSites(table);

        TsTvRatio k = new TsTvRatio(1);
        Omega w = new Omega(1);
        BaseFrequencies pi = new BaseFrequencies(new double[]{1, 1, 1, 1});
        Fitness F = new Fitness(
                new int[]{7, 0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19},
                new double[20]
        );

        model = new FMutSel0(k, w, pi, F);
        model.build();

        double[] x = Mapper.getOptimisable(model.getParameters());
        int n = x.length;

        int m = 5; // no. of corrections
        double[] diag = new double[n];
        double eps, xtol, gtol, t1, t2, stpmin, stpmax;
        int iprint [ ] , iflag[] = new int[1], icall,  mp, lp, j;
        iprint = new int [ 2 ];
        boolean diagco;
        iprint [ 1 -1] = 1;
        iprint [ 2 -1] = 0;
        diagco= false;
        eps= Constants.VALUE_CONVERGENCE_TOL;
        xtol= 1.0e-16;
        icall=0;
        iflag[0]=0;


        r.setRunnerTree(t);

        do {

            Mapper.setOptimisable(model.getParameters(), x);
           // System.out.printf("kappa=%s\nomega=%s\npi=%s\nF=%s\n", k.get(), w.get(), Doubles.join(",", pi.get()), Doubles.join(",", F.get()));
            model.build();

            double f = -CoreUtils.sum(r.getLogLikelihood(model).values());
            //System.out.printf("f=%s\n\n", f);


            double[] g = MathUtils.gradient(new MultivariateFunction() {
                @Override
                public double value(double[] point) {
                    Mapper.setOptimisable(model.getParameters(), point);
                    model.build();
                    return CoreUtils.sum(r.getLogLikelihood(model).values());
                }
            }, x);


            LBFGS.lbfgs(n, m, x, f, g, diagco, diag, iprint, eps, xtol, iflag);


            icall++;

        } while (iflag[0] != 0 && icall <= 500 );


        r.shutdown();


    }

    private void run() {
        double x [ ] , g [ ] , diag [ ] , w [ ];
        x = new double [ ndim ]; // first guess, then filled with best estimate
        g = new double [ ndim ]; // gradient at point x
        diag = new double [ ndim ];
        w = new double [ nwork ];

        double f, eps, xtol, gtol, t1, t2, stpmin, stpmax;
        int iprint [ ] , iflag[] = new int[1], icall, n, m, mp, lp, j;
        iprint = new int [ 2 ]; // output level
        boolean diagco;

        n=100; // number of parameters in the problem
        m=5; // number of corrections
        iprint [ 1 -1] = 1;
        iprint [ 2 -1] = 0;
        diagco= false;
        eps= 1.0e-5;
        xtol= 1.0e-16;
        icall=0;
        iflag[0]=0;

        for ( j = 1 ; j <= n ; j += 2 )
        {
            x [ j -1] = - 1.2e0;
            x [ j + 1 -1] = 1.e0;
        }

        do
        {
            f= 0;
            for ( j = 1 ; j <= n ; j += 2 )
            {
                t1 = 1.e0 - x [ j -1];
                t2 = 1.e1 * ( x [ j + 1 -1] - x [ j -1] * x[j-1] );
                g [ j + 1 -1] = 2.e1 * t2;
                g [ j -1] = - 2.e0 * ( x [ j -1] * g [ j + 1 -1] + t1 );
                f= f+t1*t1+t2*t2;
            }

            try
            {
                LBFGS.lbfgs ( n , m , x , f , g , diagco , diag , iprint , eps , xtol , iflag );
            }
            catch (LBFGS.ExceptionWithIflag e)
            {
                System.err.println( "Sdrive: lbfgs failed.\n"+e );
                return;
            }

            icall += 1;
        }
        while ( iflag[0] != 0 && icall <= 200 );

    }
}
