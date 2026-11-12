package swmutsel.model;

import org.junit.Before;
import org.junit.Test;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessFDS;
import swmutsel.utils.GeneticCode;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import static org.junit.Assert.assertEquals;

/**
 * Each PadeApproximantCalculator must depend only on its own rate matrix, whatever other Padé calculators exist,
 * and must give the same results when used from several threads at once.
 *
 * Expected values are computed with each calculator on its own, before any other calculator is created, so the
 * comparisons are exact (same arithmetic, same inputs).
 */
public class PadeApproximantCalculatorTest {
    private static final double SCALING = 1.5;
    private static final double[] FITNESS = {0.0, 0.5, -1, 1, 0, 0.2, 0.3, -0.5, 0, 1.5, 0.7, -1, 0.1, 0.4, 2, 0, -0.2, 0.3, 0.9, 1.1};
    private static final double[] BRANCH_LENGTHS = {0.01, 0.1, 0.5, 1, 2, 5};

    private static final int THREADS = 4;
    private static final int REPEATS = 200;

    private int[] sense;
    private double[] ldown;

    @Before
    public void setUp() {
        GeneticCode.setCode(GeneticCode.STANDARD_CODE);
        sense = GeneticCode.getInstance().getSenseCodons();

        // An arbitrary conditional likelihood vector for the partial likelihood calculation
        ldown = new double[GeneticCode.CODON_STATES];
        for (int i : sense) ldown[i] = 0.01 + (i % 7) / 7.0;
    }

    @Test
    public void creatingASecondCalculatorDoesNotChangeTheFirst() {
        TransProbCalculator a = TransProbCalcFactory.getPadeCalculator(fdsQ(0.5), sense, SCALING);
        Results expectedA = compute(a);

        TransProbCalculator b = TransProbCalcFactory.getPadeCalculator(fdsQ(10.0), sense, SCALING);
        compute(b);

        assertSame("calculator A after creating and using B", expectedA, compute(a));
    }

    @Test
    public void differentCalculatorsInDifferentThreads() throws Exception {
        TransProbCalculator a = TransProbCalcFactory.getPadeCalculator(fdsQ(0.5), sense, SCALING);
        Results expectedA = compute(a);
        TransProbCalculator b = TransProbCalcFactory.getPadeCalculator(fdsQ(10.0), sense, SCALING);
        Results expectedB = compute(b);

        List<TransProbCalculator> calculators = new ArrayList<TransProbCalculator>();
        List<Results> expected = new ArrayList<Results>();
        for (int i = 0; i < THREADS; i++) {
            calculators.add(i % 2 == 0 ? a : b);
            expected.add(i % 2 == 0 ? expectedA : expectedB);
        }

        assertEquals("mismatched results", 0, runConcurrently(calculators, expected));
    }

    @Test
    public void oneCalculatorSharedByThreads() throws Exception {
        TransProbCalculator a = TransProbCalcFactory.getPadeCalculator(fdsQ(2.0), sense, SCALING);
        Results expectedA = compute(a);

        List<TransProbCalculator> calculators = new ArrayList<TransProbCalculator>();
        List<Results> expected = new ArrayList<Results>();
        for (int i = 0; i < THREADS; i++) {
            calculators.add(a);
            expected.add(expectedA);
        }

        assertEquals("mismatched results", 0, runConcurrently(calculators, expected));
    }

    /** Runs calculator i in thread i, REPEATS times, and counts results that differ from expected i. */
    private int runConcurrently(final List<TransProbCalculator> calculators, final List<Results> expected) throws Exception {
        ExecutorService pool = Executors.newFixedThreadPool(calculators.size());
        final CountDownLatch start = new CountDownLatch(1);
        List<Future<Integer>> futures = new ArrayList<Future<Integer>>();

        try {
            for (int i = 0; i < calculators.size(); i++) {
                final TransProbCalculator calculator = calculators.get(i);
                final Results want = expected.get(i);
                futures.add(pool.submit(new Callable<Integer>() {
                    @Override
                    public Integer call() throws Exception {
                        start.await();
                        int mismatches = 0;
                        for (int r = 0; r < REPEATS; r++) {
                            if (!compute(calculator).equalTo(want)) mismatches++;
                        }
                        return mismatches;
                    }
                }));
            }

            start.countDown();

            int mismatches = 0;
            for (Future<Integer> f : futures) mismatches += f.get();
            return mismatches;
        } finally {
            pool.shutdownNow();
        }
    }

    /** P(t) and a partial likelihood, for every branch length. */
    private Results compute(TransProbCalculator calculator) {
        Results r = new Results();
        for (int k = 0; k < BRANCH_LENGTHS.length; k++) {
            double[][] P = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
            calculator.getTransitionProbabilities(P, BRANCH_LENGTHS[k]);
            r.pt[k] = P;

            double[] lup = new double[GeneticCode.CODON_STATES];
            java.util.Arrays.fill(lup, 1.0);
            calculator.calculatePartialLikelihood(ldown, lup, BRANCH_LENGTHS[k]);
            r.partial[k] = lup;
        }
        return r;
    }

    private double[] fdsQ(double z) {
        return new SwMutSel(new SwMut(2.0, SCALING, new double[]{0.3, 0.2, 0.25, 0.25}),
                new Fitness(FITNESS.clone()), new FitnessFDS(z, 1)).getQ();
    }

    private static void assertSame(String message, Results expected, Results actual) {
        if (!actual.equalTo(expected)) {
            throw new AssertionError(message + ": results differ (max difference " + actual.maxDiff(expected) + ")");
        }
    }

    private static final class Results {
        final double[][][] pt = new double[BRANCH_LENGTHS.length][][];
        final double[][] partial = new double[BRANCH_LENGTHS.length][];

        boolean equalTo(Results other) {
            return maxDiff(other) == 0;
        }

        double maxDiff(Results other) {
            double max = 0;
            for (int k = 0; k < BRANCH_LENGTHS.length; k++) {
                for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                    for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                        max = Math.max(max, Math.abs(pt[k][i][j] - other.pt[k][i][j]));
                    }
                    max = Math.max(max, Math.abs(partial[k][i] - other.partial[k][i]));
                }
            }
            return max;
        }
    }
}
