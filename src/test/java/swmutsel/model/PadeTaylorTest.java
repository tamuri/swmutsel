package swmutsel.model;

import org.junit.Before;
import org.junit.Test;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessFDS;
import swmutsel.utils.GeneticCode;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * The -fds model has a non-reversible rate matrix, so its transition probabilities P(t) = exp(Q c t) are computed
 * with a Padé approximant (PadeApproximantCalculator) instead of an eigendecomposition. TaylorExpansionCalculator
 * is an independent way of computing the same thing; it is kept only for this test.
 *
 * 1. For a reversible matrix (no FDS), both are checked against the eigendecomposition, which is accurate.
 * 2. For FDS matrices, Padé is checked against Taylor (within the range where Taylor is valid).
 * 3. Padé's P(t) are checked to be probabilities, including beyond Taylor's range.
 *
 * Comparisons use the sense-codon block only: the calculators fill the stop-codon rows differently.
 *
 * NOTE: PadeApproximantCalculator keeps its matrices in static fields, so creating a second Padé calculator
 * overwrites the first. Each Padé calculator here is used before the next one is created.
 */
public class PadeTaylorTest {
    private static final double SCALING = 1.5;
    private static final double[] PI = {0.3, 0.2, 0.25, 0.25};
    private static final double KAPPA = 2.0;
    private static final double[] FITNESS = {0.0, 0.5, -1, 1, 0, 0.2, 0.3, -0.5, 0, 1.5, 0.7, -1, 0.1, 0.4, 2, 0, -0.2, 0.3, 0.9, 1.1};

    private static final double[] BRANCH_LENGTHS = {0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10};
    private static final double[] FDS_Z = {0.5, 2.0, 10.0};

    // Taylor diverges when ||Q c t||_1 is large (see TaylorExpansionCalculator); only compare below this
    private static final double TAYLOR_MAX_NORM = 100;

    // Observed maximum differences (sense block): Padé vs exact 6.3e-8, Taylor vs exact 7.0e-7, Padé vs Taylor 6.8e-7
    private static final double PADE_VS_EXACT_TOLERANCE = 2e-7;
    private static final double TAYLOR_TOLERANCE = 2e-6;
    private static final double PROBABILITY_TOLERANCE = 1e-12;

    private int[] sense;

    @Before
    public void setUp() {
        GeneticCode.setCode(GeneticCode.STANDARD_CODE);
        sense = GeneticCode.getInstance().getSenseCodons();
    }

    @Test
    public void reversibleModelPadeAndTaylorMatchEigendecomposition() {
        SwMutSel model = new SwMutSel(mutation(), new Fitness(FITNESS.clone()));
        double[] Q = model.getQ();

        TransProbCalculator exact = model.getPtCalculator(); // eigendecomposition for a reversible model
        TransProbCalculator taylor = TransProbCalcFactory.getTaylorCalculator(Q, sense, SCALING);
        TransProbCalculator pade = TransProbCalcFactory.getPadeCalculator(Q, sense, SCALING);

        for (double t : BRANCH_LENGTHS) {
            double[][] E = pt(exact, t);
            assertTrue("Padé vs exact at t=" + t + ": " + maxDiff(pt(pade, t), E), maxDiff(pt(pade, t), E) <= PADE_VS_EXACT_TOLERANCE);
            assertTrue("Taylor vs exact at t=" + t + ": " + maxDiff(pt(taylor, t), E), maxDiff(pt(taylor, t), E) <= TAYLOR_TOLERANCE);
        }
    }

    @Test
    public void fdsModelPadeMatchesTaylor() {
        int compared = 0;

        for (double z : FDS_Z) {
            double[] Q = new SwMutSel(mutation(), new Fitness(FITNESS.clone()), new FitnessFDS(z, 1)).getQ();
            double norm = norm1(Q) * SCALING;

            TransProbCalculator taylor = TransProbCalcFactory.getTaylorCalculator(Q, sense, SCALING);
            TransProbCalculator pade = TransProbCalcFactory.getPadeCalculator(Q, sense, SCALING);

            for (double t : BRANCH_LENGTHS) {
                if (norm * t > TAYLOR_MAX_NORM) continue;
                double diff = maxDiff(pt(pade, t), pt(taylor, t));
                assertTrue("Padé vs Taylor at Z=" + z + ", t=" + t + ": " + diff, diff <= TAYLOR_TOLERANCE);
                compared++;
            }
        }

        // 8 + 8 + 6 branch lengths are within Taylor's range for the three values of Z
        assertEquals("number of comparisons", 22, compared);
    }

    @Test
    public void fdsModelPadeGivesProbabilities() {
        for (double z : FDS_Z) {
            double[] Q = new SwMutSel(mutation(), new Fitness(FITNESS.clone()), new FitnessFDS(z, 1)).getQ();
            TransProbCalculator pade = TransProbCalcFactory.getPadeCalculator(Q, sense, SCALING);

            for (double t : BRANCH_LENGTHS) {
                double[][] P = pt(pade, t);
                for (int i : sense) {
                    double sum = 0;
                    for (int j : sense) {
                        assertTrue("P[" + i + "][" + j + "] out of range at Z=" + z + ", t=" + t,
                                P[i][j] >= -PROBABILITY_TOLERANCE && P[i][j] <= 1 + PROBABILITY_TOLERANCE);
                        sum += P[i][j];
                    }
                    assertEquals("row " + i + " sum at Z=" + z + ", t=" + t, 1.0, sum, PROBABILITY_TOLERANCE);
                }
            }
        }
    }

    private static SwMut mutation() {
        return new SwMut(KAPPA, SCALING, PI);
    }

    private static double[][] pt(TransProbCalculator calculator, double t) {
        double[][] P = new double[GeneticCode.CODON_STATES][GeneticCode.CODON_STATES];
        calculator.getTransitionProbabilities(P, t);
        return P;
    }

    private double maxDiff(double[][] a, double[][] b) {
        double max = 0;
        for (int i : sense) for (int j : sense) max = Math.max(max, Math.abs(a[i][j] - b[i][j]));
        return max;
    }

    /** Maximum absolute column sum of the sense-codon rate matrix (row-major). */
    private double norm1(double[] Q) {
        int n = sense.length;
        double max = 0;
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int i = 0; i < n; i++) sum += Math.abs(Q[i * n + j]);
            max = Math.max(max, sum);
        }
        return max;
    }
}
