package swmutsel;

import org.junit.After;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import swmutsel.cli.Arguments;
import swmutsel.cli.ArgumentsProcessor;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Test for the (single class) FMutSel0 analysis, which alternates branch length optimisation with
 * model parameter optimisation.
 */
public class FMutSel0OptimisationTest {
    private static final String TREE = "((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05,E:0.2);";
    private static final double EXPECTED_LNL = -2020.5400173;

    // BOBYQA stops on a trust-region radius, not on a change in lnL, so allow a little more slack than elsewhere
    private static final double LNL_TOLERANCE = 1e-2;

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    @After
    public void restoreConvergenceTolerance() {
        // FMutSelHandler tightens this static in each round, which would leak into other tests in the same JVM
        Constants.VALUE_CONVERGENCE_TOL = 5.0e-6;
    }

    @Test
    public void optimisedLogLikelihood() throws Exception {
        File resources = new File(System.getProperty("test.resources", "src/test/resources"));

        PrintStream original = System.out;
        ByteArrayOutputStream captured = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(captured, true, "UTF-8"));
            Arguments parsed = ArgumentsProcessor.parse(new String[]{
                    "-fmutsel0",
                    "-n", new File(tmp.getRoot(), "fmutsel0").getPath(),
                    "-s", new File(resources, "sim300.phy").getPath(),
                    "-t", TREE,
                    "-gc", "standard",
                    "-pi", "0.25,0.25,0.25,0.25",
                    "-T", "2"
            });
            parsed.getHandler().invoke();
        } finally {
            System.setOut(original);
        }

        String output = captured.toString("UTF-8");
        Matcher m = Pattern.compile("Maximum log-likelihood: (-?[0-9.]+)").matcher(output);
        assertTrue("no final log-likelihood in output:\n" + output, m.find());

        assertEquals("maximum log-likelihood", EXPECTED_LNL, Double.parseDouble(m.group(1)), LNL_TOLERANCE);
    }
}
