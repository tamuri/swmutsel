package swmutsel;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import swmutsel.cli.Arguments;
import swmutsel.cli.ArgumentsProcessor;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Tests for site-wise fitness optimisation (-sites 1-3, with the mutation parameters and branch
 * lengths fixed) on the simulated dataset described in {@link FixAllTest}, with and without a penalty.
 *
 * Starting fitnesses are the defaults and single optimiser run per site, so search is deterministic. 
 * With only 5 taxa the fitnesses are poorly determined, so we only compare the maximised per-site
 * log-likelihoods.
 */
public class FitnessOptimisationTest {
    private static final String TREE = "((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05,E:0.2);";

    // Per-site values are printed to 7 d.p. The optimiser itself is only converged to about 1e-5.
    private static final double LNL_TOLERANCE = 1e-4;

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    @Test
    public void unpenalisedFitnessOptimisation() throws Exception {
        Map<Integer, Double> lnL = run("unpenalised");

        assertEquals("site 1", -1.6180196, lnL.get(1), LNL_TOLERANCE);
        assertEquals("site 2", -6.5081828, lnL.get(2), LNL_TOLERANCE);
        assertEquals("site 3", -1.6181335, lnL.get(3), LNL_TOLERANCE);
    }

    @Test
    public void mvnPenalisedFitnessOptimisation() throws Exception {
        Map<Integer, Double> lnL = run("penalised", "-p", "mvn,1");

        // Includes the penalty, so these are lower than the unpenalised values
        assertEquals("site 1", -4.1809559, lnL.get(1), LNL_TOLERANCE);
        assertEquals("site 2", -8.5945892, lnL.get(2), LNL_TOLERANCE);
        assertEquals("site 3", -4.1809417, lnL.get(3), LNL_TOLERANCE);
    }

    /**
     * Runs the analysis for sites 1-3 and returns the per-site log-likelihoods for those sites from the MLE file.
     */
    private Map<Integer, Double> run(String label, String... extraArgs) throws Exception {
        File resources = new File(System.getProperty("test.resources", "src/test/resources"));
        String name = new File(tmp.getRoot(), label).getPath();

        String[] fixed = new String[]{
                "-n", name,
                "-s", new File(resources, "sim300.phy").getPath(),
                "-t", TREE,
                "-gc", "standard",
                "-k", "2",
                "-pi", "0.25,0.25,0.25,0.25",
                "-c", "1",
                "-sites", "1-3",
                "-T", "2"
        };

        String[] args = new String[fixed.length + extraArgs.length];
        System.arraycopy(fixed, 0, args, 0, fixed.length);
        System.arraycopy(extraArgs, 0, args, fixed.length, extraArgs.length);

        Arguments parsed = ArgumentsProcessor.parse(args);
        parsed.getHandler().invoke();

        File mle = new File(name + "_MLE.txt");
        assertTrue("output file " + mle + " was not written", mle.exists());

        List<String> lines = Files.readAllLines(mle.toPath(), StandardCharsets.UTF_8);

        // The per-site table follows the argument block and starts with this header
        int header = lines.indexOf("site\tlnL\trel_w\trel_g");
        assertTrue("no per-site table in " + mle, header >= 0);

        // The table also lists sites that share a pattern with sites 1-3, so only keep sites 1-3
        Map<Integer, Double> lnL = new HashMap<Integer, Double>();
        for (String line : lines.subList(header + 1, lines.size())) {
            if (line.trim().isEmpty()) continue;
            String[] f = line.split("\t");
            assertEquals("bad line: " + line, 4, f.length);
            int site = Integer.parseInt(f[0]);
            if (site <= 3) lnL.put(site, Double.parseDouble(f[1]));
        }

        assertEquals("expected results for sites 1-3", 3, lnL.size());
        return lnL;
    }
}
