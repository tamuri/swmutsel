package swmutsel;

import swmutsel.cli.Arguments;
import swmutsel.cli.ArgumentsProcessor;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Small e2e regression test for "-fix all" (log-likelihood only) on a simulated dataset.
 *
 * sim300.phy simulated with swMutSel: standard genetic code, kappa = 2, pi = 0.25 each, scaling = 1 
 * and single fitness vector in sim300.fitness.csv (repeated for every site). 
 *
 * Run with "ant test".
 */
public class FixAllTest {
    private static final String TREE = "((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05,E:0.2);";

    // Values recorded from the current implementation (Java 17)
    private static final double EXPECTED_TOTAL_LNL = -2037.9994596;
    private static final int EXPECTED_SITES = 300;
    private static final double EXPECTED_REL_W = 0.7561650; // same for every site
    private static final double EXPECTED_REL_G = 1.0663048;

    // The per-site values in the output file are rounded to 7 d.p., so their sum differs slightly from the total
    private static final double LNL_TOLERANCE = 1e-5;
    private static final double RATIO_TOLERANCE = 1e-6;

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    @Test
    public void logLikelihoodAtGeneratingParametersIsReproduced() throws Exception {
        File resources = new File(System.getProperty("test.resources", "src/test/resources"));
        String name = new File(tmp.getRoot(), "fixall").getPath();

        Arguments parsed = ArgumentsProcessor.parse(new String[]{
                "-n", name,
                "-s", new File(resources, "sim300.phy").getPath(),
                "-t", TREE,
                "-gc", "standard",
                "-k", "2",
                "-pi", "0.25,0.25,0.25,0.25",
                "-c", "1",
                "-f", new File(resources, "sim300.fitness.csv").getPath(),
                "-fix", "all"
        });
        parsed.getHandler().invoke();

        File out = new File(name + "_CALCULATE.txt");
        assertTrue("output file " + out + " was not written", out.exists());

        List<String> lines = Files.readAllLines(out.toPath(), StandardCharsets.UTF_8);
        assertEquals("site\tlnL\trel_w\trel_g", lines.get(0));

        int sites = 0;
        double total = 0;
        for (String line : lines.subList(1, lines.size())) {
            if (line.trim().isEmpty()) continue;
            String[] f = line.split("\t");
            assertEquals("bad line: " + line, 4, f.length);
            sites++;
            total += Double.parseDouble(f[1]);
            assertEquals("rel_w at site " + f[0], EXPECTED_REL_W, Double.parseDouble(f[2]), RATIO_TOLERANCE);
            assertEquals("rel_g at site " + f[0], EXPECTED_REL_G, Double.parseDouble(f[3]), RATIO_TOLERANCE);
        }

        assertEquals(EXPECTED_SITES, sites);
        assertEquals("total log-likelihood", EXPECTED_TOTAL_LNL, total, LNL_TOLERANCE);
    }
}
