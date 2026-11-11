package swmutsel;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Regression test for the persistent positive selection model (-fds; Tamuri and dos Reis, 2022) on one site of the
 * simulated dataset described in {@link FixAllTest}, with the mutation parameters fixed at their generating values.
 *
 * FDSHandler writes its results to "swMutSel_FDS.&lt;name&gt;" in the current working directory, so we run it in a
 * separate JVM whose working directory is a temporary folder. The expected values were recorded from the current
 * implementation (Java 17) and match the original fds branch exactly. Takes about 3 seconds.
 */
public class FDSTest {
    private static final String TREE = "((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05,E:0.2);";
    private static final String NAME = "fdstest";

    // The results file is a single line of 32 space-separated fields:
    //  1    site
    //  2-3  neutral non-synonymous and synonymous rates
    //  4-7  no-selection model (all fitnesses 0): lnL, non-syn, syn, expected subs
    //  8    lnL of the optimised FDS model
    //  9-28 optimised fitnesses (canonical amino acid order)
    //  29   optimised Z
    //  30-32 FDS model: non-syn, syn, expected subs
    private static final int FIELDS = 32;
    private static final double EXPECTED_NO_SELECTION_LNL = -4.805478;
    private static final double EXPECTED_FDS_LNL = -1.618020;
    private static final double EXPECTED_Z = 0.055673;

    // Values are printed to 6 d.p. Z is weakly determined by one site, so allow it a little more slack.
    private static final double LNL_TOLERANCE = 1e-5;
    private static final double Z_TOLERANCE = 1e-4;

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    @Test
    public void fdsSiteOptimisation() throws Exception {
        File resources = new File(System.getProperty("test.resources", "src/test/resources")).getAbsoluteFile();

        ProcessBuilder pb = new ProcessBuilder(
                new File(System.getProperty("java.home"), "bin/java").getPath(),
                "-cp", System.getProperty("java.class.path"),
                "swmutsel.Run",
                "-fds",
                "-n", NAME,
                "-s", new File(resources, "sim300.phy").getPath(),
                "-t", TREE,
                "-gc", "standard",
                "-k", "2",
                "-pi", "0.25,0.25,0.25,0.25",
                "-c", "1",
                "-site", "3");
        pb.directory(tmp.getRoot());
        pb.redirectErrorStream(true);
        Process process = pb.start();

        ByteArrayOutputStream out = new ByteArrayOutputStream();
        InputStream in = process.getInputStream();
        byte[] buffer = new byte[4096];
        int n;
        while ((n = in.read(buffer)) != -1) out.write(buffer, 0, n);
        String log = out.toString("UTF-8");

        assertEquals("exit code, output:\n" + log, 0, process.waitFor());

        File results = new File(tmp.getRoot(), "swMutSel_FDS." + NAME);
        assertTrue("results file was not written, output:\n" + log, results.exists());

        String[] f = new String(Files.readAllBytes(results.toPath()), StandardCharsets.UTF_8).trim().split("\\s+");
        assertEquals("number of fields", FIELDS, f.length);
        assertEquals("site", "3", f[0]);

        double noSelectionLnL = Double.parseDouble(f[3]);
        double fdsLnL = Double.parseDouble(f[7]);
        double z = Double.parseDouble(f[28]);

        assertEquals("no-selection log-likelihood", EXPECTED_NO_SELECTION_LNL, noSelectionLnL, LNL_TOLERANCE);
        assertEquals("FDS log-likelihood", EXPECTED_FDS_LNL, fdsLnL, LNL_TOLERANCE);
        assertEquals("Z", EXPECTED_Z, z, Z_TOLERANCE);
        assertTrue("optimised model should fit at least as well as the no-selection model", fdsLnL >= noSelectionLnL);
    }
}
