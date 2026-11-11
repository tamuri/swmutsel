package swmutsel;

import org.junit.Test;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.InputStream;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * "-help" prints the description, synopsis and option help, then exits. Because it calls System.exit() we run it
 * in a separate JVM, using the same classpath as the tests.
 */
public class HelpTest {

    @Test
    public void helpPrintsUsageAndExitsCleanly() throws Exception {
        ProcessBuilder pb = new ProcessBuilder(
                new File(System.getProperty("java.home"), "bin/java").getPath(),
                "-cp", System.getProperty("java.class.path"),
                "swmutsel.Run", "-help");
        pb.redirectErrorStream(true);
        Process process = pb.start();

        ByteArrayOutputStream out = new ByteArrayOutputStream();
        InputStream in = process.getInputStream();
        byte[] buffer = new byte[4096];
        int n;
        while ((n = in.read(buffer)) != -1) out.write(buffer, 0, n);

        int exitCode = process.waitFor();
        String text = out.toString("UTF-8");

        assertFalse("help threw an exception:\n" + text, text.contains("Exception"));
        assertEquals("exit code, output:\n" + text, 0, exitCode);
        assertTrue("no description in output:\n" + text, text.contains("This is swMutSel"));
        assertTrue("no synopsis in output:\n" + text, text.contains("SYNOPSIS"));
        assertTrue("no option help in output:\n" + text, text.contains("-sequences"));
    }
}
