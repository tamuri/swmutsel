package swmutsel;

import org.junit.internal.TextListener;
import org.junit.runner.Description;
import org.junit.runner.JUnitCore;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;

import java.io.PrintStream;

/**
 * Runs the given JUnit test classes like org.junit.runner.JUnitCore, but prints the name of each test class and
 * test as it starts (instead of a dot), so you can see which test is running. The summary at the end (failures,
 * "OK (n tests)" or "Tests run: n, Failures: m") is JUnit's usual one. Used by "ant test".
 *
 * Usage: java swmutsel.VerboseTestRunner TestClass [TestClass ...]
 */
public class VerboseTestRunner {

    public static void main(String[] args) throws Exception {
        Class<?>[] classes = new Class<?>[args.length];
        for (int i = 0; i < args.length; i++) classes[i] = Class.forName(args[i]);

        JUnitCore core = new JUnitCore();
        core.addListener(new NamesListener(System.out));
        Result result = core.run(classes);

        System.exit(result.wasSuccessful() ? 0 : 1);
    }

    private static final class NamesListener extends TextListener {
        private final PrintStream out;
        private String currentClass = null;

        NamesListener(PrintStream out) {
            super(out);
            this.out = out;
        }

        @Override
        public void testStarted(Description description) {
            if (!description.getClassName().equals(currentClass)) {
                currentClass = description.getClassName();
                out.println();
                out.println("=== " + currentClass);
            }
            out.println("--- " + description.getMethodName());
            out.flush();
        }

        @Override
        public void testFailure(Failure failure) {
            out.println("--- FAILED: " + failure.getDescription().getMethodName());
            out.flush();
        }

        @Override
        public void testIgnored(Description description) {
            out.println("--- IGNORED: " + description.getMethodName());
            out.flush();
        }
    }
}
