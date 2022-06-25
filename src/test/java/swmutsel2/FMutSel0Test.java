package swmutsel2;


import com.google.common.collect.Lists;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;
import swmutsel.model.FMutSel0;
import swmutsel.model.SubstitutionModel;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.utils.GeneticCode;
import swmutsel2.alignment.AlignmentUtils;
import swmutsel2.alignment.SitePatterns;
import swmutsel2.likelihood.ClassicCalculator;
import swmutsel2.likelihood.DataController;
import swmutsel2.traversal.TipPartials;
import swmutsel2.traversal.TraversalInfo;
import swmutsel2.tree.Tree;
import swmutsel2.tree.TreeReader;

import java.io.File;
import java.net.URL;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import static org.junit.Assert.assertEquals;

/**
 * FMutSel0 Tester.
 *
 * Some end-to-end testing.
 *
 * @author Asif Tamuri
 * @version 1.0
 */
@RunWith(Parameterized.class)
public class FMutSel0Test {
    private static final double DELTA = 1.0e-6;

    @Parameters
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {"mit4s.tree", "mit4s.phylip", GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE, -42461.784344},
                {"data_401x759.tree", "data_401x759.phylip", GeneticCode.STANDARD_CODE, -102277.157401}
        });
    }

    @Parameter(value = 0)
    public String treePath;

    @Parameter(value = 1)
    public String alignmentPath;

    @Parameter(value = 2)
    public String geneticCode;

    @Parameter(value = 3)
    public double expected;

    private List<TraversalInfo> traversalData;

    private SitePatterns patterns;
    private SubstitutionModel model;
    
    @Before
    public void before() throws Exception {
        GeneticCode.setCode(geneticCode);

        Tree tree = TreeReader.read(getAbsolutePath(treePath));
        tree.addVirtualRoot();

        traversalData = tree.getTraversalData();
        patterns = AlignmentUtils.loadSequencePatterns(getAbsolutePath(alignmentPath), tree);

        TsTvRatio kappa = new TsTvRatio(2);
        model = new FMutSel0(kappa,
                new Omega(2),
                new BaseFrequencies(BaseFrequencies.getDefault()),
                new Fitness());
        model.build();
    }

    @After
    public void after() throws Exception {
    }

    @Test
    public void testSingleThreaded() throws Exception {
        final int threads = 1;

        final byte[][][] splitPattern = new byte[threads][][];
        patterns.split(splitPattern);
        final DataController controller = new DataController(splitPattern);
        
        final double[][] tipPartials = TipPartials.forCodons();

        ClassicCalculator calculator = new ClassicCalculator(controller.get(0).getConditionals(), splitPattern[0], tipPartials, model);
        final double[] likelihoods = calculator.getLogLikelihood(traversalData);

        final double[][] lnls = new double[threads][];
        lnls[0] = likelihoods;

        final double sum = getWeightedLogLikelihood(patterns, lnls);

        assertEquals(expected, sum, DELTA);
    }

    @Test
    public void testThreeThreads() throws Exception {
        final int threads = 3;

        final byte[][][] splitPattern = new byte[threads][][];
        patterns.split(splitPattern);
        final DataController controller = new DataController(splitPattern);

        final double[][] tipPartials = TipPartials.forCodons();

        ExecutorService service = Executors.newFixedThreadPool(threads);

        class MyCallable implements Callable<double[]> {
            private final List<TraversalInfo> traversalData;
            private final ClassicCalculator calculator;

            public MyCallable(double[][][] partials, byte[][] states, double[][] tips, SubstitutionModel model, List<TraversalInfo> traversalData) {
                this.traversalData = traversalData;
                this.calculator = new ClassicCalculator(partials, states, tips, model);
            }

            @Override
            public double[] call() throws Exception {
                return calculator.getLogLikelihood(traversalData);
            }
        }

        MyCallable[] mc = new MyCallable[threads];
        for (int i = 0; i < threads; i++) {
            final double[][][] p = controller.get(i).getConditionals();
            final byte[][] s = splitPattern[i];

            mc[i] = new MyCallable(p, s, tipPartials, model, traversalData);
        }

        List<Future<double[]>> results = Lists.newArrayList();
        for (int i = 0; i < threads; i++) {
            results.add(service.submit(mc[i]));
        }

        final double[][] likelihoods = new double[threads][];
        int pos = 0;
        for (Future<double[]> r : results) {
            likelihoods[pos++] = r.get();
        }

        final double sum = getWeightedLogLikelihood(patterns, likelihoods);

        assertEquals(expected, sum, DELTA);
    }

    private double getWeightedLogLikelihood(SitePatterns pattern, double[][] lnls) {
        int offset = 0;
        int[] weighting = pattern.getWeighting();

        double sum = 0;
        for (int i = 0; i < lnls.length; i++) {
            for (int j = 0; j < lnls[i].length; j++) {
                sum += weighting[offset] * lnls[i][j];
                offset++;
            }
        }
        return sum;
    }

    private String getAbsolutePath(String file) {
        URL url = this.getClass().getResource("/" + file);
        return new File(url.getFile()).getAbsolutePath();
    }
}
