package swmutsel.model;

import org.junit.Test;
import org.junit.Before;
import org.junit.After;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.utils.GeneticCode;

import java.io.*;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

/**
 * TransProbCalcFactory Tester.
 *
 * @author <Authors name>
 * @since <pre>Jun 18, 2015</pre>
 * @version 1.0
 */
public class TransProbCalcFactoryTest {

    private static final double DELTA = 1.0e-6;

    @Before
    public void before() throws Exception {
    }

    @After
    public void after() throws Exception {
    }

    /**
     *
     * Method: calculatePartialLikelihoodPair(final double[] parent, final double[] pCond,
     *      final double pLen, final double[] qCond, final double qLen, final double[] pTemp,
     *      final double[] qTemp)
     *
     */
    @Test
    public void testCalculatePartialLikelihoodPair() throws Exception {
        SubstitutionModel model = new FMutSel0(new TsTvRatio(1),
                new Omega(1),
                new BaseFrequencies(BaseFrequencies.getDefault()),
                new Fitness());
        model.build();

        double[] parent = new double[GeneticCode.CODON_STATES];
        model.getPtCalculator().calculatePartialLikelihoodPair(parent,
                scanVector("conditional/pCond.dat", 64),
                scanVector("conditional/pLen.dat", 1)[0],
                scanVector("conditional/qCond.dat", 64),
                scanVector("conditional/qLen.dat", 1)[0],
                new double[GeneticCode.CODON_STATES],
                new double[GeneticCode.CODON_STATES]);

        double[] result = scanVector("conditional/result.dat", 64);

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            assertEquals(result[i], parent[i], DELTA);
        }
    }

    private double[] scanVector(String file, int size) {
        double[] scan = new double[size];
        List<Double> numbers = readNumbers(file);

        for (int i = 0; i < numbers.size(); i++) {
            scan[i] = numbers.get(i);
        }
        return scan;
    }

    private double[][] scanMatrix(String file, int size_r, int size_c) {
        double[][] scan = new double[size_r][size_c];

        List<Double> numbers = readNumbers(file);

        int offset = 0;
        for (int i = 0; i < numbers.size(); i++) {
            for (int j = 0; j < numbers.size(); j++) {
                scan[i][j] = numbers.get(offset++);
            }
        }

        return scan;
    }

    private List<Double> readNumbers(String filePath) {
        List<Double> list = new ArrayList<>();
        File file = new File(getAbsolutePath(filePath));
        BufferedReader reader = null;

        try {
            reader = new BufferedReader(new FileReader(file));
            String text;

            while ((text = reader.readLine()) != null) {
                if (text.trim().length() == 0) continue;
                list.add(Double.parseDouble(text));
            }
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (reader != null) {
                    reader.close();
                }
            } catch (IOException ignored) {
            }
        }
        return list;
    }

    private String getAbsolutePath(String file) {
        URL url = this.getClass().getResource("/" + file);
        return new File(url.getFile()).getAbsolutePath();
    }
} 
