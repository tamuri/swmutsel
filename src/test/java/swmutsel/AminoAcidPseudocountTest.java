package swmutsel;

import com.google.common.collect.ArrayTable;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import swmutsel.cli.Arguments;
import swmutsel.cli.ArgumentsProcessor;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.PhyloUtils;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * When an amino acid is absent from the alignment, its observed frequency is 0. Starting fitnesses are the log of
 * the frequencies, so they would be -Infinity, giving NaN in the model and a crash in the eigen solver. -fmutsel0
 * therefore adds a pseudocount to the amino acid counts.
 */
public class AminoAcidPseudocountTest {
    private static final String TREE = "((A:0.1,B:0.1):0.05,(C:0.1,D:0.1):0.05,E:0.2);";

    // Log-likelihood at the starting FMutSel0 parameters (kappa = 1, omega = 1, observed pi, pseudocount fitnesses)
    private static final double EXPECTED_LNL = -357.9116843;

    @Rule
    public TemporaryFolder tmp = new TemporaryFolder();

    /** One sequence of 10 codons: 6 x TTT (Phe), 4 x GCT (Ala) -> only 2 of the 20 amino acids are observed. */
    private Table<String, Integer, Byte> twoAminoAcidTable() {
        GeneticCode.setCode(GeneticCode.STANDARD_CODE);
        byte phe = (byte) GeneticCode.getInstance().getCodonIndexFromNucleotides("TTT".toCharArray());
        byte ala = (byte) GeneticCode.getInstance().getCodonIndexFromNucleotides("GCT".toCharArray());

        Table<String, Integer, Byte> table = ArrayTable.create(Lists.newArrayList("seq"), Lists.newArrayList(1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
        for (int site = 1; site <= 10; site++) table.put("seq", site, site <= 6 ? phe : ala);
        return table;
    }

    @Test
    public void defaultPseudocountIsZero() {
        Table<String, Integer, Byte> table = twoAminoAcidTable();
        double[] expected = PhyloUtils.getObservedAminoAcidFrequencies(table, 0);
        double[] actual = PhyloUtils.getObservedAminoAcidFrequencies(table);

        int ala = GeneticCode.getInstance().getAminoAcidIndexByChar('A');
        int phe = GeneticCode.getInstance().getAminoAcidIndexByChar('F');

        for (int i = 0; i < actual.length; i++) assertEquals("amino acid " + i, expected[i], actual[i], 0);
        assertEquals(0.4, actual[ala], 1e-12);
        assertEquals(0.6, actual[phe], 1e-12);

        int unobserved = 0;
        for (double f : actual) if (f == 0) unobserved++;
        assertEquals(18, unobserved);
    }

    @Test
    public void pseudocountGivesEveryAminoAcidAPositiveFrequency() {
        Table<String, Integer, Byte> table = twoAminoAcidTable();
        double[] f = PhyloUtils.getObservedAminoAcidFrequencies(table, 0.1);

        int ala = GeneticCode.getInstance().getAminoAcidIndexByChar('A');
        int phe = GeneticCode.getInstance().getAminoAcidIndexByChar('F');

        double sum = 0;
        for (double x : f) {
            assertTrue("frequency should be positive", x > 0);
            sum += x;
        }
        assertEquals("frequencies sum to 1", 1.0, sum, 1e-12);

        // counts are 6 and 4 out of 10, plus 0.1 for each of the 20 amino acids: total 12
        assertEquals((4 + 0.1) / 12.0, f[ala], 1e-12);
        assertEquals((6 + 0.1) / 12.0, f[phe], 1e-12);
        assertEquals(0.1 / 12.0, f[GeneticCode.getInstance().getAminoAcidIndexByChar('W')], 1e-12);
    }

    /**
     * sim50.phy (5 taxa, 50 codon sites) does not contain every amino acid. Before the pseudocount, creating the
     * single class FMutSel0 model for it crashed with an ArrayIndexOutOfBoundsException.
     */
    @Test
    public void fmutsel0RunsOnAlignmentWithMissingAminoAcids() throws Exception {
        File resources = new File(System.getProperty("test.resources", "src/test/resources"));

        PrintStream original = System.out;
        ByteArrayOutputStream captured = new ByteArrayOutputStream();
        try {
            System.setOut(new PrintStream(captured, true, "UTF-8"));
            Arguments parsed = ArgumentsProcessor.parse(new String[]{
                    "-fmutsel0",
                    "-n", new File(tmp.getRoot(), "missing").getPath(),
                    "-s", new File(resources, "sim50.phy").getPath(),
                    "-t", TREE,
                    "-gc", "standard",
                    "-fix", "all"
            });
            parsed.getHandler().invoke();
        } finally {
            System.setOut(original);
        }

        String output = captured.toString("UTF-8");
        Matcher m = Pattern.compile("Total log-likelihood = (-?[0-9.]+)").matcher(output);
        assertTrue("no log-likelihood in output:\n" + output, m.find());

        // Recorded from the current implementation
        assertEquals("log-likelihood", EXPECTED_LNL, Double.parseDouble(m.group(1)), 1e-4);
    }
}
