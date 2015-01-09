package swmutsel.results;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Charsets;
import com.google.common.io.Files;
import swmutsel.Constants;
import swmutsel.cli.jc.GeneticCodeConverter;
import swmutsel.utils.GeneticCode;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;

/**
 * Given the model files produced by ModelWriter (S.txt etc.), this parses those files and produced the distribution
 * of selection coefficients for mutations and substitutions (all & non-synonymous only for both).
 * <p/>
 * TODO: PB2 plots significant vs. other sites (i.e. the black + red histogram in paper)
 * TODO: Really needs a less hard-coded way to correct construct the histogram bins, allowing some args etc.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 * @see ModelWriter
 * @see FitnessExtractor
 */
public class DistributionWriter {

    @Parameter(names = "-gc", converter = GeneticCodeConverter.class, required = true)
    private GeneticCode geneticCode;

    private static final String ROOT_DIR = "./";

    private static final double HI = 21D;
    private static final double LOW = -21D;
    private static final int MUTS_BINS = 168; // Means each bin is (21 - -21) / 168 = 0.25
    private static final int SUBS_BINS = 167; // So that the bins are symmetric around zero for substitutions (reversibility)
    private static final int S_LIMIT = 10;

    private static final String MUTATIONS_FILENAME = "distribution.mutations.csv";
    private static final String SUBSTITUTIONS_FILENAME = "distribution.substitutions.csv";

    public static void main(String[] args) throws Exception {
        DistributionWriter dw = new DistributionWriter();
        JCommander j = new JCommander(dw);
        j.parse(args);
        dw.run();
    }

    public void run() throws Exception {
        // We need to know the genetic code to determine non-synonymous changes
        //char[] aaCode = GeneticCode.CURRENT_CODE.toCharArray();

        double[] mutsAll = new double[MUTS_BINS];
        double[] mutsNonSyn = new double[MUTS_BINS];
        double[] subsAll = new double[SUBS_BINS];
        double[] subsNonSyn = new double[SUBS_BINS];

        // Get the neutral mutation matrix, Q0
        BufferedReader Q0Reader = new BufferedReader(new FileReader(ROOT_DIR + Constants.Q0_FILENAME));
        String[] Q0Parts = Q0Reader.readLine().split("\\s");
        double[] Q0 = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];
        for (int i = 0; i < (GeneticCode.CODON_STATES * GeneticCode.CODON_STATES); i++) {
            Q0[i] = Double.parseDouble(Q0Parts[i]);
        }
        Q0Reader.close();

        BufferedReader SReader = new BufferedReader(new FileReader(ROOT_DIR + Constants.S_FILENAME));
        BufferedReader QSReader = new BufferedReader(new FileReader(ROOT_DIR + Constants.QS_FILENAME));
        BufferedReader PiReader = new BufferedReader(new FileReader(ROOT_DIR + Constants.PI_FILENAME));

        double mutsAllDenom = 0, mutsNonSynDenom = 0, subsAllDenom = 0, subsNonSynDenom = 0;

        String SLine;
        while ((SLine = SReader.readLine()) != null) {
            String QSLine = QSReader.readLine();
            String PiLine = PiReader.readLine();

            String[] SParts = SLine.split("\\s"); // 64 * 64 = 4096 entries
            String[] QSParts = QSLine.split("\\s"); // 64 * 64 = 4096 entries
            String[] PiParts = PiLine.split("\\s"); // 64 entries

            double mutsAllSiteSum = 0, mutsNonSynSiteSum = 0, subsAllSiteSum = 0, subsNonSynSiteSum = 0;

            for (int i = 0; i < (GeneticCode.CODON_STATES * GeneticCode.CODON_STATES); i++) {

                int codon = (i - i % GeneticCode.CODON_STATES) / GeneticCode.CODON_STATES;
                boolean isCodonChange = codon != (i % GeneticCode.CODON_STATES);
                // boolean isNonSynChange = aaCode[codon] != aaCode[i % GeneticCode.CODON_STATES];
                boolean isNonSynChange = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(codon) !=
                        GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i % GeneticCode.CODON_STATES);

                double piValue = Double.parseDouble(PiParts[codon]);
                double qValue = Double.parseDouble(QSParts[i]);
                double deltaS = Double.parseDouble(SParts[i]);

                deltaS = Math.max(deltaS, -S_LIMIT);
                deltaS = Math.min(deltaS, S_LIMIT);

                double val = deltaS - LOW;
                int mutsBin = (int) (MUTS_BINS * (val / (HI - LOW))); // this just drops any decimal places i.e. (int) Math.floor(double)
                int subsBin = (int) (SUBS_BINS * (val / (HI - LOW))); // perhaps it should be (int) Math.round(SUBS_BIN * ... - LOW) etc.

                if (isNonSynChange) {
                    mutsNonSynSiteSum += piValue * Q0[i];
                    subsNonSynSiteSum += piValue * qValue;

                    mutsNonSyn[mutsBin] += piValue * Q0[i];
                    subsNonSyn[subsBin] += piValue * qValue;
                }

                if (isCodonChange) {
                    mutsAllSiteSum += piValue * Q0[i];
                    subsAllSiteSum += piValue * qValue;

                    mutsAll[mutsBin] += piValue * Q0[i];
                    subsAll[subsBin] += piValue * qValue;
                }
            }

            mutsAllDenom += mutsAllSiteSum;
            mutsNonSynDenom += mutsNonSynSiteSum;
            subsAllDenom += subsAllSiteSum;
            subsNonSynDenom += subsNonSynSiteSum;
        }

        SReader.close();
        QSReader.close();
        PiReader.close();

        int mutsStart = (int) (MUTS_BINS * ((-S_LIMIT - LOW) / (HI - LOW)));
        int mutsEnd = (int) (MUTS_BINS * ((S_LIMIT - LOW) / (HI - LOW))) + 1;
        double bin = -S_LIMIT;
        BufferedWriter mutsWriter = Files.newWriter(new File(MUTATIONS_FILENAME), Charsets.US_ASCII);
        for (int i = mutsStart; i < mutsEnd; i++) {
            mutsAll[i] /= mutsAllDenom;
            mutsNonSyn[i] /= mutsNonSynDenom;
            mutsWriter.write(String.format("%s,%s,%s\n", bin, mutsAll[i], mutsNonSyn[i]));
            bin += (HI - LOW) / MUTS_BINS;
        }
        mutsWriter.close();

        int subsStart = (int) (SUBS_BINS * ((-S_LIMIT - LOW) / (HI - LOW)));
        int subsEnd = (int) (SUBS_BINS * ((S_LIMIT - LOW) / (HI - LOW))) + 1;
        bin = -S_LIMIT;
        BufferedWriter subsWriter = Files.newWriter(new File(SUBSTITUTIONS_FILENAME), Charsets.US_ASCII);
        for (int i = subsStart; i < subsEnd; i++) {
            subsAll[i] /= subsAllDenom;
            subsNonSyn[i] /= subsNonSynDenom;
            subsWriter.write(String.format("%s,%s,%s\n", bin, subsAll[i], subsNonSyn[i]));
            bin += (HI - LOW) / MUTS_BINS; // NOTE: it's still 0.25!
        }
        subsWriter.close();
    }

}
