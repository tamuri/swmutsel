package swmutsel.results;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import com.google.common.base.Charsets;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.io.Files;
import com.google.common.io.LineProcessor;
import com.google.common.primitives.Doubles;
import swmutsel.Constants;
import swmutsel.model.SwMut;
import swmutsel.model.SwMutSel;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessFDS;
import swmutsel.utils.Functions;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.PhyloUtils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

/**
 * Takes a list of fitnesses (usually made by running FitnessExtractor) and swMutSel0 global parameters and
 * writes the following files:
 * <p/>
 * 1. Q0.txt - the neutral substitution rate matrix for each site
 * 2. QS.txt - the substitution rate matrix, with selection
 * 3. S.txt - the selection coefficient (S_ij) matrix
 * 4. PiS.txt - the codon frequencies
 * 5. PiAA.txt - the amino acid frequencies
 * <p/>
 * DistributionWriter then uses these to create a file containing the distribution of selection coefficients.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 * @see FitnessExtractor
 * @see DistributionWriter
 */
public class ModelWriter {
    All.Options o;
    SwMut mutation;

    FileWriter outS, outPiS, outQS, outPiAA;
    List<Integer> aminoAcids = ImmutableList.copyOf(Lists.<Integer>newArrayList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19));


    public ModelWriter(All.Options o) throws Exception {
        this.o = o;
    }

    public void run() throws Exception {

        mutation = new SwMut(o.kappa, o.scaling, o.pi);
        mutation.build();

        // Neutral Q
        DoubleMatrix2D q0 = DoubleFactory2D.dense.make(GeneticCode.CODON_STATES, GeneticCode.CODON_STATES);

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                if (i == j) q0.setQuick(i, j, 0);
                else q0.setQuick(i, j, mutation.getMutationRate(i, j));
            }
        }

        // Q0 rows sum to 0
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            double sum = 0;
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                sum += q0.getQuick(i, j);
            }
            q0.setQuick(i, i, -sum);
        }

        FileWriter outQ0 = new FileWriter(new File(Constants.Q0_FILENAME));
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            outQ0.write(Doubles.join(" ", q0.viewRow(i).toArray()));
            if (i < GeneticCode.CODON_STATES - 1) outQ0.write(" ");
        }

        outQ0.close();

        // S, Q with selection and codon pi files
        outS = new FileWriter(new File(Constants.S_FILENAME));
        outPiS = new FileWriter(new File(Constants.PI_FILENAME));
        outQS = new FileWriter(new File(Constants.QS_FILENAME));
        outPiAA = new FileWriter(new File(Constants.PIAA_FILENAME));

        Files.readLines(new File(Constants.F_FILENAME), Charsets.UTF_8, new FitnessProcessor());

        outS.close();
        outPiS.close();
        outQS.close();
        outPiAA.close();
    }

    class FitnessProcessor implements LineProcessor<Object> {
        @Override
        public boolean processLine(String line) throws IOException {

            List<Double> fitnesses = Lists.transform(Arrays.asList(line.split(" ")), Functions.stringToDouble());

            SwMutSel model;
            double[] S;

            if (fitnesses.size() == 20) {
                Fitness f = new Fitness(Doubles.toArray(fitnesses));
                model = new SwMutSel(mutation, f);
                S = getS(f);
            } else if (fitnesses.size() == 21) {
                // 20 fitnesses and one Z parameter

                double[] fitness = new double[20];
                for (int i = 0; i < 20; i++) {
                    fitness[i] = fitnesses.get(i);
                }

                Fitness f = new Fitness(fitness);
                FitnessFDS fds = new FitnessFDS(fitnesses.get(20), 1);
                model = new SwMutSel(mutation, f, fds);
                S = getS(f, fds);
            } else {
                throw new RuntimeException("Line has " + fitnesses.size() + " numbers?");
            }

            model.build();

            double[] QS = model.getFullQ();
            double[] PiS = model.getCodonFrequencies();

            outS.write(String.format("%s\n", Doubles.join(" ", S)));
            outPiS.write(String.format("%s\n", Doubles.join(" ", PiS)));
            outQS.write(String.format("%s\n", Doubles.join(" ", QS)));
            outPiAA.write(String.format("%s\n", Doubles.join(" ", PhyloUtils.getAminoAcidFrequencies(model.getCodonFrequencies()))));

            return true;
        }

        @Override
        public Object getResult() {
            return null;
        }

        private double[] getS(Fitness f) {
            double[] fullS = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    if (i == j) {
                        fullS[i * GeneticCode.CODON_STATES + j] = 0;
                    } else {
                        int aa_from = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                        int aa_to = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);
                        if (aa_to < 0) {
                            // going to STOP codon
                            fullS[i * GeneticCode.CODON_STATES + j] = -Constants.FITNESS_BOUND;
                        } else if (aa_from < 0) {
                            // going from a STOP codon to any non-STOP codon
                            fullS[i * GeneticCode.CODON_STATES + j] = Constants.FITNESS_BOUND;
                        } else {
                            // both amino acids that occur at this site
                            fullS[i * GeneticCode.CODON_STATES + j] = f.get()[aa_to] - f.get()[aa_from];
                        }
                    }
                }
            }
            return fullS;
        }

        private double[] getS(Fitness fitness, FitnessFDS fds) {
            double[] fullS = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

            for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
                for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                    if (i == j) {
                        fullS[i * GeneticCode.CODON_STATES + j] = 0;
                    } else {
                        int aa_from = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                        int aa_to = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);
                        if (aa_to < 0) {
                            // going to STOP codon
                            fullS[i * GeneticCode.CODON_STATES + j] = -Constants.FITNESS_BOUND;
                        } else if (aa_from < 0) {
                            // going from a STOP codon to any non-STOP codon
                            fullS[i * GeneticCode.CODON_STATES + j] = Constants.FITNESS_BOUND;
                        } else {
                            // both amino acids that occur at this site
                            fullS[i * GeneticCode.CODON_STATES + j] = fitness.get()[aa_to] - fitness.get()[aa_from] + fds.get();
                        }
                    }
                }
            }
            return fullS;
        }

    }



}
