package swmutsel.model;

import cern.colt.matrix.DoubleMatrix2D;
import com.google.common.primitives.Ints;
import swmutsel.Constants;
import swmutsel.model.parameters.Fitness;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.Pair;

/**
 * Using ZY's notation from CME pp.68
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 15:24
 */
final public class SwMutSel extends SubstitutionModel {
    private static final long serialVersionUID = 3332802987971499184L;
    private final SwMut mutation;
    private final Fitness fitness;

    private final int[] senseCodons;
    private final int numSenseCodons;

    private final double[] senseCodonFrequencies;
    private final double[] lambda;
    private final double[][] U;
    private final double[][] Uinv;

    private volatile transient boolean changed = true;
    private volatile transient TransProbCalculator calculator;

    private double expectedSubstitutions;
    private Pair<Double, Double> proportionsSynAndNonSyn;

    public SwMutSel(SwMut mutation, Fitness fitness) {
        this.mutation = mutation;
        this.fitness = fitness;
        super.clearParameters();
        super.addParameters(fitness); // We're assuming mutation component is fixed here

        this.senseCodons = GeneticCode.getInstance().getSenseCodons();
        this.numSenseCodons = this.senseCodons.length;

        this.senseCodonFrequencies = new double[this.numSenseCodons];
        this.U = new double[this.numSenseCodons][this.numSenseCodons];
        this.Uinv = new double[this.numSenseCodons][this.numSenseCodons];
        this.lambda = new double[this.numSenseCodons];

        build();
    }


    @Override
    final public void build() {
        setSenseCodonFrequencies();
        setEigenFactors(makeB(makeQ(), numSenseCodons, senseCodonFrequencies));
        changed = true;
    }

    @Override
    final public boolean parametersValid() {
        return true;
    }

    private void setSenseCodonFrequencies() {
        double z = 0;
        for (int i = 0; i < numSenseCodons; i++) {
            senseCodonFrequencies[i] = this.mutation.getCodonFrequency(senseCodons[i]) *
                    Math.exp(fitness.get()[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(senseCodons[i])]);
            z += senseCodonFrequencies[i];
        }

        for (int i = 0; i < numSenseCodons; i++) {
            senseCodonFrequencies[i] /= z;
        }
    }

    private void setEigenFactors(DoubleMatrix2D B) {
        setEigenFactors(B, lambda, numSenseCodons, senseCodonFrequencies, U, Uinv);
        for (int i = 0; i < numSenseCodons; i++) {
            lambda[i] *= mutation.getBranchScaling().get();
        }
    }

    private double[] makeQ() {
        double[] F = this.fitness.get();
        double[] Q = new double[numSenseCodons * numSenseCodons];

        double expectedSubstitutions = 0;
        double totalSynonymous = 0;
        double totalNonSynonymous = 0;

        for (int i = 0; i < numSenseCodons; i++) {
            for (int j = 0; j < numSenseCodons; j++) {

                if (i == j) continue;

                // If mutation rate is 0, substitution rate is 0 (e.g. multiple nucleotide changes)
                double mutationRate = this.mutation.getMutationRate(senseCodons[i], senseCodons[j]);
                //if (mutationRate == 0) continue;

                int Ai = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(senseCodons[i]);
                int Aj = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(senseCodons[j]);

                double Fi = F[Ai];
                double Fj = F[Aj];

                double hS = getRelativeFixationProbability(Fj - Fi);

                Q[i * numSenseCodons + j] = mutationRate * hS;

                double rate = senseCodonFrequencies[i] * Q[i * numSenseCodons + j];
                if (Ai == Aj) {
                    totalSynonymous += rate;
                } else {
                    totalNonSynonymous += rate;
                }
            }
        }

        this.proportionsSynAndNonSyn = Pair.of(totalNonSynonymous, totalSynonymous);

        // Each row sums to 0
        for (int row = 0; row < numSenseCodons; row++) {
            double rowSum = 0;
            for (int column = 0; column < numSenseCodons; column++) {
                rowSum += Q[row * numSenseCodons + column];
            }
            Q[row * numSenseCodons + row] = -rowSum;
            expectedSubstitutions += rowSum * senseCodonFrequencies[row];
        }

        this.expectedSubstitutions = expectedSubstitutions;

        return Q;
    }

    public final double getExpectedSubsPerSite() {
        return expectedSubstitutions;
    }

    public final Pair<Double, Double> getRhoNonSynAndSyn() {
        return proportionsSynAndNonSyn;
    }

    public final double[] getFullQ() {
        double[] Q = makeQ();

        double[] fullQ = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

        // the S to go from a STOP codon to any sense codon = "very high"
        double stopS = getRelativeFixationProbability(Constants.FITNESS_BOUND - -Constants.FITNESS_BOUND);

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                int aa_from = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                int aa_to = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);
                // if going to a stop codon
                if (aa_to < 0) {
                    // rate is 0
                    fullQ[i * GeneticCode.CODON_STATES + j] = 0;
                    // if coming from a stop codon and going to a non-stop codon
                } else if (aa_from < 0) {
                    fullQ[i * GeneticCode.CODON_STATES + j] = mutation.getMutationRate(i, j) * stopS;
                } else {
                    fullQ[i * GeneticCode.CODON_STATES + j] = Q[Ints.indexOf(this.senseCodons, i) * this.numSenseCodons + Ints.indexOf(this.senseCodons, j)];
                }
            }
        }

        for (int row = 0; row < GeneticCode.CODON_STATES; row++) {
            double sum = 0;
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {
                sum += fullQ[row * GeneticCode.CODON_STATES + j];
            }
            fullQ[row * GeneticCode.CODON_STATES + row] = -sum;
        }

        return fullQ;
    }

    @Override
    public final double[] getCodonFrequencies() {
        // return frequencies for all 64 codons - compare to setSenseCodonFrequencies()
        double[] frequencies = new double[GeneticCode.CODON_STATES];
        for (int i = 0; i < numSenseCodons; i++) {
            frequencies[senseCodons[i]] = senseCodonFrequencies[i];
        }
        return frequencies;
    }

    @Override
    public TransProbCalculator getPtCalculator() {
        if (changed || calculator == null) {
            synchronized (this) {
                if (changed || calculator == null) {
                    this.calculator = TransProbCalcFactory.getPtCalculator(lambda, U, Uinv, senseCodons, false);
                    changed = false;
                }
            }
        }
        return this.calculator;
    }

}
