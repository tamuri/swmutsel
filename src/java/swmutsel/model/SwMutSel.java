package swmutsel.model;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.jet.math.Functions;
import com.google.common.primitives.Ints;
import swmutsel.Constants;
import swmutsel.MatrixArrayPool;
import swmutsel.model.parameters.Fitness;
import swmutsel.utils.GeneticCode;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 15:24
 */
public class SwMutSel extends SubstitutionModel {
    private static final long serialVersionUID = 3332802987971499184L;
    private final SwMut mutation;
    private final Fitness fitness;

    private final double[] senseCodonFrequencies;
    private DoubleMatrix1D lambda; // eigenvalue of B
    private final double[] U;
    private final double[] UInv;

    private final int[] senseCodons;

    public SwMutSel(SwMut mutation, Fitness fitness) {
        this.mutation = mutation;
        this.fitness = fitness;
        setParameters(fitness); // We're assuming mutation component is fixed here

        this.senseCodons = GeneticCode.getInstance().getSenseCodons();

        this.senseCodonFrequencies = new double[this.senseCodons.length];
        this.U = new double[this.senseCodons.length * this.senseCodons.length];
        this.UInv = new double[this.senseCodons.length * this.senseCodons.length];

        build();
    }

    @Override
    public void build() {
        setSenseCodonFrequencies();
        setEigenFactors(makeB(makeQ()));
    }

    @Override
    public boolean parametersValid() {
        return true;
    }

    private void setSenseCodonFrequencies() {
        double z = 0;
        for (int i = 0; i < senseCodons.length; i++) {
            senseCodonFrequencies[i] = this.mutation.getCodonFrequency(senseCodons[i]) *
                    Math.exp(fitness.get()[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(senseCodons[i])]);
            z += senseCodonFrequencies[i];
        }

        for (int i = 0; i < senseCodonFrequencies.length; i++) {
            senseCodonFrequencies[i] /= z;
        }
    }

    private void setEigenFactors(DoubleMatrix2D B) {
        EigenvalueDecomposition evdB = new EigenvalueDecomposition(B);

        lambda = DoubleFactory2D.dense.diagonal(evdB.getD());
        // we scale branch length by global parameter mu here so we only have do it once
        lambda.assign(Functions.mult(mutation.getBranchScaling().get()));

        DoubleMatrix2D R = evdB.getV();

        for (int i = 0; i < senseCodons.length; i++) {
            double piSqrt = Math.sqrt(senseCodonFrequencies[i]);
            double piInvSqrt = 1 / piSqrt;
            for (int j = 0; j < senseCodons.length; j++) {
                U[i * senseCodons.length + j] = piInvSqrt * R.getQuick(i, j);
                UInv[j * senseCodons.length + i] = piSqrt * R.getQuick(i, j); // inverse(R) == transpose(R)
            }
        }
    }

    private DoubleMatrix2D makeB(double[] Q) {
        DoubleMatrix2D B = DoubleFactory2D.dense.make(senseCodons.length, senseCodons.length);
        for (int i = 0; i < senseCodons.length; i++) {
            for (int j = 0; j < senseCodons.length; j++) {
                B.setQuick(i, j, Q[i * senseCodons.length + j] * Math.sqrt(senseCodonFrequencies[i]) / Math.sqrt(senseCodonFrequencies[j]));
            }
        }
        return B;
    }

    private double[] makeQ() {
        double[] F = this.fitness.get();
        double[] Q = new double[senseCodons.length * senseCodons.length];

        for (int i = 0; i < senseCodons.length; i++) {
            for (int j = 0; j < senseCodons.length; j++) {

                if (i == j) continue;

                // If mutation rate is 0, substitution rate is 0 (e.g. multiple nucleotide changes)
                double mutationRate = this.mutation.getMutationRate(senseCodons[i], senseCodons[j]);
                //if (mutationRate == 0) continue;

                int Ai = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(senseCodons[i]);
                int Aj = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(senseCodons[j]);

                double Fi = F[Ai];
                double Fj = F[Aj];

                double hS = getRelativeFixationProbability(Fj - Fi);

                Q[i * senseCodons.length + j] = mutationRate * hS;
            }
        }

        // Each row sums to 0
        for (int row = 0; row < senseCodons.length; row++) {
            double sum = 0;
            for (int column = 0; column < senseCodons.length; column++) {
                sum += Q[row * senseCodons.length + column];
            }
            Q[row * senseCodons.length + row] = -sum;
        }

        return Q;
    }

    public double getExpectedSubsPerSite() {
        double[] Q = makeQ();
        double sum = 0;
        for (int i = 0; i < senseCodons.length; i++) {
            sum += senseCodonFrequencies[i] * Q[i * senseCodons.length + i];
        }
        return -sum;
    }

    private double getRelativeFixationProbability(double S) {
        if (S == 0) return 1;
        if (S < -1e3) return 0;
        if (S > 1e3) return S;
        return S / -Math.expm1(-S);

//      return (S == 0) ? 1 : S / (1 - Math.exp(-S));
    }

    public double[] getFullQ() {
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
                    fullQ[i * GeneticCode.CODON_STATES + j] = Q[Ints.indexOf(this.senseCodons, i) * this.senseCodons.length + Ints.indexOf(this.senseCodons, j)];
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
    public double[] getCodonFrequencies() {
        // return frequencies for all 64 codons - compare to setSenseCodonFrequencies()
        double[] frequencies = new double[GeneticCode.CODON_STATES];
        for (int i = 0; i < senseCodons.length; i++) {
            frequencies[senseCodons[i]] = senseCodonFrequencies[i];
        }
        return frequencies;
    }


    @Override
    public void getTransitionProbabilities(final double[] Pt, final double branchLength) {
        // TODO: Pool temp matrix
        double[] calc = MatrixArrayPool.popCodonMatrix();

        // NOTE: We store matrices in row-major order: column i, row j
        for (int i = 0; i < senseCodons.length; i++) { // for each column
            double temp = Math.exp(branchLength * lambda.getQuick(i)); // get the diagonal of lambda (= column)
            for (int j = 0; j < senseCodons.length; j++) { // for each row
                calc[j * senseCodons.length + i] = U[j * senseCodons.length + i] * temp;
            }
        }

        for (int j = 0; j < senseCodons.length; j++) {
            for (int i = 0; i < senseCodons.length; i++) {
                double p = 0;
                for (int k = 0; k < senseCodons.length; k++) {
                    p += calc[i * senseCodons.length + k] * UInv[k * senseCodons.length + j]; // TODO: should UInv be stored in column-major style?
                }
                if (p < 0) p = 0;
                Pt[senseCodons[i] * GeneticCode.CODON_STATES + senseCodons[j]] = p;
            }
        }

        MatrixArrayPool.pushCodonMatrix(calc);
    }

}
