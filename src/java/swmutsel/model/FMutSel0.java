package swmutsel.model;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import swmutsel.MatrixArrayPool;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.utils.GeneticCode;

import java.io.Serializable;
import java.util.concurrent.ExecutionException;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 15:24
 */
public class FMutSel0 extends SubstitutionModel implements Serializable {

    private static final long serialVersionUID = 7703558036907403610L;

    // The parameters of the model
    private final TsTvRatio kappa;
    private final Omega omega;
    private final BaseFrequencies pi;
    private final Fitness fitness;

    private final int[] senseCodons;
    private final double[] senseCodonFrequencies;

    private DoubleMatrix1D lambda; // eigenvalue of B
    private final double[] U;
    private final double[] UInv;

    private final LoadingCache<Double, double[]> cachedTransitionP;

    public FMutSel0(double kappa, double omega, double[] pi, double[] fitness) {
        this(new TsTvRatio(kappa), new Omega(omega), new BaseFrequencies(pi), new Fitness(fitness));
    }

    public FMutSel0(TsTvRatio kappa, Omega omega, BaseFrequencies pi, Fitness fitness) {
        this.kappa = kappa;
        this.omega = omega;
        this.pi = pi;
        this.fitness = fitness;
        setParameters(pi, kappa, omega, fitness);

        this.senseCodons = GeneticCode.getInstance().getSenseCodons();
        this.senseCodonFrequencies = new double[this.senseCodons.length];

        this.U = new double[this.senseCodons.length * this.senseCodons.length];
        this.UInv = new double[this.senseCodons.length * this.senseCodons.length];

        this.cachedTransitionP = CacheBuilder.newBuilder().build(new CachedPtLoader());

        build();
    }

    @Override
    public void build() {
        setSenseCodonFrequencies();
        setEigenFactors(makeB(makeQ()));
        this.cachedTransitionP.invalidateAll();
    }

    @Override
    public boolean parametersValid() {

        for (int i = 0; i < this.fitness.get().length; i++) {
            if (this.fitness.get()[i] > 29 || this.fitness.get()[i] < -29) return false;
        }

        return (this.kappa.get() > 0 && this.kappa.get() < 10) && (this.omega.get() > 0 && this.omega.get() < 10);
    }

    @Override
    public void getTransitionProbabilities(double[] Pt, double branchLength) {
        try {
            double[] temp = cachedTransitionP.get(branchLength);
            System.arraycopy(temp, 0, Pt, 0, Pt.length);
        } catch (ExecutionException e) {
            throw new RuntimeException(e);
        }
    }

    private void setSenseCodonFrequencies() {
        double sum = 0;

        for (int i = 0; i < this.senseCodons.length; i++) {
            int codon = senseCodons[i];
            char[] nuc = GeneticCode.getInstance().getNucleotidesFromCodonIndex(codon);

            double prod = 1;
            for (char n : nuc) {
                prod *= this.pi.get()[GeneticCode.getInstance().getNucleotideIndexByChar(n)];
            }

            senseCodonFrequencies[i] = prod * Math.exp(this.fitness.get()[GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(codon)]);
            sum += senseCodonFrequencies[i];
        }

        for (int i = 0; i < senseCodonFrequencies.length; i++) {
            senseCodonFrequencies[i] /= sum;
        }
    }


    private double[] makeQ() {
        double[] Q = new double[this.senseCodons.length * this.senseCodons.length];

        // Fill off-diagonal entries
        for (int i = 0; i < senseCodons.length; i++) {
            col: for (int j = 0; j < senseCodons.length; j++) {

                if (i == j) continue;

                int codonI = senseCodons[i];
                int codonJ = senseCodons[j];

                int changes = 0;
                boolean transition = false;
                double piTo = 0;

                char[] ni = GeneticCode.getInstance().getNucleotidesFromCodonIndex(codonI);
                char[] nj = GeneticCode.getInstance().getNucleotidesFromCodonIndex(codonJ);

                for (int k = 0; k < 3; k++) {
                    if (ni[k] != nj[k]) {
                        changes++;

                        if (changes > 1) continue col; // don't allow multiple substitutions

                        piTo = pi.get()[GeneticCode.getInstance().getNucleotideIndexByChar(nj[k])];

                        if (GeneticCode.getInstance().isTransitionByChar(ni[k], nj[k])) {
                            transition = true;
                        }
                    }
                }

                double rate = piTo;

                if (transition) rate *= kappa.get();

                int Ai = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(codonI);
                int Aj = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(codonJ);

                // If non-synonymous change
                if (Ai != Aj) rate *= omega.get();

                double Fi = this.fitness.get()[Ai];
                double Fj = this.fitness.get()[Aj];

                double hS = getRelativeFixationProbability(Fj - Fi);

                rate *= hS;

                Q[i * senseCodons.length + j] = rate;
            }
        }

        // Normalise so average rate is 1
        double sum = 0;
        for (int i = 0; i < senseCodons.length; i++) {
            for (int j = 0; j < senseCodons.length; j++) {
                if (i == j ) continue;
                sum += Q[i * senseCodons.length + j] * senseCodonFrequencies[i];
            }
        }

        for (int i = 0; i < Q.length; i++) Q[i]  = Q[i] / sum;

        // Each row sums to 0
        for (int i = 0; i < senseCodons.length; i++) {
            sum = 0;
            for (int j = 0; j < senseCodons.length; j++) {
                sum += Q[i * senseCodons.length + j];
            }
            Q[i * senseCodons.length + i] = -sum;
        }

        return Q;
    }

    private double getRelativeFixationProbability(double S) {
        return (S == 0) ? 1 : S / (1 - Math.exp(-S));
    }

    private DoubleMatrix2D makeB(double[] Q) {
        DoubleMatrix2D B = DoubleFactory2D.dense.make(senseCodons.length, senseCodons.length);
        for (int i = 0; i < senseCodons.length; i++) {
            for (int j = 0; j < senseCodons.length; j++) {
                B.setQuick(i, j, Q[i * senseCodons.length + j] * Math.sqrt(this.senseCodonFrequencies[i]) / Math.sqrt(this.senseCodonFrequencies[j]));
            }
        }

        return B;
    }

    private void setEigenFactors(DoubleMatrix2D B) {
        EigenvalueDecomposition evdB = new EigenvalueDecomposition(B);
        lambda = DoubleFactory2D.dense.diagonal(evdB.getD());

        DoubleMatrix2D R = evdB.getV();

        for (int i = 0; i < senseCodons.length; i++) {
            double piSqrt = Math.sqrt(this.senseCodonFrequencies[i]);
            double piInvSqrt = 1 / piSqrt;
            for (int j = 0; j < senseCodons.length; j++) {
                U[i * senseCodons.length + j] = piInvSqrt * R.getQuick(i, j);
                UInv[j * senseCodons.length + i] = piSqrt * R.getQuick(i, j); // inverse(R) == transpose(R)
            }
        }
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

    private class CachedPtLoader extends CacheLoader<Double, double[]> {
        @Override
        public double[] load(Double branchLength) throws Exception {
            double[] Pt = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

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
                        p += calc[i * senseCodons.length + k] * UInv[k * senseCodons.length + j];
                    }
                    if (p < 0) p = 0;
                    Pt[senseCodons[i] * GeneticCode.CODON_STATES + senseCodons[j]] = p;
                }
            }
            return Pt;
        }
    }

}
