package swmutsel.model;

import cern.colt.matrix.DoubleMatrix2D;
import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.Omega;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;

import java.io.Serializable;

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
    private final int numSenseCodons;
    
    private final double[] senseCodonFrequencies;

    private final double[] lambda; // eigenvalue of B
    private final double[][] U;
    private final double[][] UInv;

    private volatile transient TransProbCalculator calculator;

    public FMutSel0(double kappa, double omega, double[] pi, double[] fitness) {
        this(new TsTvRatio(kappa), new Omega(omega), new BaseFrequencies(pi), new Fitness(fitness));
    }

    public void setScale(double scale) {
        this.scale = scale;
    }

    private double scale = -1;

    public FMutSel0(TsTvRatio kappa, Omega omega, BaseFrequencies pi, Fitness fitness) {
        this.kappa = kappa;
        this.omega = omega;
        this.pi = pi;
        this.fitness = fitness;

        super.clearParameters();
        super.addParameters(pi, kappa, omega, fitness);

        this.senseCodons = GeneticCode.getInstance().getSenseCodons();
        this.numSenseCodons = this.senseCodons.length;
        
        this.senseCodonFrequencies = new double[this.numSenseCodons];

        this.U = new double[this.numSenseCodons][this.numSenseCodons];
        this.UInv = new double[this.numSenseCodons][this.numSenseCodons];
        this.lambda = new double[this.numSenseCodons];

        build();
    }

    @Override
    public void build() {
        setSenseCodonFrequencies();
        setEigenFactors(makeB(makeQ(), numSenseCodons, senseCodonFrequencies));
        // this.calculator = TransProbCalcFactory.getDecompositionCalculator(lambda, U, UInv, senseCodons, true);
        this.calculator = TransProbCalcFactory.getDecompositionCalculator(lambda, U, UInv, senseCodons, false);
    }

    @Override
    public boolean parametersValid() {

        for (int i = 0; i < this.fitness.get().length; i++) {
            if (this.fitness.get()[i] > 29 || this.fitness.get()[i] < -29) return false;
        }

        return true;
    }

    @Override
    public TransProbCalculator getPtCalculator() {
        /*if (changed || calculator == null) {
            synchronized (this) {
                if (changed || calculator == null) {
                    this.calculator = TransProbCalcFactory.getDecompositionCalculator(lambda, U, UInv, senseCodons, true);
                    changed = false;
                }
            }
        }*/
        return this.calculator;
    }

    public void setSenseCodonFrequencies() {
        double sum = 0;

        for (int i = 0; i < this.numSenseCodons; i++) {
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


    public double[] makeQ() {
        double[] Q = new double[this.numSenseCodons * this.numSenseCodons];

        // Fill off-diagonal entries
        for (int i = 0; i < numSenseCodons; i++) {
            col: for (int j = 0; j < numSenseCodons; j++) {

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

                Q[i * numSenseCodons + j] = rate;
            }
        }

        double sum = 0;
        for (int i = 0; i < numSenseCodons; i++)
            for (int j = 0; j < numSenseCodons; j++)
                if (i != j)
                    sum += Q[i * numSenseCodons + j] * senseCodonFrequencies[i];

        if (scale == -1) {
            // Normalise so average rate is 1
            for (int i = 0; i < Q.length; i++) Q[i] = Q[i] / sum;
        } else if (scale > 0) {
            for (int i = 0; i < Q.length; i++) Q[i] = Q[i] / scale;
        } else {
            // don't normalise at all! (i.e. scale = 0)
        }

        // Each row sums to 0
        for (int i = 0; i < numSenseCodons; i++) {
            sum = 0;
            for (int j = 0; j < numSenseCodons; j++) {
                sum += Q[i * numSenseCodons + j];
            }
            Q[i * numSenseCodons + i] = -sum;
        }

        return Q;
    }

    public double getAverageRate() {
        double[] Q = makeQ();
        double sum = 0;
        for (int i = 0; i < numSenseCodons; i++)
            for (int j = 0; j < numSenseCodons; j++)
                if (i != j)
                    sum += Q[i * numSenseCodons + j] * senseCodonFrequencies[i];

        return sum;


    }


    protected void setEigenFactors(DoubleMatrix2D B) {
        setEigenFactors(B, this.lambda, this.numSenseCodons, this.senseCodonFrequencies, this.U, this.UInv);
    }

    @Override
    public double[] getCodonFrequencies() {
        // return frequencies for all 64 codons - compare to setSenseCodonFrequencies()
        double[] frequencies = new double[GeneticCode.CODON_STATES];
        for (int i = 0; i < numSenseCodons; i++) {
            frequencies[senseCodons[i]] = senseCodonFrequencies[i];
        }
        return frequencies;
    }

    @Override
    public String toString() {
        return String.format("FMutSel0{ -kappa %.7f -omega %.7f -pi %s -fitness %s }",
                kappa.get(),
                omega.get(),
                CoreUtils.join("%.7f", ",", pi.get()),
                CoreUtils.join("%.7f", ",", fitness.get()));
    }

    public TsTvRatio getKappa() {
        return kappa;
    }

    public Omega getOmega() {
        return omega;
    }

    public BaseFrequencies getPi() {
        return pi;
    }
}
