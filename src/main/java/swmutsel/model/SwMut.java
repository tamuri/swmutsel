package swmutsel.model;

import swmutsel.model.parameters.BaseFrequencies;
import swmutsel.model.parameters.BranchScaling;
import swmutsel.model.parameters.Tau;
import swmutsel.model.parameters.TsTvRatio;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.Pair;

import java.io.Serializable;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 16:28
 */
public class SwMut extends SubstitutionModel implements Serializable {

    private static final long serialVersionUID = -6456754520710073776L;

    private final TsTvRatio kappa;
    public static final double DEFAULT_KAPPA = 2.0;

    private final BranchScaling c;
    public static final double DEFAULT_SCALING = 1.0;

    private final BaseFrequencies pi;
    public static final double[] DEFAULT_BASE_FREQUENCIES = new double[]{0.25, 0.25, 0.25, 0.25};

    private final Tau tau;

    private final double[] codonFrequencies = new double[GeneticCode.CODON_STATES];
    private final double[] mutationMatrix = new double[GeneticCode.CODON_STATES * GeneticCode.CODON_STATES];

    public SwMut(TsTvRatio kappa, BranchScaling c, BaseFrequencies pi) {
        this.kappa = kappa;
        this.c = c;
        this.pi = pi;
        this.tau = new Tau(Tau.getDefault());
        super.clearParameters();
        super.addParameters(this.kappa, this.pi, this.c);
        build();
    }

    public SwMut(double kappa, double scaling, double[] pi) {
        this(new TsTvRatio(kappa), new BranchScaling(scaling), new BaseFrequencies(pi));
    }

    public SwMut() {
        this(new TsTvRatio(DEFAULT_KAPPA), new BranchScaling(DEFAULT_SCALING), new BaseFrequencies(DEFAULT_BASE_FREQUENCIES));
    }

    public void setTau(double tau) {
        CoreUtils.msg("WARNING: You're using the -tau parameter!\n");
        this.tau.set(tau);
    }

    @Override
    public boolean parametersValid() {
        return true;
    }

    /*
    @Override
    public void getTransitionProbabilities(double[][] matrix, double branchlength) {
        throw new RuntimeException("ERROR: Not implemented. Use SwMutSel instead!");
    }
*/

    @Override
    public TransProbCalculator getPtCalculator() {
        throw new RuntimeException("ERROR: Not implemented!");
    }

    @Override
    public double[] getCodonFrequencies() {
        return this.codonFrequencies;
    }

    public void build() {
        setCodonFrequencies();
        setMutationRate();
    }

    public Pair<Double, Double> getRhoNonSynAndSyn() {
        double totalSynonymous = 0;
        double totalNonSynonymous = 0;

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {

                if (i == j) continue;

                int aaI = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(i);
                int aaJ = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(j);

                double rate = getCodonFrequency(i) * getMutationRate(i, j);

                // NOTE: we are including rates to/from STOP codons
                if (aaI == aaJ) {
                    totalSynonymous += rate;
                } else {
                    totalNonSynonymous += rate;
                }
            }
        }

        return Pair.of(totalNonSynonymous, totalSynonymous);
    }

    public double getMutationRate(int i, int j) {
        return mutationMatrix[i * GeneticCode.CODON_STATES + j];
    }

    // TODO: can be simplified only 576 entries need to be filled...?
    private void setMutationRate() {
        double scaling = 0.0; // Tamuri et al (2012): nu = 1 / sum_ij(mu_ij * pi_i)

        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            for (int j = 0; j < GeneticCode.CODON_STATES; j++) {

                if (i == j) {
                    mutationMatrix[i * GeneticCode.CODON_STATES + j] = 0;
                    continue;
                }

                char[] ni = GeneticCode.getInstance().getNucleotidesFromCodonIndex(i);
                char[] nj = GeneticCode.getInstance().getNucleotidesFromCodonIndex(j);

                double prodPi = 1.0;
                int changes = 0;
                int transitions = 0;


                for (int k = 0; k < 3; k++) {
                    if (ni[k] != nj[k]) {
                        changes++;
                        prodPi *= this.pi.get()[GeneticCode.getInstance().getNucleotideIndexByChar(nj[k])];
                    }
                    if (GeneticCode.getInstance().isTransitionByChar(ni[k], nj[k])) {
                        transitions++;
                    }
                }

                if (changes == 1) {
                    mutationMatrix[i * GeneticCode.CODON_STATES + j] = Math.pow(kappa.get(), transitions) * prodPi;
                    scaling += mutationMatrix[i * GeneticCode.CODON_STATES + j] * codonFrequencies[i];
                } else {
                    // No multiple substitutions allowed
                    //mutationMatrix[i * GeneticCode.CODON_STATES + j] = 0;


                /* Implementation with tau:                 */

                    mutationMatrix[i * GeneticCode.CODON_STATES + j] = Math.pow(tau.get(), changes - 1) * Math.pow(kappa.get(), transitions) * prodPi;
                    scaling += mutationMatrix[i * GeneticCode.CODON_STATES + j] * codonFrequencies[i];
                }

            }
        }

        for (int i = 0; i < mutationMatrix.length; i++) {
            mutationMatrix[i] /= scaling; // 'nu' in the 2012 paper.
        }
    }

    private void setCodonFrequencies() {
        double[] freq = this.pi.get();
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            char[] n = GeneticCode.getInstance().getCodonTLA(i).toCharArray();
            codonFrequencies[i] =
                freq[GeneticCode.getInstance().getNucleotideIndexByChar(n[0])] *
                freq[GeneticCode.getInstance().getNucleotideIndexByChar(n[1])] *
                freq[GeneticCode.getInstance().getNucleotideIndexByChar(n[2])];
        }
    }

    public double getCodonFrequency(int codon) {
        return codonFrequencies[codon];
    }

    public TsTvRatio getKappa() {
        return kappa;
    }

    public BranchScaling getBranchScaling() {
        return c;
    }

    public BaseFrequencies getPi() {
        return pi;
    }

    @Override
    public String toString() {
        return String.format("%s %.7f %s %.7f,%.7f,%.7f,%.7f %s %.7f",
                kappa.getArgument(), kappa.get(),
                pi.getArgument(), pi.get()[0], pi.get()[1], pi.get()[2], pi.get()[3],
                c.getArgument(), c.get());
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        SwMut swMut = (SwMut) o;

        if (!c.equals(swMut.c)) return false;
        if (!kappa.equals(swMut.kappa)) return false;
        if (!pi.equals(swMut.pi)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = kappa.hashCode();
        result = 31 * result + c.hashCode();
        result = 31 * result + pi.hashCode();
        return result;
    }
}

