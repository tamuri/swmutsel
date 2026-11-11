package swmutsel.model;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import com.google.common.primitives.Ints;
import swmutsel.Constants;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessFDS;
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
    private final boolean hasFitnessFDS;
    private final FitnessFDS fitnessFDS;

    private final int[] senseCodons;
    private final int numSenseCodons;

    private final DoubleMatrix2D Bcodon;
    private final DoubleMatrix2D tQ;

    private final double[] senseCodonFrequencies;
    private final double[] lambda;
    private final double[][] U;
    private final double[][] Uinv;
    private double[] Q;

    private final double[] outFrequencies = new double[GeneticCode.CODON_STATES];

    private volatile transient boolean changed = true;
    private volatile transient TransProbCalculator calculator;

    private double expectedSubstitutions;

    public SwMutSel(SwMut mutation, Fitness fitness) {
        this(mutation, fitness, null);
    }

    public SwMutSel(SwMut mutation, Fitness fitness, FitnessFDS fitnessFDS) {
        this.mutation = mutation;
        this.fitness = fitness;

        super.clearParameters();

        if (fitnessFDS == null) {
            this.hasFitnessFDS = false;
            this.fitnessFDS = null;
            super.addParameters(fitness);
        } else {
            this.hasFitnessFDS = true;
            this.fitnessFDS = fitnessFDS;
            super.addParameters(fitness, fitnessFDS);
        }

//        super.addParameters(fitness, fitnessFDS); // We're assuming mutation component is fixed here

        this.senseCodons = GeneticCode.getInstance().getSenseCodons();
        this.numSenseCodons = this.senseCodons.length;
        this.Bcodon = DoubleFactory2D.dense.make(numSenseCodons, numSenseCodons);
        this.tQ = DoubleFactory2D.dense.make(numSenseCodons, numSenseCodons);

        this.senseCodonFrequencies = new double[this.numSenseCodons];
        this.U = new double[this.numSenseCodons][this.numSenseCodons];
        this.Uinv = new double[this.numSenseCodons][this.numSenseCodons];
        this.lambda = new double[this.numSenseCodons];

        build();
    }


    @Override
    final public void build() {
        this.Q = makeQ();
        setSenseCodonFrequencies();
        if (!hasFitnessFDS || fitnessFDS.getModel() == 2) {
            setEigenFactors(makeB(this.Q, numSenseCodons, senseCodonFrequencies));
        }
        //setEigenFactorsQ(this.Q);
        changed = true;
    }

    private void getCallingMethodName() {

        StackTraceElement[] f = Thread.currentThread().getStackTrace();
        for (int i = 2; i < f.length; i++) {
            if (i > 2) System.out.print("\t");
            System.out.println(f[i].toString());
        }

    }

    private void setEigenFactorsQ(double[] Q) {
        DoubleMatrix2D qq = DoubleFactory2D.dense.make(Q, numSenseCodons);

        for (int i = 0; i < numSenseCodons; i++) {
            for (int j = 0; j < numSenseCodons; j++) {
                qq.setQuick(i, j, Q[i * numSenseCodons + j]);
            }
        }

        EigenvalueDecomposition eigen = new EigenvalueDecomposition(qq);


//getCallingMethodName();
//        System.exit(0);

        DoubleMatrix1D I = eigen.getImagEigenvalues();
        if (I.zSum() > 0) {
            System.out.println(I.toString());
            System.exit(1);
        }

        DoubleMatrix1D l = eigen.getRealEigenvalues();
        for (int i = 0; i < numSenseCodons; i++) {
            lambda[i] = l.get(i) * mutation.getBranchScaling().get();
        }

        DoubleMatrix2D evec = eigen.getV();
        DoubleMatrix2D ievec = Algebra.DEFAULT.inverse(evec);

        for (int i = 0; i < numSenseCodons; i++) {
            for (int j = 0; j < numSenseCodons; j++) {
                U[i][j] = evec.getQuick(i, j);
                Uinv[i][j] = ievec.getQuick(i, j);
            }
        }
    }

    @Override
    final public boolean parametersValid() {
        return true;
    }


    private void setSenseCodonFrequencies() {
        if (hasFitnessFDS && fitnessFDS.getModel() == 1) {
            for (int i = 0; i < numSenseCodons; i++) {
                for (int j = 0; j < numSenseCodons; j++) {
                    tQ.set(i, j, Q[j * numSenseCodons + i]);
                }
            }
            for (int i = 0; i < numSenseCodons; i++) {
                tQ.setQuick(0, i, 1.0);
            }
            Bcodon.setQuick(0, 0, 1.0);
            DoubleMatrix2D out = Algebra.DEFAULT.solve(tQ, Bcodon);
            for (int i = 0; i < numSenseCodons; i++) {
                senseCodonFrequencies[i] = out.getQuick(i, 0);
            }
        } else {
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
        // Q = new double[numSenseCodons * numSenseCodons];

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

                double hS;
                if (Ai == Aj) {
                    hS = getRelativeFixationProbability(Fj - Fi);
                } else {
                    //hS = getRelativeFixationProbability(Fj - Fi + (hasFitnessFDS ? fitnessFDS.get() : 0));
                    if (hasFitnessFDS) {
                        if (fitnessFDS.getModel() == 1) {
                            hS = getRelativeFixationProbability(Fj - Fi + fitnessFDS.get());
                        } else if (fitnessFDS.getModel() == 2) {
                            double Z = fitnessFDS.get();
                            if (Z == 0) {
                                hS = 1;
                            } else {
                                hS = Z / (1 - Math.exp(-Z));
                            }
                        } else {
                            throw new RuntimeException("unknown FDS model");
                        }
                    } else {
                        hS = getRelativeFixationProbability(Fj - Fi);
                    }
                }

                Q[i * numSenseCodons + j] = mutationRate * hS;

            }
        }


        // Each row sums to 0
        for (int row = 0; row < numSenseCodons; row++) {
            double rowSum = 0;
            for (int column = 0; column < numSenseCodons; column++) {
                rowSum += Q[row * numSenseCodons + column];
            }
            Q[row * numSenseCodons + row] = -rowSum;
        }


        return Q;
    }

    public final double getExpectedSubsPerSite() {
        double expectedSubstitutions = 0;
        for (int i = 0; i < numSenseCodons; i++) {
            expectedSubstitutions += -Q[i * numSenseCodons + i] * senseCodonFrequencies[i];
        }
        return expectedSubstitutions;
    }

    public final Pair<Double, Double> getRhoNonSynAndSyn() {
        double totalSynonymous = 0;
        double totalNonSynonymous = 0;
        for (int i = 0; i < numSenseCodons; i++) {
            for (int j = 0; j < numSenseCodons; j++) {
                if (i == j) {
                    continue;
                }
                double rate = senseCodonFrequencies[i] * Q[i * numSenseCodons + j];
                int Ai = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(senseCodons[i]);
                int Aj = GeneticCode.getInstance().getAminoAcidIndexFromCodonIndex(senseCodons[j]);
                if (Ai == Aj) {
                    totalSynonymous += rate;
                } else {
                    totalNonSynonymous += rate;
                }
            }
        }
        return Pair.of(totalNonSynonymous, totalSynonymous);
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
        // double[] frequencies = new double[GeneticCode.CODON_STATES];
        for (int i = 0; i < numSenseCodons; i++) {
            outFrequencies[senseCodons[i]] = senseCodonFrequencies[i];
        }
        return outFrequencies;
    }

    @Override
    public TransProbCalculator getPtCalculator() {
        if (changed || calculator == null) {
            synchronized (this) {
                if (changed || calculator == null) {
                    if (hasFitnessFDS && fitnessFDS.getModel() == 1) {
                        this.calculator = TransProbCalcFactory.getPadeCalculator(Q, senseCodons, mutation.getBranchScaling().get());
                    } else {
                        this.calculator = TransProbCalcFactory.getDecompositionCalculator(lambda, U, Uinv, senseCodons, false);
                    }
                    changed = false;
                }
            }
        }
        return this.calculator;
    }

    public boolean hasFitnessFDS() {
        return hasFitnessFDS;
    }

    private boolean doPenalty = false;

    public void setDoPenalty(boolean doPenalty) {
        this.doPenalty = doPenalty;
    }

    public void setPenaltyShape(double alpha, double beta) {
        this.alpha = alpha;
        this.beta = beta;
    }

    private double alpha = 1.0;
    private double beta = 1.0;

    public double getFDSPenalty() {
        if (doPenalty) {
            double parameter = fitnessFDS.get();
            if (fitnessFDS.getModel() == 1) {
                // return Math.log(Math.pow(parameter, alpha - 1) * Math.exp(-beta * parameter));
                return -beta * parameter;
            } else if (fitnessFDS.getModel() == 2) {
                double theta = Math.exp(parameter) / (1 + Math.exp(parameter));
                return Math.log(Math.pow(theta, 0.01 - 1.0) * Math.pow(1 - theta, 0.01 - 1.0));
            }
            return 0;
        } else {
            return 0;
        }
    }

}
