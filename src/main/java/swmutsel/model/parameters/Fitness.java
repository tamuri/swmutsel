package swmutsel.model.parameters;

import swmutsel.Constants;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 22:28
 */
public class Fitness extends Parameter {

    private static final long serialVersionUID = 1479750983891431808L;
    private double[] fitness;
    private double[] fitnessForOptimiser = new double[GeneticCode.AMINO_ACID_STATES - 1]; // One amino acid is fixed to 0
    private int[] fitnessOrderForOptimiser = new int[GeneticCode.AMINO_ACID_STATES];

    public Fitness() {
        this(CoreUtils.seq(0, 19));
    }

    public Fitness(double[] fitness) {
        this(CoreUtils.seq(0, 19), fitness);
    }

    public Fitness(int[] order) {
        this(order, new double[20]);
    }

    public Fitness(int[] order, double[] fitness) {
        setArgument("-fitness");
        this.fitnessOrderForOptimiser = order;
        set(fitness);
    }

    public Fitness copy() {
        return new Fitness(fitnessOrderForOptimiser.clone(), fitness.clone());
    }

    public double[] get() {
        return fitness;
    }

    public void set(double[] fitness) {
        this.fitness = fitness;
    }

    @Override
    public double[] getOptimisable() {
        int pos = 0;
        for (int i = 1; i < GeneticCode.AMINO_ACID_STATES; i++) {
            this.fitnessForOptimiser[pos] = this.fitness[this.fitnessOrderForOptimiser[i]];
            pos++;
        }

        return this.fitnessForOptimiser;
    }

    @Override
    public void setOptimisable(double[] params) {
        this.fitnessForOptimiser = params;
        this.fitness[this.fitnessOrderForOptimiser[0]] = Constants.FITNESS_FIXED_FOR_RELATIVE;

        int pos = 0;
        for (int i = 1; i < 20; i++) {
            this.fitness[this.fitnessOrderForOptimiser[i]] = params[pos];
            pos++;
        }
    }

    @Override
    public int getOptimisableCount() {
        return 19;
    }

    @Override
    public String toString() {
        return "Fitness{" +
                "fitness=" + CoreUtils.join("%.7f", ", ", fitness) +
                '}';
    }
}
