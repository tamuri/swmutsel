package swmutsel.model.parameters;

import swmutsel.Run;

import java.util.Objects;

/**
 * A fitness parameter used for frequency-dependent selection
 *
 * Created by atamuri on 29/06/15.
 */
public class FitnessFDS extends Parameter {
    private double fitness;
    private int model;

    public FitnessFDS() {
        setArgument("-fitnessfds");
    }

    public int getModel() {
        return this.model;
    }

    public FitnessFDS(double fitness, int model) {
        this();
        this.fitness = fitness;
        this.model = model;
    }

    public static double getDefault() { return 0.0; }

    @Override
    public double[] getOptimisable() {
        return new double[]{get()};
    }

    @Override
    public void setOptimisable(double[] params) {
        if (model == 1) {
            // this is swMutSel-F Zk, positive only
            //set(Math.pow(params[0], 2));
            set(Math.abs(params[0]));
        } else if (model == 2) {
            // this is NY03
            set(params[0]);
        } else {
            throw new RuntimeException("FDS model not known");
        }
    }

    private void set(double fitness) {
        this.fitness = fitness;
    }

    public double get() {
        return this.fitness;
    }

    @Override
    public int getOptimisableCount() {
        return 1;
    }

    @Override
    public String toString() {
        return "FitnessFDS{" +
                "fitness=" + fitness +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        FitnessFDS that = (FitnessFDS) o;
        return Objects.equals(fitness, that.fitness);
    }

    @Override
    public int hashCode() {
        return Objects.hash(fitness);
    }
}
