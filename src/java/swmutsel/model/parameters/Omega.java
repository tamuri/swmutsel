package swmutsel.model.parameters;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 18/10/2013 11:16
 */
public class Omega extends Parameter {
    private static final long serialVersionUID = 8875392200058405668L;
    private double omega;

    public Omega() {
        setArgument("-omega");
    }

    public Omega(double omega) {
        this();
        this.omega = omega;
    }

    public static double getDefault() {
        return 1.0;
    }

    public double get() {
        return omega;
    }

    public double set(double omega) {
        return this.omega = omega;
    }

    @Override
    public double[] getOptimisable() {
        return new double[]{Math.log(omega)};
    }

    @Override
    public void setOptimisable(double[] params) {
        set(Math.exp(params[0]));
    }

    @Override
    public int getOptimisableCount() {
        return 1;
    }

    @Override
    public String toString() {
        return String.format("Omega{omega=%.7f}", omega);
    }


}
