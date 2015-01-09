package swmutsel.model.parameters;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 21:27
 */
public class TsTvRatio extends Parameter {

    private static final long serialVersionUID = -5812168500083535060L;
    private double kappa;

    public TsTvRatio() {
        setArgument("-kappa");
    }

    public TsTvRatio(double kappa) {
        this();
        this.kappa = kappa;
    }

    public static double getDefault() {
        return 1.0;
    }

    @Override
    public double[] getOptimisable() {
        return new double[]{Math.log(get())};
    }

    @Override
    public void setOptimisable(double[] params) {
        set(Math.exp(params[0]));
    }

    @Override
    public int getOptimisableCount() {
        return 1;
    }

    public double get() {
        return this.kappa;
    }

    public void set(double kappa) {
        this.kappa = kappa;
    }

    @Override
    public String toString() {
        return "TsTvRatio{" +
                "kappa=" + kappa +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        TsTvRatio tsTvRatio = (TsTvRatio) o;

        if (Double.compare(tsTvRatio.kappa, kappa) != 0) return false;

        return true;
    }

    @Override
    public int hashCode() {
        long temp = Double.doubleToLongBits(kappa);
        return (int) (temp ^ (temp >>> 32));
    }
}
