package swmutsel.model.parameters;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 21:38
 */
public class BranchScaling extends Parameter {
    private static final long serialVersionUID = -703874549036930982L;
    private double c;

    public BranchScaling() {
        setArgument("-scaling");
    }

    public BranchScaling(double scaling) {
        this();
        this.c = scaling;
    }

    public static double getDefault() {
        return 1.0;
    }

    public double get() {
        return c;
    }

    public void set(double c) {
        this.c = c;
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


    @Override
    public String toString() {
        return "BranchScaling{" +
                "c=" + c +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        BranchScaling that = (BranchScaling) o;

        if (Double.compare(that.c, c) != 0) return false;

        return true;
    }

    @Override
    public int hashCode() {
        long temp = Double.doubleToLongBits(c);
        return (int) (temp ^ (temp >>> 32));
    }
}
