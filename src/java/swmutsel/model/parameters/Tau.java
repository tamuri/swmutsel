package swmutsel.model.parameters;

/**
 * NOTE: Tau is fixed to 0
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 23:18
 */
public class Tau extends Parameter {
    private static final long serialVersionUID = 2332885343496674578L;
    private final double tau = 0.0;

    public Tau() {
        setArgument("-tau");
    }

    public double get() {
        return tau;
    }

    public void set() {
        // do nothing;
    }

    @Override
    public double[] getOptimisable() {
        return null;
    }

    @Override
    public void setOptimisable(double[] params) {
        // do nothing
    }

    @Override
    public int getOptimisableCount() {
        return 0;
    }
}
