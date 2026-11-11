package swmutsel.model.parameters;

import java.io.Serializable;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 21:25
 */
public abstract class Parameter implements Serializable {
    private static final long serialVersionUID = 8578887523887828536L;

    protected String argument;
    private boolean optimise = true;

    public abstract double[] getOptimisable();
    public abstract void setOptimisable(double[] params);
    public abstract int getOptimisableCount();

    public String getArgument() {
        return argument;
    }

    public void setArgument(String argument) {
        this.argument = argument;
    }

    public boolean isOptimisable() {
        return optimise;
    }

    public void setOptimisable(boolean optimise) {
        this.optimise = optimise;
    }
}
