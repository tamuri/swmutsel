package swmutsel.model;

import java.io.Serializable;

public interface Penalty extends Serializable {
    public static final long serialVersionUID = -4158732377735288388L;

    public double calculate(final double[] parameters);
}
