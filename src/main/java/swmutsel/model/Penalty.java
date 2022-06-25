package swmutsel.model;

import java.io.Serializable;

public interface Penalty extends Serializable {
    long serialVersionUID = -4158732377735288388L;
    double calculate(final double[] parameters);
}
