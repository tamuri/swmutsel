package swmutsel.model;

import com.beust.jcommander.internal.Lists;
import swmutsel.model.parameters.Parameter;

import java.io.Serializable;
import java.util.List;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 15:23
 */
public abstract class SubstitutionModel implements Serializable {
    private static final long serialVersionUID = 6408557902068035305L;
    private List<Parameter> parameters;

    public abstract void getTransitionProbabilities(double[] matrix, double branchlength);


    public abstract double[] getCodonFrequencies();
    public abstract void build();
    public abstract boolean parametersValid();

    public void setParameters(Parameter... parameters) {
        this.parameters = Lists.newArrayList(parameters);
    }

    public List<Parameter> getParameters() {
        return this.parameters;
    }
}
