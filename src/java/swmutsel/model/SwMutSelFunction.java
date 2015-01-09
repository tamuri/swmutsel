package swmutsel.model;

import com.google.common.collect.Lists;
import org.apache.commons.math3.analysis.MultivariateFunction;
import swmutsel.Constants;
import swmutsel.model.parameters.Mapper;
import swmutsel.model.parameters.Parameter;

import java.util.List;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 31/10/2013 15:21
 */
public class SwMutSelFunction implements MultivariateFunction {

    private final LikelihoodCalculator calculator;
    private final List<SubstitutionModel> models;
    private final List<Parameter> parameters;

    public SwMutSelFunction(LikelihoodCalculator calculator, SubstitutionModel model) {
        this(calculator, Lists.newArrayList(model));
    }

    public SwMutSelFunction(LikelihoodCalculator calculator, List<SubstitutionModel> models) {
        this.calculator = calculator;
        this.models = models;

        this.parameters = Lists.newArrayList();

        for (SubstitutionModel sm : models) {
            this.parameters.addAll(sm.getParameters());
        }
    }

    public double[] getCurrentParameters() {
        return Mapper.getOptimisable(this.parameters);
    }

    @Override
    public double value(double[] point) {

        // todo: looks like we can ignore bounds when using lbfgs
        for (double d : point) {
            if (d < -(Constants.FITNESS_BOUND + 1) || d > (Constants.FITNESS_BOUND + 1)) {
                return Constants.VERY_BAD_LIKELIHOOD;
            }
        }

        Mapper.setOptimisable(parameters, point);
        for (SubstitutionModel sm : models) sm.build();

        return calculator.getLogLikelihood();

    }
}
