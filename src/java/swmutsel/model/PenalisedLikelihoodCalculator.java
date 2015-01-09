package swmutsel.model;

import pal.tree.Tree;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.Parameter;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;

import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;

/**
 * Wraps LikelihoodCalculator and adds penalty on Fitness parameters to log-likelihood.
 *
 * TODO: possibly have a list of penalties and penalised parameters...?
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 18/10/2013 21:22
 */
public class PenalisedLikelihoodCalculator extends LikelihoodCalculator {
    private final Parameter penalisedParameter;
    private final Penalty penalty;

    public PenalisedLikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, LinkedHashMap<String, SubstitutionModel> models, Penalty penalty, Parameter penalisedParameter) {
        super(tree, siteStates, models);
        this.penalisedParameter = penalisedParameter;
        this.penalty = penalty;
    }

    @SuppressWarnings("unchecked")
    public PenalisedLikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, SubstitutionModel model, Penalty penalty, Parameter penalisedParameter) {
        this(tree, siteStates, CoreUtils.getLinkedHashMap(Pair.of("ALL", model)), penalty, penalisedParameter);
    }

    @Override
    public double getLogLikelihood() {
        return super.getLogLikelihood() + calculatePenalty();
    }

    @Override
    public double updateTree(Tree tree, Map<Integer, Integer> newToOldNumbering, Set<Integer> updateNodes) {
        return super.updateTree(tree, newToOldNumbering, updateNodes) + calculatePenalty();
    }

    @Override
    public double getLogLikelihoodForBranch(int node, double branchLength) {
        double out = super.getLogLikelihoodForBranch(node, branchLength);
        out += calculatePenalty();
        return out;
    }

    private double calculatePenalty() {
        double penalty = 0.0;
        if (this.penalisedParameter instanceof Fitness) {
            penalty += this.penalty.calculate(penalisedParameter.getOptimisable());
        }
        return penalty;
    }
}
