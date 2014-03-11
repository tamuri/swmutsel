package swmutsel.runner.distributed;

import pal.tree.Tree;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.model.Penalty;
import swmutsel.model.SwMut;
import swmutsel.utils.Pair;

import java.util.List;
import java.util.Map;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 30/10/2013 11:34
 */
public interface SlaveAPI {
    public void setSites(int[] sites);

    public void setRunnerPenalty(Penalty penalty);

    public void setRunnerMutation(SwMut mutation);

    public void setRunnerTree(Tree tree);

    public void setRunnerFitnesses(FitnessStore fitnesses);

    public Map<Integer, Double> getLogLikelihood(Tree tree, SwMut mutation, FitnessStore fitnesses, Penalty penalty, boolean save);

    public double getLogLikelihoodForBranch(int node, double branchLength);

    public Pair<Map<Integer,Double>, List<FitnessStore>> optimiseFitness(Tree tree, SwMut mutation, FitnessStore fitnesses, Penalty penalty, List<String> cladeModel, int numberOfOptimRestarts);

    public void setBranchLength(int node, double branchLength);

    public void shutdown();

    public void setPatterns(Map<Integer,Integer> patternSiteMap, Map<Integer,Integer> patternWeight);

    public Map<Integer,Double> getLogLikelihoodForTree(Tree tree, boolean save);

    public Map<Integer,Double> getLogLikelihoodForMutation(SwMut mutation, boolean save);
}
