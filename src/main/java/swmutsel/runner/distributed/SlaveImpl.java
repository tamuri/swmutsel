package swmutsel.runner.distributed;

import com.google.common.collect.Table;
import com.google.common.primitives.Ints;
import pal.tree.Tree;
import swmutsel.ArrayPool;
import swmutsel.model.Penalty;
import swmutsel.model.SubstitutionModel;
import swmutsel.model.SwMut;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.runner.MultiThreadedRunner;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GuavaUtils;
import swmutsel.utils.Pair;

import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 30/10/2013 11:34
 */
public class SlaveImpl implements SlaveAPI  {
    private final Table<String, Integer, Byte> allSites;
    private final MultiThreadedRunner runner;

    public SlaveImpl(Table<String, Integer, Byte> allSites, int threads) {
        this.allSites = allSites;
        this.runner = new MultiThreadedRunner(threads);
    }

    @Override
    public double getLogLikelihoodForBranch(int node, double branchLength) {
        return this.runner.getLogLikelihoodForBranch(node, branchLength);
    }

    @Override
    public void setBranchLength(int node, double branchLength) {
        this.runner.setBranchLength(node, branchLength);
    }

    @Override
    public Map<Integer, Double> getLogLikelihood(Tree tree, SwMut mutation, FitnessStore fitnesses, Penalty penalty, boolean save) {
        return this.runner.getLogLikelihood(tree, mutation, fitnesses, penalty, save);
    }

    @Override
    public Pair<Map<Integer,Double>, List<FitnessStore>> optimiseFitness(Tree tree, SwMut mutation, FitnessStore fitnesses, Penalty penalty, List<String> cladeModels, int numberOfOptimRestarts) {
        return this.runner.optimiseFitness(tree, mutation, fitnesses, penalty, cladeModels, numberOfOptimRestarts);
    }

    @Override
    public void setSites(int[] sites) {
        Table<String, Integer, Byte> siteTable = GuavaUtils.getColumnSubset(allSites, Ints.asList(sites));
        this.runner.setSites(siteTable);
/*

        // TODO: this should be in a better place!!!
        // In unrooted binary tree, there are n-2 internal nodes
        ArrayPool.setTreeSize(allSites.rowKeySet().size());
*/

        CoreUtils.msg("Local runner sites set. Sites are %s\n", this.runner.getSites().columnKeySet());
    }

    @Override
    public void setRunnerPenalty(Penalty penalty) {
        this.runner.setRunnerPenalty(penalty);
    }

    @Override
    public void setRunnerMutation(SwMut mutation) {
        this.runner.setRunnerMutation(mutation);
    }

    @Override
    public void setRunnerTree(Tree tree) {
        ArrayPool.setTreeSize(tree.getExternalNodeCount());
        this.runner.setRunnerTree(tree);
    }

    @Override
    public void setRunnerFitnesses(FitnessStore fitnesses) {
        this.runner.setRunnerFitnesses(fitnesses);
    }

    @Override
    public void shutdown() {
        this.runner.shutdown();
        Slave.shutdown();
    }

    @Override
    public void setPatterns(Map<Integer, Integer> patternSiteMap, Map<Integer, Integer> patternWeight) {
        this.runner.setPatterns(patternSiteMap, patternWeight);
    }

    @Override
    public Map<Integer, Double> getLogLikelihoodForTree(final Tree tree, final boolean save) {
        return this.runner.getLogLikelihood(tree, save);
    }

    @Override
    public Map<Integer, Double> getLogLikelihoodForMutation(final SwMut mutation, final boolean save) {
        return this.runner.getLogLikelihood(mutation, save);
    }

    @Override
    public Map<Integer, Double> getLogLikelihoodForSubsmodel(final SubstitutionModel model, final boolean save) {
        model.build();
        return this.runner.getLogLikelihood(model, save);
    }

    @Override
    public double updateTree(final Tree t, final Map<Integer, Integer> first, final Set<Integer> second) {
        return this.runner.updateTree(t, first, second);
    }

}
