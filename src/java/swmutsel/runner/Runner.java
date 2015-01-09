package swmutsel.runner;


import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.collect.Table;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.univariate.*;
import pal.tree.Node;
import pal.tree.Tree;
import swmutsel.ArrayPool;
import swmutsel.Constants;
import swmutsel.model.*;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.model.parameters.Mapper;
import swmutsel.optim.NativeLBFGSBMultivariateOptimizer;
import swmutsel.trees.RerootedTreeIterator;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;
import swmutsel.utils.PhyloUtils;

import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public abstract class Runner {
// ------------------------------ FIELDS ------------------------------

    private Table<String, Integer, Byte> sites;
    private Map<Integer, Integer> patternSiteMap;
    private Map<Integer, Integer> patternWeight;

    private Tree runnerTree;
    private FitnessStore runnerFitnesses;
    private SwMut runnerMutation;
    private Penalty runnerPenalty;

// --------------------- GETTER / SETTER METHODS ---------------------

    protected Map<Integer,Integer> getPatternSiteMap() {
        return patternSiteMap;
    }

    protected Map<Integer, Integer> getPatternWeight() {
        return patternWeight;
    }

    public FitnessStore getRunnerFitnesses() {
        return runnerFitnesses;
    }

    public void setRunnerFitnesses(FitnessStore runnerFitnesses) {
        this.runnerFitnesses = runnerFitnesses;
    }

    public SwMut getRunnerMutation() {
        return runnerMutation;
    }

    public void setRunnerMutation(SwMut runnerMutation) {
        this.runnerMutation = runnerMutation;
    }

    public Penalty getRunnerPenalty() {
        return runnerPenalty;
    }

    public void setRunnerPenalty(Penalty runnerPenalty) {
        this.runnerPenalty = runnerPenalty;
    }

    public Tree getRunnerTree() {
        return runnerTree;
    }

    public void setRunnerTree(Tree runnerTree) {
        ArrayPool.setTreeSize(runnerTree.getExternalNodeCount());
        this.runnerTree = runnerTree;
    }

    public Table<String, Integer, Byte> getSites() {
        return sites;
    }

    public void setSites(Table<String, Integer, Byte> sites) {
        this.sites = sites;
    }

// -------------------------- OTHER METHODS --------------------------

    public double getAverageSubstitutionRate(final SwMut mutation, final FitnessStore fitnesses) {
        double lambda = 0;
        double totalSites = 0;
        for (final int site : getSites().columnKeySet()) {
            Fitness f = fitnesses.get(site);
            SwMutSel model = new SwMutSel(mutation, f);
            lambda += model.getExpectedSubsPerSite() * getPatternWeight().get(site);
            totalSites += getPatternWeight().get(site);
        }
        return lambda / totalSites;
    }

    public Map<Integer, Pair<Double, Double>> getDnDs(final SwMut mutation, final FitnessStore fitnesses) {
        Map<Integer, Pair<Double, Double>> siteDnDs = Maps.newHashMap();

        Pair<Double, Double> neutralDnDs = mutation.getRhoNonSynAndSyn();
        double neutralDn = neutralDnDs.first;
        double neutralDs = neutralDnDs.second;
        CoreUtils.msg("ϱn⁰ / ϱs⁰  = %.7f / %.7f\n", neutralDn, neutralDs);

        for (final int site : getSites().columnKeySet()) {
            Fitness f = fitnesses.get(site);
            SwMutSel model = new SwMutSel(mutation, f);
            Pair<Double, Double> dnDs = model.getRhoNonSynAndSyn();
            siteDnDs.put(site, Pair.of(dnDs.first / neutralDn, dnDs.second / neutralDs));
        }

        // For each pattern -> site mapping
        for (Map.Entry<Integer, Integer> e : getPatternSiteMap().entrySet()) {
            int pattern = e.getKey();
            int site = e.getValue();
            // If the pattern is not site AND we have the relevant site (for DistributedRunner)
            if (pattern != site) {
                siteDnDs.put(e.getKey(), siteDnDs.get(e.getValue()));
            }
        }

        return siteDnDs;
    }

    public Map<Integer, Double> getLogLikelihood(SubstitutionModel model) {
        return getLogLikelihood(model, false);
    }

    public Map<Integer, Double> getLogLikelihood(final SwMut mutation, boolean save) {
        return getLogLikelihood(getRunnerTree(), mutation, getRunnerFitnesses(), getRunnerPenalty(), save);
    }

    public Map<Integer, Double> getLogLikelihood(final Tree tree, boolean save) {
        return getLogLikelihood(tree, getRunnerMutation(), getRunnerFitnesses(), getRunnerPenalty(), save);
    }

    public Map<Integer, Double> getLogLikelihood(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty) {
        return getLogLikelihood(tree, mutation, fitnesses, penalty, false);
    }

    private Pair<Double, Tree> optimiseBranchLengths(final Tree tree, final BLOptCallback blSetup) {
        double postOptima = Double.NEGATIVE_INFINITY, preOptima = Double.NEGATIVE_INFINITY;

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);
        //ConvergenceChecker<UnivariatePointValuePair> branchConvergenceChecker = new SimpleUnivariateValueChecker(-1, Constants.MIN_BRANCH_LENGTH);
        ConvergenceChecker<UnivariatePointValuePair> branchConvergenceChecker = new SimpleUnivariateValueChecker(-1, Constants.CONVERGENCE_TOL);

        int iteration = 0;

        RerootedTreeIterator rti = new RerootedTreeIterator(tree);
        Iterator<Tree> it = rti.iterator();
        Tree startTree = it.next();

        CoreUtils.msg("Optimising tree branch lengths: ");
        double lnl = blSetup.init(startTree);
        System.out.printf("%.6f ", lnl);

        Map<Node, Pair<Integer, Set<Node>>> nodeChildMapping = Maps.newHashMap();
        PhyloUtils.saveNodeRelationships(startTree, nodeChildMapping);

        boolean converged = false;
        int smallChanges = 0;

        while (!converged) {
            // Save branches we've already optimised so we don't do them again
            Set<Set<Node>> optimisedBranches = Sets.newHashSet();

            iteration++;
            it = rti.iterator();

            // loop over every rerooted tree
            while (it.hasNext()) {
                final Tree t = it.next();

                Pair<Map<Integer,Integer>, Set<Integer>> changes = PhyloUtils.getNodeChanges(t, nodeChildMapping);
                PhyloUtils.saveNodeRelationships(t, nodeChildMapping);

                setRunnerTree(t);
                double lnL = updateTree(t, changes.first, changes.second);

                if (Double.isInfinite(preOptima)) preOptima = lnL;

                // get the root for this tree
                Node root = t.getRoot();

                // loop through each child node of the root
                for (int i = 0; i < root.getChildCount(); i++) {
                    final Node child = root.getChild(i);

                    Set<Node> branch = Sets.newHashSet(root, child);

                    // If we haven't optimised this branch yet
                    if (!optimisedBranches.contains(branch)) {
                        optimisedBranches.add(branch);

                        double old = child.getBranchLength();

                        UnivariateOptimizer opt;
                        if (Constants.CONVERGENCE_TOL > 5e-6) {
                            opt = new BrentOptimizer(Constants.MIN_BRANCH_LENGTH, Constants.MIN_BRANCH_LENGTH, branchConvergenceChecker);
                        } else {
                            opt = new BrentOptimizer(Constants.BRENT_REL, Constants.BRENT_ABS, branchConvergenceChecker);
                        }

                        try {
                            double rnd = Math.random();
                            UnivariatePointValuePair run = opt.optimize(
                                    new MaxEval(Constants.MAX_EVALUATIONS),
                                    new UnivariateObjectiveFunction(new UnivariateFunction() {
                                        @Override
                                        public double value(double branchLength) {
                                            double lnL = getLogLikelihoodForBranch(PhyloUtils.getNodeNumber(t, child), branchLength);

                                            return lnL;
                                        }
                                    }),
                                    GoalType.MAXIMIZE,
                                    // TODO: what if child.getBranchLength() is < Constants.MIN_BRANCH_LENGTH ??
                                    new SearchInterval(Constants.MIN_BRANCH_LENGTH,
                                            Math.max(Constants.MAX_BRANCH_LENGTH + rnd, child.getBranchLength() + rnd),
                                            child.getBranchLength()
                                    )
                            );

                            // TODO: ^^^ A better way to let branch length gradually increase?

                            postOptima = run.getValue();
                            child.setBranchLength(run.getPoint()); // TODO: this only updates the local tree
                            // setBranchLength(child.getNumber(), opt.getResult()); // TODO: this updates the saved partials
                            setBranchLength(PhyloUtils.getNodeNumber(t, child), run.getPoint());
                        } catch (Exception e) {
                            throw new RuntimeException(e);
                        }
                    }
                }
            }

            // Running output
            if (postOptima < preOptima) {
                System.out.printf("↓"); // getting worse!
            } else {
                if ((postOptima - preOptima) > Constants.SMALL_CHANGE ) {
                    System.out.printf("↑"); // getting better
                } else {
                    System.out.printf("→"); // small improvement
                    smallChanges++;
                }
            }

            converged = convergenceChecker.converged(iteration,
                    new PointValuePair(null, preOptima),
                    new PointValuePair(null, postOptima));

            preOptima = postOptima;

            if (smallChanges > Constants.MAX_BRANCH_LENGTH_SMALL_CHANGES) {
                converged = true;
            }
        }

        System.out.println();

        return Pair.of(postOptima, rti.getOriginalRooting());
    }

    public abstract double updateTree(Tree t, Map<Integer, Integer> first, Set<Integer> second);

    public abstract double getLogLikelihoodForBranch(int node, double branchLength);

    public abstract void setBranchLength(int child, double branchLength);

    private interface BLOptCallback {
        public double init(Tree tree);
    }

    public Pair<Double, Tree> optimiseBranchLengths(final Tree tree, final SubstitutionModel model) {
        return optimiseBranchLengths(tree, new BLOptCallback() {
            @Override
            public double init(Tree tree) {
                setRunnerTree(tree);
                return CoreUtils.sum(getLogLikelihood(model, true).values());
            }
        });
    }

    public Pair<Double, Tree> optimiseBranchLengths(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty) {
        // set parameters that don't change (better for distributed runner)
        setRunnerMutation(mutation);
        setRunnerFitnesses(fitnesses);
        setRunnerPenalty(penalty);

        return optimiseBranchLengths(tree, new BLOptCallback() {
            @Override
            public double init(Tree tree) {
                setRunnerTree(tree);
                return CoreUtils.sum(getLogLikelihood(tree, true).values());
            }
        });
    }

    public Pair<Map<Integer, Double>, List<FitnessStore>> optimiseFitness(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, int numberOfOptimRestarts) {
        List<String> clades = Lists.newArrayList("ALL");
        return optimiseFitness(tree, mutation, fitnesses, penalty, clades, numberOfOptimRestarts);
    }

    public Pair<Double, SubstitutionModel> optimiseModel(final Tree tree, final SubstitutionModel model) {
        double[] startingValues = Mapper.getOptimisable(model.getParameters());
        PointValuePair optima;

        setRunnerTree(tree);

        SubstitutionModelFunction function = new SubstitutionModelFunction(model, this);

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);

        //int n = initialGuess.length;
        //MultivariateOptimizer optimiser = new BOBYQAOptimizer(2 * n + 1);
        MultivariateOptimizer optimiser = new NativeLBFGSBMultivariateOptimizer(convergenceChecker);

        //MultivariateOptimizer optimiser = new SimplexOptimizer(convergenceChecker);

        try {
            optima = optimiser.optimize(new ObjectiveFunction(function),
                    GoalType.MAXIMIZE,
                    new InitialGuess(startingValues),//,
                    new NelderMeadSimplex(startingValues.length),
                    new MaxEval(Constants.MAX_EVALUATIONS)
            );
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        System.out.println(" " + optimiser.getEvaluations() + " evaluations");

        if (optima == null)
            throw new RuntimeException("ERROR: Runner.optimiseMutationParameters was unsuccessful (optima == null).");

        double[] d = optima.getPoint();

        Mapper.setOptimisable(model.getParameters(), d);
        model.build();

        //CoreUtils.msg("Runner.optimiseMutationParameters optima: %s -> %s\n", Doubles.join(",", optima.getPoint()), optima.getValue());

        return Pair.of(optima.getValue(), model);
    }

    public Pair<Double, SwMut> optimiseMutationParameters(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty) {
        double[] startingValues = Mapper.getOptimisable(mutation.getParameters());
        PointValuePair optima;

        // To improve performance (esp. the distributed runner), we set the tree
        // and fitnesses for the runner so they do not need to be sent over the
        // wire each iteration
        setRunnerTree(tree);
        setRunnerFitnesses(fitnesses);
        setRunnerPenalty(penalty);

        CoreUtils.msg("Optimising mutation parameters: ");

        SwMutFunction function = new SwMutFunction(mutation, this);

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);

        //int n = initialGuess.length;
        //MultivariateOptimizer optimiser = new BOBYQAOptimizer(2 * n + 1);
        //MultivariateOptimizer optimiser = new NativeLBFGSBMultivariateOptimizer(convergenceChecker);

        MultivariateOptimizer optimiser = new SimplexOptimizer(convergenceChecker);

        try {
            optima = optimiser.optimize(new ObjectiveFunction(function),
                    GoalType.MAXIMIZE,
                    new InitialGuess(startingValues),//,
                    new NelderMeadSimplex(startingValues.length),
                    new MaxEval(Constants.MAX_EVALUATIONS)
            );
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        System.out.println(" " + optimiser.getEvaluations() + " evaluations");

        if (optima == null)
            throw new RuntimeException("ERROR: Runner.optimiseMutationParameters was unsuccessful (optima == null).");

        double[] d = optima.getPoint();

        Mapper.setOptimisable(mutation.getParameters(), d);
        mutation.build();

        //CoreUtils.msg("Runner.optimiseMutationParameters optima: %s -> %s\n", Doubles.join(",", optima.getPoint()), optima.getValue());

        return Pair.of(optima.getValue(), mutation);
    }

    public void setPatterns(Map<Integer, Integer> patternSiteMap, Map<Integer, Integer> patternWeight) {
        if (getSites() == null) throw new RuntimeException("You must setSites() first!");

        this.patternSiteMap = patternSiteMap;
        this.patternWeight = patternWeight;
    }

    public abstract Map<Integer, Double> getLogLikelihood(SubstitutionModel model, boolean save);

    public abstract Map<Integer, Double> getLogLikelihood(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, boolean save);

    public abstract Pair<Map<Integer,Double>, List<FitnessStore>> optimiseFitness(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, final List<String> cladeModel, int numberOfOptimRestarts);

    public abstract void shutdown();
}
