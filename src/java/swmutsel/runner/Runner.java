package swmutsel.runner;


import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.google.common.collect.Table;
import org.apache.commons.math3.analysis.MultivariateFunction;
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
import swmutsel.Constants;
import swmutsel.model.SwMutSel;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.model.Penalty;
import swmutsel.model.SwMut;
import swmutsel.model.parameters.Mapper;
import swmutsel.trees.RerootedTreeIterator;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;
import swmutsel.utils.PhyloUtils;

import java.util.List;
import java.util.Map;
import java.util.Set;

public abstract class Runner {

    private Table<String, Integer, Byte> sites;
    private Map<Integer, Integer> patternSiteMap;
    private Map<Integer, Integer> patternWeight;

    private Tree runnerTree;
    private FitnessStore runnerFitnesses;
    private SwMut runnerMutation;
    private Penalty runnerPenalty;

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
        this.runnerTree = runnerTree;
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

    public void setSites(Table<String, Integer, Byte> sites) {
        this.sites = sites;
    }

    public Table<String, Integer, Byte> getSites() {
        return sites;
    }

    public void setPatterns(Map<Integer, Integer> patternSiteMap, Map<Integer, Integer> patternWeight) {

        if (getSites() == null) throw new RuntimeException("You must setSites() first!");

        this.patternSiteMap = patternSiteMap;
        this.patternWeight = patternWeight;
    }

    protected Map<Integer,Integer> getPatternSiteMap() {
        return patternSiteMap;
    }

    protected Map<Integer, Integer> getPatternWeight() {
        return patternWeight;
    }

    public Pair<Double, SwMut> optimiseMutationParameters(
            final Tree tree,
            final SwMut mutation,
            final FitnessStore fitnesses,
            final Penalty penalty) {

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);

        double[] startingValues = Mapper.getOptimisable(mutation.getParameters());

        MultivariateOptimizer optimiser = new SimplexOptimizer(convergenceChecker);
        //int n = initialGuess.length;
        //MultivariateOptimizer optimiser = new BOBYQAOptimizer(2 * n + 1);

        PointValuePair optima;

        // To improve performance (esp. the distributed runner), we set the tree
        // and fitnesses for the runner so they do not need to be sent over the
        // wire each iteration
        setRunnerTree(tree);
        setRunnerFitnesses(fitnesses);
        setRunnerPenalty(penalty);

        CoreUtils.msg("Optimising mutation parameters: ");

        try {
            optima = optimiser.optimize(new ObjectiveFunction(new MultivariateFunction() {
                private int iteration = 0;
                private double current = Double.NEGATIVE_INFINITY;

                @Override
                public double value(double[] point) {

                    Mapper.setOptimisable(mutation.getParameters(), point);

                    if (!mutation.parametersValid()) {
                        return Constants.VERY_BAD_LIKELIHOOD;
                    } else {
                        mutation.build();
                    }

                    double lnL = CoreUtils.sum(getLogLikelihood(mutation, false).values());

                    iteration++;

                    // Running output
                    if (iteration == 1) System.out.printf("%.7f ", lnL);
                    if (iteration % 25 == 0) {
                        if (lnL < current) {
                            System.out.printf("↓"); // getting worse!
                        } else {
                            if ((lnL - current) > Constants.SMALL_CHANGE ) {
                                System.out.printf("↑"); // getting better
                            } else {
                                System.out.printf("→"); // small improvement
                            }
                        }
                        current = lnL;
                    }

                    return lnL;
                }
            }),
                    GoalType.MAXIMIZE,
                    new InitialGuess(startingValues),
                    new NelderMeadSimplex(startingValues.length),
                    //new SimpleBounds(new double[]{0.0001, -2, -2, -2, 0.0001}, new double[]{50  , 2, 2, 2, 50}),
                    new MaxEval(Constants.MAX_EVALUATIONS)
            );
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        System.out.println();

        if (optima == null)
            throw new RuntimeException("ERROR: Runner.optimiseMutationParameters was unsuccessful (optima == null).");

        double[] d = optima.getPoint();

        Mapper.setOptimisable(mutation.getParameters(), d);
        mutation.build();

        //CoreUtils.msg("Runner.optimiseMutationParameters optima: %s -> %s\n", Doubles.join(",", optima.getPoint()), optima.getValue());

        return Pair.of(optima.getValue(), mutation);
    }


    public Pair<Double, Tree> optimiseBranchLengths(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty) {

        double postOptima = Double.NEGATIVE_INFINITY, preOptima = Double.NEGATIVE_INFINITY;

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);
        ConvergenceChecker<UnivariatePointValuePair> branchConvergenceChecker = new SimpleUnivariateValueChecker(-1, Constants.MIN_BRANCH_LENGTH);

        int iteration = 0;

        RerootedTreeIterator rti = new RerootedTreeIterator(tree);

        boolean converged = false;

        // set parameters that don't change (better for distributed runner)
        setRunnerMutation(mutation);
        setRunnerFitnesses(fitnesses);
        setRunnerPenalty(penalty);

        CoreUtils.msg("Optimising branch lengths: ");

        while (!converged) {
            // Save branches we've already optimised so we don't do them again
            Set<Set<Node>> optimisedBranches = Sets.newHashSet();

            iteration++;

            boolean start = true;

            // loop over every rerooted tree
            for (final Tree t : rti) {

                // We get the current log-likelihood and also save the partial likelihoods around the root node
                // ready for getLogLikelihoodForBranch call (see below).
                // double lnL = CoreUtils.sum(getLogLikelihood(t, mutation, fitnesses, penalty, true).values());
                double lnL = CoreUtils.sum(getLogLikelihood(t, true).values());

                if (Double.isInfinite(postOptima) && Double.isInfinite(preOptima))
                    System.out.printf("%.7f ", lnL);

                // if this is the first iteration over branches, use this log-likelihood to check for convergence
                if (start) {
                    preOptima = lnL;
                    start = false;
                }

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

                        UnivariateOptimizer opt = new BrentOptimizer(Constants.BRENT_REL, Constants.BRENT_ABS, branchConvergenceChecker);

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
                                    new SearchInterval(Constants.MIN_BRANCH_LENGTH,
                                            Math.max(Constants.MAX_BRANCH_LENGTH + rnd, child.getBranchLength() + rnd),
                                            child.getBranchLength()
                                    ));

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
                }
            }

            converged = convergenceChecker.converged(iteration,
                    new PointValuePair(null, preOptima),
                    new PointValuePair(null, postOptima));

            preOptima = postOptima;


        }

        System.out.println();

        return Pair.of(postOptima, rti.getOriginalRooting());
    }

    public Map<Integer, Double> getLogLikelihood(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty) {
        return getLogLikelihood(tree, mutation, fitnesses, penalty, false);
    }

    public abstract Map<Integer, Double> getLogLikelihood(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, boolean save);

    public Map<Integer, Double> getLogLikelihood(final SwMut mutation, boolean save) {
        return getLogLikelihood(getRunnerTree(), mutation, getRunnerFitnesses(), getRunnerPenalty(), save);
    }

    public Map<Integer, Double> getLogLikelihood(final Tree tree, boolean save) {
        return getLogLikelihood(tree, getRunnerMutation(), getRunnerFitnesses(), getRunnerPenalty(), save);
    }

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

    public abstract double getLogLikelihoodForBranch(int node, double branchLength);

    public abstract void setBranchLength(int child, double branchLength);

    public abstract Pair<Map<Integer,Double>, List<FitnessStore>> optimiseFitness(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, final List<String> cladeModel, int numberOfOptimRestarts);

    public Pair<Map<Integer, Double>, List<FitnessStore>> optimiseFitness(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, int numberOfOptimRestarts) {
        List<String> clades = Lists.newArrayList("ALL");
        return optimiseFitness(tree, mutation, fitnesses, penalty, clades, numberOfOptimRestarts);
    }

    public abstract void shutdown();

}
