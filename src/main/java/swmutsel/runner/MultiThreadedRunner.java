package swmutsel.runner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.random.RandomDataGenerator;
import pal.tree.Tree;
import swmutsel.Constants;
import swmutsel.model.*;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.MathUtils;
import swmutsel.utils.Pair;
import swmutsel.utils.Triple;

import java.util.*;
import java.util.concurrent.*;

/**
 * REMEMBER: lots of cpus != lots of memory
 * imagine that you're running with 100s of cpus but <1GB memory
 * Think about a small tree with lots of sites - can't store all likelihood calculators
 */
public class MultiThreadedRunner extends Runner {
// ------------------------------ FIELDS ------------------------------

    public int totalevals = 0;
    private final ExecutorService threadPool;

    private final Map<Integer, LikelihoodCalculator> calculators = Maps.newConcurrentMap();

// --------------------------- CONSTRUCTORS ---------------------------

    // private final List<LikelihoodCalculator> calculators = Lists.newArrayList();
    public MultiThreadedRunner(int threads) {
        this.threadPool = Executors.newFixedThreadPool(threads);
    }

// -------------------------- OTHER METHODS --------------------------

    final Map<Integer, List<LikelihoodCalculator>> mixtureCalculators = Maps.newConcurrentMap();


    @Override
    public Map<Integer, Double> getLogLikelihood(final SubstitutionModel model, final boolean save) {
        totalevals++;
        List<Future<Triple<Integer, Double, LikelihoodCalculator>>> futures = Lists.newArrayList();
        mixtureCalculators.clear();
        calculators.clear();

        for (final Map.Entry<Integer, Map<String, Byte>> site : getSites().columnMap().entrySet()) {
            Future<Triple<Integer, Double, LikelihoodCalculator>> future = threadPool.submit(new Callable<Triple<Integer, Double, LikelihoodCalculator>>() {
                @Override
                public Triple<Integer, Double, LikelihoodCalculator> call() throws Exception {

                    double result;

                    if (model instanceof MixtureModel) {

                        List<LikelihoodCalculator> saved = Lists.newArrayList();

                        MixtureModel mixtureModel = (MixtureModel) model;
                        SubstitutionModel[] models = mixtureModel.getModels();
                        double[] each = new double[models.length];

                        for (int i = 0; i < models.length; i++) {
                            LikelihoodCalculator calculator = new LikelihoodCalculator(getRunnerTree(), site.getValue(), models[i]);

                            if (save) saved.add(calculator);

                            if (!save) calculator.getStorage();
                            each[i] = calculator.getLikelihood(mixtureModel.getWeight(i));
                            if (!save) calculator.releaseStorage();
                        }

                        if (save) {
                            mixtureCalculators.put(site.getKey(), saved);
                        }

                        result = MathUtils.sumLogs(each);

                        // TODO: we actually don't need this likelihoodcalculator - fix the way getLikelihoodForBranch() works!!
                        return Triple.of(site.getKey(), result, null);


                    } else {
                        LikelihoodCalculator calculator = new LikelihoodCalculator(getRunnerTree(), site.getValue(), model);

                        if (!save) calculator.getStorage();
                        result = calculator.getLogLikelihood();
                        if (!save) calculator.releaseStorage();

                        if (save) {
                            return Triple.of(site.getKey(), result, calculator);
                        } else {
                            return Triple.of(site.getKey(), result, null);
                        }
                    }

                }
            });

            futures.add(future);
        }

        Map<Integer, Double> siteLogLikelihood = Maps.newHashMap();


        for (Triple<Integer, Double, LikelihoodCalculator> p : CoreUtils.getFutureResults(futures)) {
            siteLogLikelihood.put(p.first, p.second);
            if (!(model instanceof MixtureModel) && save) calculators.put(p.first, p.third);
        }

        // For each pattern -> site mapping
        for (Map.Entry<Integer, Integer> e : getPatternSiteMap().entrySet()) {
            int pattern = e.getKey();
            int site = e.getValue();
            // If the pattern is not site AND we have the relevant site (for DistributedRunner)
            if (pattern != site && siteLogLikelihood.containsKey(e.getValue())) {
                siteLogLikelihood.put(e.getKey(), siteLogLikelihood.get(e.getValue()));
            }
        }

        return siteLogLikelihood;
    }

    @Override
    public Map<Integer, Double> getLogLikelihood(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, final boolean saveCalculators) {
        List<Future<Triple<Integer, Double, LikelihoodCalculator>>> futures = Lists.newArrayList();

        for (final Map.Entry<Integer, Map<String, Byte>> site : getSites().columnMap().entrySet()) {
            Future<Triple<Integer, Double, LikelihoodCalculator>> future = threadPool.submit(new Callable<Triple<Integer, Double, LikelihoodCalculator>>() {
                @Override
                public Triple<Integer, Double, LikelihoodCalculator> call() throws Exception {
                    Fitness f = fitnesses.get(site.getKey());
                    SubstitutionModel model = new SwMutSel(mutation, f);

                    LikelihoodCalculator calculator;

                    if (penalty == null) {
                        calculator = new LikelihoodCalculator(tree, site.getValue(), model);
                    } else {
                        calculator = new PenalisedLikelihoodCalculator(tree, site.getValue(), model, penalty, f);
                    }

                    if (!saveCalculators) calculator.getStorage();
                    double result = calculator.getLogLikelihood();
                    if (!saveCalculators) calculator.releaseStorage();

                    if (saveCalculators) {
                        return Triple.of(site.getKey(), result, calculator);
                    } else {
                        return Triple.of(site.getKey(), result, null);
                    }
                }
            });

            futures.add(future);
        }

        Map<Integer, Double> siteLogLikelihood = Maps.newHashMap();
        calculators.clear();

        for (Triple<Integer, Double, LikelihoodCalculator> p : CoreUtils.getFutureResults(futures)) {
            siteLogLikelihood.put(p.first, p.second);
            if (saveCalculators) {
                calculators.put(p.first, p.third);
            }
        }

        // For each pattern -> site mapping
        for (Map.Entry<Integer, Integer> e : getPatternSiteMap().entrySet()) {
            int pattern = e.getKey();
            int site = e.getValue();
            // If the pattern is not site AND we have the relevant site (for DistributedRunner)
            if (pattern != site && siteLogLikelihood.containsKey(e.getValue())) {
                siteLogLikelihood.put(e.getKey(), siteLogLikelihood.get(e.getValue()));
            }
        }

        return siteLogLikelihood;
    }

    @Override
    public double updateTree(final Tree t, final Map<Integer, Integer> newNodeMapping, final Set<Integer> changedNodes) {
        final List<Future<Pair<Integer, Double>>> futures = Lists.newArrayList();

        for (final int site : getSites().columnKeySet()) {

            Future<Pair<Integer, Double>> future = threadPool.submit(new Callable<Pair<Integer, Double>>() {
                @Override
                public Pair<Integer, Double> call() throws Exception {

                    if (!mixtureCalculators.isEmpty()) {

                        int mixtures = mixtureCalculators.entrySet().iterator().next().getValue().size();
                        List<LikelihoodCalculator> mixture = mixtureCalculators.get(site);

                        double[] each = new double[mixtures];
                        int pos = 0;
                        for (LikelihoodCalculator l : mixture) {
                            each[pos++] = l.updateTree(t, newNodeMapping, changedNodes);
                        }

                        double lnL = MathUtils.sumLogs(each);

                        return Pair.of(site, lnL);

                    } else {
                        return Pair.of(site, calculators.get(site).updateTree(t, newNodeMapping, changedNodes));
                    }


                }
            });
            futures.add(future);

        }

        double sum = 0;
        for (Pair<Integer, Double> siteLogLikelihood : CoreUtils.getFutureResults(futures)) {
            sum += siteLogLikelihood.second * getPatternWeight().get(siteLogLikelihood.first);
        }

        return sum;
    }

    @Override
    public double getLogLikelihoodForBranch(final int node, final double branchLength) {
        final List<Future<Pair<Integer, Double>>> futures = Lists.newArrayList();

        for (final int site : getSites().columnKeySet()) {
        /*for (final Map.Entry<Integer, LikelihoodCalculator> e : calculators.entrySet()) {*/
            Future<Pair<Integer, Double>> future = threadPool.submit(new Callable<Pair<Integer, Double>>() {
                @Override
                public Pair<Integer, Double> call() throws Exception {

                    if (!mixtureCalculators.isEmpty()) {

                        int mixtures = mixtureCalculators.entrySet().iterator().next().getValue().size();
                        List<LikelihoodCalculator> mixture = mixtureCalculators.get(site);

                        double[] each = new double[mixtures];
                        int pos = 0;
                        for (LikelihoodCalculator l : mixture) {
                            each[pos++] = l.getLogLikelihoodForBranch(node, branchLength);
                        }

                        double lnL = MathUtils.sumLogs(each);




                        return Pair.of(site, lnL);

                    } else {
                        return Pair.of(site, calculators.get(site).getLogLikelihoodForBranch(node, branchLength));
                    }


                }
            });
            futures.add(future);
        }

        double sum = 0;
        for (Pair<Integer, Double> siteLogLikelihood : CoreUtils.getFutureResults(futures)) {
            sum += siteLogLikelihood.second * getPatternWeight().get(siteLogLikelihood.first);
        }

        return sum;
    }

    @Override
    public Pair<Map<Integer, Double>, List<FitnessStore>> optimiseFitness(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, final List<String> cladeModel, final int numberOfOptimRestarts) {
        final List<Future<Triple<Integer, Double, ArrayList<Fitness>>>> futures = Lists.newArrayList();

        final RandomDataGenerator randomData = new RandomDataGenerator();

        CompletionService<Triple<Integer, Double, ArrayList<Fitness>>> completionService =
                new ExecutorCompletionService<Triple<Integer, Double, ArrayList<Fitness>>>(threadPool);

        for (final Map.Entry<Integer, Map<String, Byte>> site : getSites().columnMap().entrySet() ) {
            Future<Triple<Integer, Double, ArrayList<Fitness>>> future = completionService.submit(() -> {
                // We allow for multiple attempts to optimise fitness
                List<Pair<Double, ArrayList<Fitness>>> optimals = Lists.newArrayList();

                for (int run = 0; run < numberOfOptimRestarts; run++) {
                    ArrayList<Fitness> siteFitness = Lists.newArrayList();
                    LinkedHashMap<String, SubstitutionModel> models = Maps.newLinkedHashMap();

                    for (String clade : cladeModel) {
                        Fitness f = fitnesses.get(site.getKey()).copy();

                        // If we're anything but the very first optimisation attempt
                        if (run > 0) {
                            // starting fitnesses are drawn from a uniform distribution
                            double[] randf = new double[19];
                            for (int j = 0; j < randf.length; j++)
                                randf[j] = randomData.nextUniform(-Constants.RANDOM_INITIAL_FITNESS_RANGE, Constants.RANDOM_INITIAL_FITNESS_RANGE);
                            f.setOptimisable(randf);
                        }

                        siteFitness.add(f);
                        models.put(clade, new SwMutSel(mutation, f));
                    }

                    Map<String, Byte> states = site.getValue();

                    LikelihoodCalculator calculator;

                    if (penalty == null) {
                        calculator = new LikelihoodCalculator(tree, states, models);
                    } else {
                        if (cladeModel.size() == 1) {
                            calculator = new PenalisedLikelihoodCalculator(tree, site.getValue(), models, penalty, siteFitness.get(0));
                        } else {
                            throw new RuntimeException("MultiThreadedRunner.optimiseFitness - Penalised likelihood not implemented for non-homogeneous models");
                        }
                    }

                    calculator.getStorage();

                    SwMutSelFunction function = new SwMutSelFunction(calculator, Lists.newArrayList(models.values()));
                    ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.VALUE_CONVERGENCE_TOL);
                    MultivariateOptimizer optimiser = new SimplexOptimizer(convergenceChecker);

                    //MultivariateOptimizer optimiser = new LBFGSMultivariateOptimizer(convergenceChecker);
                    //MultivariateOptimizer optimiser = new NativeLBFGSBMultivariateOptimizer(convergenceChecker);

                    double[] initialGuess = function.getCurrentParameters();

                    PointValuePair optima;
                    optima = optimiser.optimize(
                            new ObjectiveFunction(function),
                            GoalType.MAXIMIZE,
                            new MaxEval(Constants.MAX_EVALUATIONS),
                            new MaxIter(Constants.MAX_ITERATIONS),
                            new InitialGuess(initialGuess),
                            new NelderMeadSimplex(initialGuess.length)
                    );

                    if (optimiser.getIterations() >= Constants.MAX_ITERATIONS) {
                        System.out.printf("[%d, %d, %d] ", site.getKey(), optimiser.getEvaluations(), optimiser.getIterations());
                    }

                    if (optima == null) return null; // ERROR!!

                    function.value(optima.getPoint());

                    calculator.releaseStorage();

                    optimals.add(Pair.of(optima.getValue(), siteFitness));
                }

                int best = 0;
                for (int i = 1; i < optimals.size(); i++)
                    if (optimals.get(i).first > optimals.get(best).first) best = i;

                return Triple.of(site.getKey(), optimals.get(best).first, optimals.get(best).second);
            });
            futures.add(future);
        }


        Map<Integer, Double> logLikelihoods = Maps.newHashMap();


        List<FitnessStore> newFitnesses = Lists.newArrayList();

        // TODO: not the same clademodel we set up in callable!!!!

        // if we're running homogenous model
        if (cladeModel.size() == 1) {
            // only one set of fitnesses
            newFitnesses.add(new FitnessStore());
        } else {
            // otherwise, we're running non-homogeneous model
            for (int i = 0; i < cladeModel.size(); i++) {
                // add as many fitness stores as there are clades
                newFitnesses.add(new FitnessStore());
            }
        }

        int totalSites = getSites().columnKeySet().size();
        int[] percentageSites = getPercentageSites(totalSites);
        CoreUtils.msg("Optimising fitness for %s sites: ", totalSites);

        for (int n = 1; n <= totalSites; n++) {
            try {
                Triple<Integer, Double, ArrayList<Fitness>> t = completionService.take().get();
                logLikelihoods.put(t.first, t.second);
                for (int i = 0; i < t.third.size(); i++) {
                    newFitnesses.get(i).put(t.first, t.third.get(i));
                }
            } catch (Exception e) {
                throw new RuntimeException(e);
            }

            int progress = Ints.indexOf(percentageSites, n);
            if (progress >= 0) {
                System.out.printf("%s%% ", (progress + 1) * 10);
            }
        }

        System.out.println();

        // TODO: check!!!
        for (Map.Entry<Integer, Integer> e : getPatternSiteMap().entrySet()) {
            int pattern = e.getKey();
            int site = e.getValue();
            if (pattern != site && logLikelihoods.containsKey(e.getValue())) {
                logLikelihoods.put(e.getKey(), logLikelihoods.get(e.getValue()));

                for (FitnessStore fs : newFitnesses) {
                    fs.put(e.getKey(), fs.get(e.getValue()));
                }
            }
        }


        //CoreUtils.msg("MultiThreadedRunner.optimiseFitness -> %s\n", total.get());

        // return CoreUtils.sum(Doubles.toArray(logLikelihoods.values()));
        return Pair.of(logLikelihoods, newFitnesses);

        // -1941.6130952000872
        // -1941.6130952000954
    }

    private int[] getPercentageSites(int totalSites) {
        int[] percentageSites = new int[10];
        for (int i = 0; i < percentageSites.length; i++) {
            double perc = 0.1 * (i + 1);
            percentageSites[i] = (int) Math.floor(perc * (double) (totalSites));
        }
        return percentageSites;
    }

    @Override
    public void setBranchLength(int child, double branchLength) {
        if (mixtureCalculators.size() > 0) {
            for (List<LikelihoodCalculator> ls : mixtureCalculators.values()) {
                for (LikelihoodCalculator l : ls) {
                    l.setBranch(child, branchLength);
                }
            }
        } else {
            for (LikelihoodCalculator l : calculators.values()) {
                l.setBranch(child, branchLength);
            }
        }
    }

    @Override
    public void shutdown() {
        this.threadPool.shutdown();
    }
}
