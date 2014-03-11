package swmutsel.runner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.random.RandomDataGenerator;
import pal.tree.Tree;
import swmutsel.Constants;
import swmutsel.LikelihoodCalculatorFunction;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.model.*;
import swmutsel.model.parameters.*;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;
import swmutsel.utils.Triple;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.*;

/**
 * REMEMBER: lots of cpus != lots of memory
 * imagine that you're running with 100s of cpus but <1GB memory
 * Think about a small tree with lots of sites - can't store all likelihood calculators
 */
public class MultiThreadedRunner extends Runner {
    private final ExecutorService threadPool;

    private final Map<Integer, LikelihoodCalculator> calculators = Maps.newConcurrentMap();
    // private final List<LikelihoodCalculator> calculators = Lists.newArrayList();

    public MultiThreadedRunner(int threads) {
        this.threadPool = Executors.newFixedThreadPool(threads);
    }

    @Override
    public void shutdown() {
        this.threadPool.shutdown();
    }

    public int totalevals = 0;

/*    class MyFunction implements pal.math.MultivariateFunction {

        private SubstitutionModel model;
        private Tree tree;

        MyFunction(SubstitutionModel model, Tree tree) {
            this.model = model;
            this.tree = tree;
        }

        @Override
            public double evaluate(double[] argument) {
                Mapper.setOptimisable(model.getParameters(), argument);

            model.build();

            *//*    if (!model.parametersValid()) {
                    return Constants.VERY_BAD_LIKELIHOOD;
                } else {

                }
*//*
                double d = CoreUtils.sum(getLogLikelihood(tree, model).values());
                return -d;
            }

            @Override
            public int getNumArguments() {
                int x = 0;
                for (Parameter p : model.getParameters()) {
                    x += p.getOptimisableCount();
                }
                return x;
            }

            @Override
            public double getLowerBound(int n) {

                int pos = 0;
                for (Parameter p : model.getParameters()) {
                    for (int i = 0; i < p.getOptimisableCount(); i++) {
                        if (pos == n) {
                            if (p instanceof BaseFrequencies) {
                                return 0.00001;
                            } else if (p instanceof TsTvRatio) {
                                return 0.00001;
                            } else if (p instanceof Omega) {
                                return 0.00001;
                            } else if (p instanceof Fitness) {
                                return -29;
                            }
                        }
                        pos++;
                    }
                }

                return 0;
            }

            @Override
            public double getUpperBound(int n) {
                int pos = 0;
                for (Parameter p : model.getParameters()) {
                    for (int i = 0; i < p.getOptimisableCount(); i++) {
                        if (pos == n) {
                            if (p instanceof BaseFrequencies) {
                                return 20;
                            } else if (p instanceof TsTvRatio) {
                                return 20;
                            } else if (p instanceof Omega) {
                                return 20;
                            } else if (p instanceof Fitness) {
                                return 29;
                            }
                        }
                        pos++;
                    }
                }
                return 0;
            }

            @Override
            public OrthogonalHints getOrthogonalHints() {
                return null;
            }
        }*/
/*

    public Pair<Map<Integer, Double>, SubstitutionModel> optimiseModel(final Tree tree, final SubstitutionModel model) {



        for (int i = 0; i < 5; i++) {

            List<Parameter> params = model.getParameters();


            for (final Parameter param : params) {

                System.out.printf("%s\n", param.getArgument());

                for (int j = 0; j < param.getOptimisable().length; j++) {

                    final int offset = j;

                    double initial = param.getOptimisable()[j];


                    double lower = 0;
                    double upper = 0;
                    if (param instanceof Fitness) {
                        lower = -25;
                        upper = 25;
                    } else if (param instanceof BaseFrequencies) {
                        System.out.printf("pi: %s\n", Doubles.join(" ", param.getOptimisable()));
                        lower = 0.0000001;
                        upper = 20;
                    } else  {
                        lower = 0.000001;
                        upper = 20;
                    }

                    //System.out.printf("%s %s %s\n", initial, lower, upper);

                    BrentOptimizer optim = new BrentOptimizer();

                    try {
                        double opt = optim.optimize(new UnivariateRealFunction() {
                            @Override
                            public double value(double x) throws FunctionEvaluationException {

                                double[] k = param.getOptimisable();
                                k[offset] = x;
                                param.setOptimisable(k);
                                model.build();

                                double d = CoreUtils.sum(getLogLikelihood(tree, model).values());
                                return d;
                            }
                        }, GoalType.MAXIMIZE, lower + Math.random(), upper + Math.random(), initial);

                        double[] k = param.getOptimisable();
                        k[offset] = opt;
                        param.setOptimisable(k);
                        model.build();
                        //  System.out.printf("%s\n", opt);

                    } catch (MaxIterationsExceededException e) {
                        e.printStackTrace();
                    } catch (FunctionEvaluationException e) {
                        e.printStackTrace();
                    }

                }

            }

            model.build();
            System.out.printf("%s\n", Joiner.on(' ').join(model.getParameters()));
            System.out.printf("-------------------------- %s\n", CoreUtils.sum(getLogLikelihood(tree, model).values()));
            System.out.printf("%s\n", totalevals);
        }




        MyFunction my = new MyFunction(model, tree);


        for (int i = 0; i < my.getNumArguments(); i++) {
            System.out.printf("%s -> %s\n", my.getLowerBound(i), my.getUpperBound(i));
        }


        ConjugateGradientSearch search = new ConjugateGradientSearch();
        search.prin = 2;
        search.conjugateGradientStyle = 1;

        search.optimize(my, Mapper.getOptimisable(model.getParameters()), 1e-5, 1e-5);


        */
/*


        //DirectSearchOptimizer optimiser = new MultiDirectional(10 / (i + 1), 0.1 / (i + 1));
        DirectSearchOptimizer optimiser = new NelderMead();

        //optimiser.setConvergenceChecker(new SimpleScalarValueChecker(-1, 1e-1 / (Math.pow(10,i+ 1))));
        optimiser.setConvergenceChecker(new SimpleScalarValueChecker(-1, Constants.SCALING_THRESHOLD));

        RealPointValuePair optima;

        //final Map<Integer, Double> siteLogLikelihoods = Maps.newTreeMap();

        try {
            optima = optimiser.optimize(new MultivariateRealFunction() {
                private int iteration = 0;

                @Override
                public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {

                    iteration++;

                    //siteLogLikelihoods.clear();

                    Mapper.setOptimisable(model.getParameters(), point);



                    if (!model.parametersValid()) {
                        return Constants.VERY_BAD_LIKELIHOOD;
                    } else {
                        model.build();
                    }


                    double d = CoreUtils.sum(getLogLikelihood(tree, model).values());
                    System.out.printf("%s\n", d);

                    return d;
                }
            }, GoalType.MAXIMIZE, Mapper.getOptimisable(model.getParameters()));
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        if (optima == null) throw new RuntimeException("ERROR: Runner.optimiseMutationParameters was unsuccessful (optima == null).");

        double[] d = optima.getPoint();

        Mapper.setOptimisable(model.getParameters(), d);
        model.build();


*//*


*/
/*

            *//*

*/
/*



                for (final Parameter param : params) {
            if (param.getOptimisable().length == 1) {
                // univariate


                double initial = param.getOptimisable()[0];
                System.out.printf("%s = %s\n", param.getArgument(), initial);

                BrentOptimizer optim = new BrentOptimizer();

                try {
                    double opt = optim.optimize(new UnivariateRealFunction() {
                        @Override
                        public double value(double x) throws FunctionEvaluationException {
                            param.setOptimisable(new double[]{x});
                            model.build();
                            double d = CoreUtils.sum(getLogLikelihood(tree, model).values());
                            return d;
                        }
                    }, GoalType.MAXIMIZE, 0.00001 + Math.random(), 20 + Math.random(), initial);

                    param.setOptimisable(new double[]{opt});
                    model.build();
                    System.out.printf("%s\n", opt);

                } catch (MaxIterationsExceededException e) {
                    e.printStackTrace();
                } catch (FunctionEvaluationException e) {
                    e.printStackTrace();
                }

            } else {
                // multivariate
                double[] initial = param.getOptimisable();

                for (int j = 0; j < initial.length; j++) {
                    boolean addOrMinus = Math.random() > 0.5;

                    if (addOrMinus) {
                        initial[j] += Math.random();
                    } else {
                        initial[j] -= Math.random();
                    }
                }

                System.out.printf("%s = %s\n",param.getArgument(), Doubles.join(" ", initial));

                DirectSearchOptimizer optim = new MultiDirectional();
                // optim.setConvergenceChecker(new SimpleScalarValueChecker(-1, 1e-2 / Math.pow(10, i + 1)));
                optim.setConvergenceChecker(new SimpleScalarValueChecker(-1, 1e-4));


                try {
                    RealPointValuePair opt = optim.optimize(new MultivariateRealFunction() {
                        @Override
                        public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {
                            param.setOptimisable(point);
                            model.build();
                            double d = CoreUtils.sum(getLogLikelihood(tree, model).values());
                            return d;
                        }
                    }, GoalType.MAXIMIZE, initial);


                    param.setOptimisable(opt.getPoint());
                    model.build();

                    System.out.printf("%s -> %s\n", Doubles.join(" ", opt.getPoint()), opt.getValue());


                } catch (FunctionEvaluationException e) {
                    e.printStackTrace();
                } catch (OptimizationException e) {
                    e.printStackTrace();
                }

            }
        }
*//*
*/
/*



        //}







*//*





*/
/*


        double[] startingPoint = Mapper.getOptimisable(model.getParameters());

             int interpolations = ((startingPoint.length + 1) * (startingPoint.length + 2)) / 2;
            //int interpolations = 2 * startingPoint.length + 1;
        //BOBYQAOptimizer bobyqaOptimizer = new BOBYQAOptimizer(interpolations, 2.0, (1e-6 / Math.pow(10, i + 1)));
        BOBYQAOptimizer bobyqaOptimizer = new BOBYQAOptimizer(interpolations);


        PointValuePair result = bobyqaOptimizer.optimize(new MaxEval(20000),
                new SimpleBounds( // pi = 3, kappa = 1, omega = 1, fitness = 19
                        new double[]{0.0001, 0.0001, 0.0001, 1e-2, 1e-2, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29, -29},
                        new double[]{100, 100, 100, 100, 100, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29}
                ),
                org.apache.commons.math3.optim.nonlinear.scalar.GoalType.MAXIMIZE, SimpleBounds.unbounded(startingPoint.length), new InitialGuess(startingPoint), new ObjectiveFunction(new MultivariateFunction() {
            @Override
            public double value(double[] point) {
                Mapper.setOptimisable(model.getParameters(), point);


                if (!model.parametersValid()) {
                    return Constants.VERY_BAD_LIKELIHOOD;
                } else {
                    model.build();
                }

                double d = CoreUtils.sum(getLogLikelihood(tree, model).values());

                System.out.printf("%s\n", d);
                return d;
            }
        }));


        if (result == null)
            throw new RuntimeException("ERROR: Runner.optimiseMutationParameters was unsuccessful (optima == null).");

        double[] dd = result.getPoint();

            System.out.printf("%s point: %s\n", result.getValue(), Doubles.join(" ", dd));
        Mapper.setOptimisable(model.getParameters(), dd);
        model.build();
*//*



*/
/*

            //DirectSearchOptimizer optimiser = new MultiDirectional(10 / (i + 1), 0.1 / (i + 1));
            DirectSearchOptimizer optimiser = new NelderMead();

            //optimiser.setConvergenceChecker(new SimpleScalarValueChecker(-1, 1e-1 / (Math.pow(10,i+ 1))));
            //optimiser.setConvergenceChecker(new SimpleScalarValueChecker(-1, Constants.SCALING_THRESHOLD));

            RealPointValuePair optima;

            //final Map<Integer, Double> siteLogLikelihoods = Maps.newTreeMap();

            try {
                optima = optimiser.optimize(new MultivariateRealFunction() {
                    private int iteration = 0;

                    @Override
                    public double value(double[] point) throws FunctionEvaluationException, IllegalArgumentException {

                        iteration++;

                        //siteLogLikelihoods.clear();

                        Mapper.setOptimisable(model.getParameters(), point);



                        if (!model.parametersValid()) {
                            return Constants.VERY_BAD_LIKELIHOOD;
                        } else {
                            model.build();
                        }


                        double d = CoreUtils.sum(getLogLikelihood(tree, model).values());
                        System.out.printf("%s\n", d);

                    return d;
                }
            }, GoalType.MAXIMIZE, Mapper.getOptimisable(model.getParameters()));
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }

        if (optima == null) throw new RuntimeException("ERROR: Runner.optimiseMutationParameters was unsuccessful (optima == null).");

        double[] d = optima.getPoint();

        Mapper.setOptimisable(model.getParameters(), d);
        model.build();

            //System.out.printf("%s mds point: %s\n", i, Doubles.join(" ", d));

*//*





        // }
        //CoreUtils.msg("Runner.optimiseMutationParameters optima: %s -> %s\n", Doubles.join(",", optima.getPoint()), optima.getValue());
        return Pair.of(getLogLikelihood(tree, model), model);
    }
*/


    public Map<Integer, Double> getLogLikelihood(final Tree tree, final SubstitutionModel model) {
        totalevals++;
        List<Future<Triple<Integer, Double, LikelihoodCalculator>>> futures = Lists.newArrayList();

        for (final Map.Entry<Integer, Map<String, Byte>> site : getSites().columnMap().entrySet()) {
            Future<Triple<Integer, Double, LikelihoodCalculator>> future = threadPool.submit(new Callable<Triple<Integer, Double, LikelihoodCalculator>>() {

                @Override
                public Triple<Integer, Double, LikelihoodCalculator> call() throws Exception {
                    // SubstitutionModel model = new FMutSel0(new TsTvRatio(7.96527), new Omega(0.05327), new BaseFrequencies(new double[]{0.17572,0.28418,0.45954,0.08056}));

                    LikelihoodCalculator calculator;

                    calculator = new LikelihoodCalculator(tree, site.getValue(), model);


                    calculator.getStorage();
                    double result = calculator.getLogLikelihood();
                    calculator.releaseStorage();

                    return Triple.of(site.getKey(), result, null);
                }
            });

            futures.add(future);
        }

        Map<Integer, Double> siteLogLikelihood = Maps.newHashMap();
        calculators.clear();

        for (Triple<Integer, Double, LikelihoodCalculator> p : CoreUtils.getFutureResults(futures)) {
            siteLogLikelihood.put(p.first, p.second);
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

                    calculator.getStorage();
                    double result = calculator.getLogLikelihood();
                    calculator.releaseStorage();

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
            if (saveCalculators) calculators.put(p.first, p.third);
        }

        // For each pattern -> site mapping
        for (Map.Entry<Integer, Integer> e : getPatternSiteMap().entrySet()) {
            // If the pattern is not site AND we have the relevant site (for DistributedRunner)
            if (e.getKey() != e.getValue() && siteLogLikelihood.containsKey(e.getValue())) {
                siteLogLikelihood.put(e.getKey(), siteLogLikelihood.get(e.getValue()));
            }
        }

        return siteLogLikelihood;
    }

    @Override
    public double getLogLikelihoodForBranch(final int node, final double branchLength) {
        final List<Future<Pair<Integer, Double>>> futures = Lists.newArrayList();
        for (final Map.Entry<Integer, LikelihoodCalculator> e : calculators.entrySet()) {
            Future<Pair<Integer, Double>> future = threadPool.submit(new Callable<Pair<Integer, Double>>() {
                @Override
                public Pair<Integer, Double> call() throws Exception {
                    return Pair.of(e.getKey(), e.getValue().getLogLikelihoodForBranch(node, branchLength));
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
    public void setBranchLength(int child, double branchLength) {
        for (LikelihoodCalculator l : calculators.values()) {
            l.setBranch(child, branchLength);
        }
    }

    @Override
    public Pair<Map<Integer, Double>, List<FitnessStore>> optimiseFitness(
            final Tree tree,
            final SwMut mutation,
            final FitnessStore fitnesses,
            final Penalty penalty,
            final List<String> cladeModel,
            final int numberOfOptimRestarts) {

        final List<Future<Triple<Integer, Double, ArrayList<Fitness>>>> futures = Lists.newArrayList();

        final RandomDataGenerator randomData = new RandomDataGenerator();

        CompletionService<Triple<Integer, Double, ArrayList<Fitness>>> completionService =
                new ExecutorCompletionService<Triple<Integer, Double, ArrayList<Fitness>>>(threadPool);

        for (final Map.Entry<Integer, Map<String, Byte>> site : getSites().columnMap().entrySet() ) {
            Future<Triple<Integer, Double, ArrayList<Fitness>>> future = completionService.submit(new Callable<Triple<Integer, Double, ArrayList<Fitness>>>() {
                @Override
                public Triple<Integer, Double, ArrayList<Fitness>> call() throws Exception {

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

                        LikelihoodCalculatorFunction function = new LikelihoodCalculatorFunction(calculator, Lists.newArrayList(models.values()));

                        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);
                        SimplexOptimizer optimiser = new SimplexOptimizer(convergenceChecker);

                        double[] initialGuess = function.getCurrentParameters();

                        PointValuePair optima;
                        optima = optimiser.optimize(
                                new ObjectiveFunction(function),
                                GoalType.MAXIMIZE,
                                new MaxEval(Constants.MAX_EVALUATIONS),
                                new InitialGuess(initialGuess),
                                new NelderMeadSimplex(initialGuess.length)
                        );

                        if (optimiser.getEvaluations() >= optimiser.getMaxEvaluations()) {
                            // Did not converge - use the last attempted point and assume it's close to the optima...!
                            optima = new PointValuePair(optima.getPoint(), function.value(optima.getPoint()));
                            CoreUtils.msg("MultiThreadedRunner.optimiseFitness - Site %s exceeded maximum number of evaluations (%s).%n", site, optimiser.getEvaluations());
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
                }
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
                    newFitnesses.get(i).set(t.first, t.third.get(i));
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
            if (e.getKey() != e.getValue() && logLikelihoods.containsKey(e.getValue())) {
                logLikelihoods.put(e.getKey(), logLikelihoods.get(e.getValue()));

                for (int i = 0; i < newFitnesses.size(); i++) {
                    newFitnesses.get(i).set(e.getKey(), newFitnesses.get(i).get(e.getValue()));
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

}
