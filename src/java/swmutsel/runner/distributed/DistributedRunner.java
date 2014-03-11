package swmutsel.runner.distributed;

import com.caucho.hessian.client.HessianProxyFactory;
import com.google.common.collect.Iterators;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Table;
import com.google.common.primitives.Ints;
import pal.tree.Tree;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.model.Penalty;
import swmutsel.model.SwMut;
import swmutsel.runner.Runner;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;

import java.net.MalformedURLException;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 30/10/2013 11:34
 */
public class DistributedRunner extends Runner {
    private List<SlaveAPI> slaveService = Lists.newArrayList();
    private List<List<Integer>> slaveSites = Lists.newArrayList();
    private final ExecutorService threadPool;

    public DistributedRunner(List<String> hosts) {
        for (String slave : hosts) {
            try {
                HessianProxyFactory factory = new HessianProxyFactory();
                SlaveAPI service = (SlaveAPI) factory.create(SlaveAPI.class, "http://" + slave + Slave.PATH);
                this.slaveService.add(service);
            } catch (MalformedURLException e) {
                throw new RuntimeException(e);
            }
        }

        this.threadPool = Executors.newFixedThreadPool(this.slaveService.size());
    }


    @Override
    public void setRunnerPenalty(Penalty p) {
        super.setRunnerPenalty(p);
        for (SlaveAPI slave : this.slaveService) slave.setRunnerPenalty(p);
    }

    @Override
    public void setRunnerTree(Tree t) {
        super.setRunnerTree(t);
        for (SlaveAPI slave : this.slaveService) slave.setRunnerTree(t);
    }

    @Override
    public void setRunnerFitnesses(FitnessStore f) {
        super.setRunnerFitnesses(f);
        for (SlaveAPI slave : this.slaveService) slave.setRunnerFitnesses(f);
    }

    @Override
    public void setRunnerMutation(SwMut m) {
        super.setRunnerMutation(m);
        for (SlaveAPI slave : this.slaveService) slave.setRunnerMutation(m);
    }

    @Override
    public void setSites(Table<String, Integer, Byte> sites) {
        super.setSites(sites);

        // Now assign sites to slaves
        if (this.slaveService.isEmpty()) {
            throw new RuntimeException("ERROR: (swmutsel.runner.distributed.DistributedRunner) No slaves found to assign sites!");
        }

        // create an empty site list for each slave
        this.slaveSites.clear();
        for (int i = 0; i < this.slaveService.size(); i++) {
            this.slaveSites.add(Lists.<Integer>newArrayList());
        }

        // add sites to each slave
        Iterator<List<Integer>> slaves = Iterators.cycle(this.slaveSites);
        for (int site : sites.columnKeySet()) {
            slaves.next().add(site);
        }

        // tell each slave the sites it has
        for (int i = 0; i < this.slaveService.size(); i++) {
            this.slaveService.get(i).setSites(Ints.toArray(this.slaveSites.get(i)));
        }
    }

    @Override
    public void setPatterns(Map<Integer, Integer> patternSiteMap, Map<Integer, Integer> patternWeight) {
        super.setPatterns(patternSiteMap, patternWeight);

        // Now assign sites to slaves
        if (this.slaveService.isEmpty()) {
            throw new RuntimeException("ERROR: (swmutsel.runner.distributed.DistributedRunner) No slaves found to assign sites!");
        }

        for (SlaveAPI slave : this.slaveService) {
            slave.setPatterns(patternSiteMap, patternWeight);
        }
    }

    @Override
    public Map<Integer, Double> getLogLikelihood(final Tree tree, final boolean save) {
        // Assumes that setRunnerMutation, setRunnerFitnesses and setRunnerPenalty have already been run!
        List<Future<Map<Integer, Double>>> futures = Lists.newArrayList();

        for (final SlaveAPI slave : this.slaveService) {
            Future<Map<Integer, Double>> future = this.threadPool.submit(new Callable<Map<Integer, Double>>() {
                @Override
                public Map<Integer, Double> call() throws Exception {
                    return slave.getLogLikelihoodForTree(tree, save);
                }
            });

            futures.add(future);
        }

        return collectLogLikelihoodResults(futures);
    }

    @Override
    public Map<Integer, Double> getLogLikelihood(final SwMut mutation, final boolean save) {
        // Assumes that setRunnerTree, setRunnerFitnesses and setRunnerPenalty have already been run!
        List<Future<Map<Integer, Double>>> futures = Lists.newArrayList();

        for (final SlaveAPI slave : this.slaveService) {
            Future<Map<Integer, Double>> future = this.threadPool.submit(new Callable<Map<Integer, Double>>() {
                @Override
                public Map<Integer, Double> call() throws Exception {
                    return slave.getLogLikelihoodForMutation(mutation, save);
                }
            });

            futures.add(future);
        }

        return collectLogLikelihoodResults(futures);
    }

    @Override
    public Map<Integer, Double> getLogLikelihood(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, final boolean save) {
        List<Future<Map<Integer, Double>>> futures = Lists.newArrayList();

        for (final SlaveAPI slave : this.slaveService) {
            Future<Map<Integer, Double>> future = this.threadPool.submit(new Callable<Map<Integer, Double>>() {
                @Override
                public Map<Integer, Double> call() throws Exception {
                    return slave.getLogLikelihood(tree, mutation, fitnesses, penalty, save);
                }
            });

            futures.add(future);
        }

        return collectLogLikelihoodResults(futures);
    }

    private Map<Integer, Double> collectLogLikelihoodResults(List<Future<Map<Integer, Double>>> futures) {
        Map<Integer, Double> allResults = Maps.newHashMap();

        for (Map<Integer, Double> slaveResult : CoreUtils.getFutureResults(futures)) {
            allResults.putAll(slaveResult);
        }

        return allResults;
    }

    @Override
    public double getLogLikelihoodForBranch(final int node, final double branchLength) {
        final List<Future<Double>> futures = Lists.newArrayList();

        for (int i = 0; i < slaveService.size(); i++) {
            final int service_i = i;

            Future<Double> future = threadPool.submit(new Callable<Double>() {
                @Override
                public Double call() throws Exception {
                    final SlaveAPI s = slaveService.get(service_i);
                    return s.getLogLikelihoodForBranch(node, branchLength);
                }
            });

            futures.add(future);
        }

        double total = 0;
        for (double x : CoreUtils.getFutureResults(futures)) total += x;

        return total;
    }

    @Override
    public void setBranchLength(final int child, final double branchLength) {
        List<Future<Void>> futures = Lists.newArrayList();
        for (int j = 0; j < slaveService.size(); j++) {
            final int service_i = j;

            Future<Void> future = threadPool.submit(new Callable<Void>() {
                @Override
                public Void call() throws Exception {
                    final SlaveAPI s = slaveService.get(service_i);
                    s.setBranchLength(child, branchLength);
                    return null;
                }
            });

            futures.add(future);
        }

        CoreUtils.getFutureResults(futures);

    }

    @Override
    public Pair<Map<Integer, Double>, List<FitnessStore>> optimiseFitness(final Tree tree, final SwMut mutation, final FitnessStore fitnesses, final Penalty penalty, final List<String> cladeModel, final int numberOfOptimRestarts) {
        List<Future<Pair<Map<Integer, Double>, List<FitnessStore>>>> futures = Lists.newArrayList();

        // For each slave
        for (int i = 0; i < slaveService.size(); i++) {

            final SlaveAPI slaveService = this.slaveService.get(i);

            // Build a reduced FitnessStore for this particular slave
            final FitnessStore slaveFitnesses = new FitnessStore();
            for (int site : slaveSites.get(i)) slaveFitnesses.put(site, fitnesses.get(site));

            // Submit a request to optimise the fitness
            Future<Pair<Map<Integer, Double>, List<FitnessStore>>> future = threadPool.submit(new Callable<Pair<Map<Integer, Double>, List<FitnessStore>>>() {
                @Override
                public Pair<Map<Integer, Double>, List<FitnessStore>> call() throws Exception {
                    return slaveService.optimiseFitness(tree, mutation, slaveFitnesses, penalty, cladeModel, numberOfOptimRestarts);
                }
            });

            futures.add(future);
        }

        List<Pair<Map<Integer, Double>, List<FitnessStore>>> results = CoreUtils.getFutureResults(futures);

        Map<Integer, Double> logLikelihoods = Maps.newHashMap();
        List<FitnessStore> newFitnesses = Lists.newArrayList();

        for (String ignored : cladeModel) {
            newFitnesses.add(new FitnessStore());
        }

        for (Pair<Map<Integer, Double>, List<FitnessStore>> result : results) {
            logLikelihoods.putAll(result.first);
            for (int i = 0; i < newFitnesses.size(); i++) {
                newFitnesses.get(i).putAll(result.second.get(i));
            }
        }


        return Pair.of(logLikelihoods, newFitnesses);
    }

    @Override
    public void shutdown() {

        for (SlaveAPI slave : this.slaveService) slave.shutdown();

        this.threadPool.shutdown();

    }
}
