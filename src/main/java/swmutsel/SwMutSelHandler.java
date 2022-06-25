package swmutsel;

import com.google.common.base.Joiner;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleValueChecker;
import pal.tree.Tree;
import swmutsel.cli.ArgumentWriter;
import swmutsel.cli.SwMutSelArguments;
import swmutsel.model.HessianCalculator;
import swmutsel.model.SwMut;
import swmutsel.model.parameters.*;
import swmutsel.runner.MultiThreadedRunner;
import swmutsel.runner.Runner;
import swmutsel.runner.distributed.DistributedRunner;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;
import swmutsel.utils.PhyloUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.List;
import java.util.Map;

/**
 * Parameter optimisation for swMutSel model
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 22/03/2014 04:06
 */
public class SwMutSelHandler implements Handler {
    private final SwMutSelArguments args;

    public SwMutSelHandler(SwMutSelArguments args) {
        this.args = args;
    }

    @Override
    public void invoke() {
        // TODO: Only need to run clademodel if args.fix.contains "mutation" and "branches" and fitnesses specified

        if (args.fix.contains("all")) {

            CoreUtils.msg("Calculating log-likelihood only.\n");
            calculate();

        } else {

            optimise();

        }

    }

    public void optimise() {
        Runner runner = getRunner();

        TsTvRatio kappa = new TsTvRatio(args.kappa);
        BranchScaling c = new BranchScaling(args.scaling);
        BaseFrequencies pi = new BaseFrequencies(args.pi);
        SwMut mutation = new SwMut(kappa, c, pi);

        // TODO: warn about tau!
        if (args.tau != Tau.getDefault()) mutation.setTau(args.tau);

        mutation.build();

        Tree tree = args.tree;
        FitnessStore fitnesses = args.fitnesses;

        Map<Integer, Double> siteLogLikelihood = Maps.newHashMap();

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);
        PointValuePair lastOptima;
        PointValuePair thisOptima = new PointValuePair(null, Double.NEGATIVE_INFINITY);
        int iteration = 0;
        boolean converged = false;
        double iterationLnL = Double.NEGATIVE_INFINITY;

        while (!converged) {
            iteration++;
            long startTime = System.currentTimeMillis();

            // Step 1 - optimise the site-invariant mutational parameters
            if (args.fix.contains("mutation")) {
                CoreUtils.msg("%s (1/3) - Fixed swMutSel mutational parameters.\n", iteration);
            } else {
                Pair<Double, SwMut> optima = runner.optimiseMutationParameters(tree, mutation, fitnesses, args.penalty);
                mutation = optima.second;
                CoreUtils.msg("%s (1/3) - Optimised swMutSel mutational parameters %s = %.7f\n", iteration, mutation.toString(), optima.first);
                iterationLnL = optima.first;
            }

            // Step 2 - optimise the branch lengths
            if (args.fix.contains("branches")) {
                CoreUtils.msg("%s (2/3) - Fixed tree branch lengths.\n", iteration);
            } else {
                Pair<Double, Tree> optima = runner.optimiseBranchLengths(tree, mutation, fitnesses, args.penalty);
                tree = optima.second;
                CoreUtils.msg("%s (2/3) - Optimised tree branch lengths %.7f = %.7f\n", iteration, PhyloUtils.getTotalTreeLength(tree), optima.first);
                iterationLnL = optima.first;
            }

            // Step 3 - optimise the site-wise fitness parameters
            if (args.fix.contains("fitness")) {
                CoreUtils.msg("%s (3/3) - Fixed fitness.\n", iteration);
            } else {
                int restarts = 1;
                if (iteration % args.optRestartInterval == 0) {
                    restarts = args.optRestart;
                    if (restarts != 1 && args.optRestartInterval != 1) {
                        CoreUtils.msg("%s runs of site-wise fitness estimation for iteration %s.\n", restarts, iteration);
                    }
                }

                Pair<Map<Integer, Double>, List<FitnessStore>> optima = runner.optimiseFitness(tree, mutation, fitnesses, args.penalty, restarts);
                siteLogLikelihood = optima.first;
                fitnesses = optima.second.get(0); // NOTE: There's only one set of fitnesses for homogeneous models!
                iterationLnL = CoreUtils.sum(optima.first.values());
                CoreUtils.msg("%s (3/3) - Optimised fitness = %.7f\n", iteration, iterationLnL);
            }


            long endTime = System.currentTimeMillis();
            CoreUtils.msg("%s iteration time = %sms / %.1fm\n", iteration, endTime - startTime, (endTime - startTime) / 60000.0);

            // Check for convergence
            if (args.fix.contains("mutation") && args.fix.contains("branches")) {
                // We only optimised fitness, so we've converged
                converged = true;
            } else if (iteration == 1) {
                // First iteration
                thisOptima = new PointValuePair(null, iterationLnL);
            } else {
                lastOptima = thisOptima;
                thisOptima = new PointValuePair(null, iterationLnL);
                converged = convergenceChecker.converged(iteration, lastOptima, thisOptima);
            }

            // Update args with MLEs
            args.kappa = mutation.getKappa().get();
            args.scaling = mutation.getBranchScaling().get();
            args.pi = mutation.getPi().get();
            args.tree = tree;
            args.fitnesses = fitnesses;

            // Write a checkpoint file if we haven't converged
            if (!converged) {
                String out = args.identifier + "_CHKPNT_" + iteration + ".txt";
                ArgumentWriter.saveStringToFile(args.getCheckpointString(), out);
                try {
                    Files.append("Total log-likelihood = " + String.format("%.7f", iterationLnL), new File(out), Charset.defaultCharset());
                } catch (IOException e) {
                    throw new RuntimeException("Run.optimise ERROR: Could not write checkpoint log-likeihood to " + out);
                }
            }
        }

        // We've converged
        CoreUtils.msg("Finished!\n");
        CoreUtils.msg("Maximum log-likelihood: %.7f\n", iterationLnL);

        String file = args.identifier + "_MLE.txt";
        ArgumentWriter.saveStringToFile(args.getCheckpointString(), file);

        if (siteLogLikelihood.isEmpty()) siteLogLikelihood = runner.getLogLikelihood(args.tree, mutation, args.fitnesses, args.penalty);

        Map<Integer, Pair<Double, Double>> siteDnDs = runner.getDnDs(mutation, args.fitnesses);
        ArgumentWriter.writeSiteInformation(args, siteLogLikelihood, siteDnDs, file);

        // TODO: optimise clade model here, if we need to
        if (args.cladeModel.size() > 1) {
            optimiseCladeModel(runner);
        }

        if (args.hessian) {
            int siteCount = 0;
            for (int n : args.patternWeight.values()) siteCount += n;
            HessianCalculator hessian = new HessianCalculator(args.tree, siteCount, runner);
            hessian.get(args.fitnesses, mutation, args.penalty);
        }

        runner.shutdown();

    }

    private void optimiseCladeModel(Runner runner) {

        // TODO: The tree needs to be rooted somewhere in the first clade - is there anyway around this?
        // Should we do this for the user or expect all internal nodes to be correctly labelled?

        // TODO: Check input - we need homogeneous fitnesses!

        SwMut mutation = new SwMut(args.kappa, args.scaling, args.pi);
        mutation.build();

        Pair<Map<Integer, Double>, List<FitnessStore>> results = runner.optimiseFitness(args.tree, mutation, args.fitnesses, args.penalty, args.cladeModel, args.optRestart);

        String out = args.identifier + "_CLADEMODEL.txt";
        try {
            BufferedWriter writer = Files.newWriter(new File(out), Charset.defaultCharset());
            writer.write("# Results for clade model ( " + Joiner.on(", ").join(args.cladeModel) + " )");
            writer.newLine();

            int clades = args.cladeModel.size();

            for (int site : args.sites.columnKeySet()) {
                writer.write(String.format("%s,%s", site, results.first.get(site)));

                for (int i = 0; i < clades; i++) {
                    writer.write("," + Doubles.join(",", results.second.get(i).get(site).get()));
                }
                writer.newLine();
            }

            writer.newLine();
            writer.close();

        } catch (IOException e) {
            throw new RuntimeException("Run.optimiseCladeModel ERROR: Could not write to " + out);
        }

    }

    public void calculate() {
        Runner runner = getRunner();
        String file = args.identifier + "_CALCULATE.txt";

        SwMut mutation = new SwMut(args.kappa, args.scaling, args.pi);
        if (args.tau != Tau.getDefault()) mutation.setTau(args.tau);
        mutation.build();

        Map<Integer, Double> siteLogLikelihood = runner.getLogLikelihood(args.tree, mutation, args.fitnesses, args.penalty);
        Map<Integer, Pair<Double, Double>> siteDnDs = runner.getDnDs(mutation, args.fitnesses);
        ArgumentWriter.writeSiteInformation(args, siteLogLikelihood, siteDnDs, file);

        double sum = CoreUtils.sum(siteLogLikelihood.values());
        CoreUtils.msg("Total log-likelihood = %.7f\n", sum);

        if (args.hessian) {
            int siteCount = 0;
            for (int n : args.patternWeight.values()) siteCount += n;
            HessianCalculator hessian = new HessianCalculator(args.tree, siteCount, runner);
            hessian.get(args.fitnesses, mutation, args.penalty);
        }

        runner.shutdown();

    }

    private Runner getRunner() {
        Runner runner;

        if (args.distributed) {
            runner = new DistributedRunner(args.hosts);
        } else {
            // TODO: There must be a better place for this!
            ArrayPool.setTreeSize(args.tree.getExternalNodeCount());
            runner = new MultiThreadedRunner(args.threads);
        }

        // todo: could be one action?
        runner.setSites(args.sites);
        runner.setPatterns(args.patternSiteMap, args.patternWeight);
        runner.setRunnerTree(args.tree);

        return runner;
    }


}
