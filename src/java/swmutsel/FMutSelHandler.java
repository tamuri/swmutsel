package swmutsel;

import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleValueChecker;
import pal.tree.Tree;
import swmutsel.cli.ArgumentWriter;
import swmutsel.cli.FMutSel0Arguments;
import swmutsel.model.*;
import swmutsel.model.parameters.*;
import swmutsel.runner.MultiThreadedRunner;
import swmutsel.runner.Runner;
import swmutsel.runner.distributed.DistributedRunner;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;
import swmutsel.utils.PhyloUtils;

import java.util.List;
import java.util.Map;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 23/03/2014 20:45
 */
public class FMutSelHandler implements Handler {
    private final FMutSel0Arguments args;

    public FMutSelHandler(FMutSel0Arguments args) {
        this.args = args;
    }

    public void invoke() {
        final Runner runner = getRunner();
        runner.setRunnerTree(args.tree);
        runner.setSites(args.sites);
        runner.setPatterns(args.patternSiteMap, args.patternWeight);

        final SubstitutionModel model;
        if (this.args.classes == 1) {
            model = getFMutSel0Model();
        } else {
            model = getKFMutSel0Model();
        }

        if (this.args.fix.contains("all")) {
            runner.setRunnerTree(args.tree);
            Map<Integer, Double> siteLogLikelihood = runner.getLogLikelihood(model);

            double sum = CoreUtils.sum(siteLogLikelihood.values());
            CoreUtils.msg("Total log-likelihood = %.7f\n", sum);

            if (args.hessian) {
                int siteCount = 0;
                for (int n : args.patternWeight.values()) siteCount += n;
                HessianCalculator hessian = new HessianCalculator(args.tree, siteCount, runner);
                hessian.get(model);
            }

            runner.shutdown();

            return;

        }

        ConvergenceChecker<PointValuePair> convergenceChecker = new SimpleValueChecker(-1, Constants.CONVERGENCE_TOL);
        PointValuePair lastOptima;
        PointValuePair thisOptima = new PointValuePair(null, Double.NEGATIVE_INFINITY);
        int iteration = 0;
        boolean converged = false;
        double iterationLnL = Double.NEGATIVE_INFINITY;

        Constants.CONVERGENCE_TOL = 5.0;

        CoreUtils.msg("Initial model parameters:\n");
        CoreUtils.msg(model.toString() + "\n");

        long start = System.currentTimeMillis();

        while (!converged) {
            iteration++;

            CoreUtils.msg("Round %d (tol = %.7f)\n", iteration, Constants.CONVERGENCE_TOL);

            // Optimise branch lengths
            CoreUtils.msg("%da - Optimising branch lengths: ", iteration);
            Pair<Double, Tree>  optimaTree = runner.optimiseBranchLengths(args.tree, model);
            args.tree = optimaTree.second;
            iterationLnL = optimaTree.first;
            CoreUtils.msg("Optimised branch lengths %.7f = %.7f\n", PhyloUtils.getTotalTreeLength(args.tree), iterationLnL);

            // Optimise model parameters
            CoreUtils.msg("%db - Optimising model parameters: ", iteration);
            Pair<Double, SubstitutionModel> optima =  runner.optimiseModel(args.tree, model);
            iterationLnL = optima.first;
            CoreUtils.msg("Optimised model parameters %s = %.7f\n", model.toString(), iterationLnL);

            CoreUtils.msg("Elapsed time : %s\n", CoreUtils.getMinSecString(System.currentTimeMillis() - start));

            // Check for convergence
            if (iteration == 1) {
                thisOptima = new PointValuePair(null, iterationLnL);
            } else {
                lastOptima = thisOptima;
                thisOptima = new PointValuePair(null, iterationLnL);
                converged = convergenceChecker.converged(iteration, lastOptima, thisOptima);
            }

            // todo: horrible - fix!
            List<Parameter> params = model.getParameters();
            for (Parameter p : params) {
                if (p instanceof TsTvRatio) args.kappa = ((TsTvRatio) p).get();
                if (p instanceof BaseFrequencies)
                    args.pi = ((BaseFrequencies) p).get();
                if (p instanceof Omega) args.omega = ((Omega) p).get();
                if (p instanceof Probabilities)
                    args.weights = ((Probabilities) p).get();
            }

            if (!converged) {
                String out = args.identifier + "_CHKPNT_" + iteration + ".txt";
                ArgumentWriter.saveStringToFile(args.getCheckpointString() + iterationLnL + "\n", out);
            }

            Constants.CONVERGENCE_TOL /= 2;

            if (Constants.CONVERGENCE_TOL < Constants.MINIMUM_TOL){
                Constants.CONVERGENCE_TOL = Constants.MINIMUM_TOL;
            }
        }

        // We've converged!
        CoreUtils.msg("Finished!\n");
        CoreUtils.msg("Maximum log-likelihood: %.7f\n", iterationLnL);

        // TODO: Save to MLE file
        CoreUtils.msg("MLE:\n%s\n", model.toString());
        CoreUtils.msg("Tree:\n%s\n", args.tree.toString());

        if (args.hessian) {
            int siteCount = 0;
            for (int n : args.patternWeight.values()) siteCount += n;
            HessianCalculator hessian = new HessianCalculator(args.tree, siteCount, runner);
            hessian.get(model);
        }

        runner.shutdown();
    }

    private SubstitutionModel getKFMutSel0Model() {

        final TsTvRatio kappa = new TsTvRatio(args.kappa);
        final Omega omega = new Omega(args.omega);
        final BaseFrequencies pi = new BaseFrequencies(args.pi);

        final Probabilities weights = new Probabilities(4);

        if (args.weights != null) {
            weights.set(args.weights);
        }

        if (args.fix.contains("weights")) {
            weights.setOptimisable(false);
        }

        // if fitnesses weren't given on the command line
        if (args.fitnessList == null) {
            // get some good starting parameters
            Pair<double[], double[][]> startingMixture = SAMMixures.getMixture(args.classes);

            // get rid of those set in ArgumentProcessor.loadFitnessStore()
            args.fitnesses.clear();

            for (int i = 0; i < args.classes; i++) {
                final Fitness f = new Fitness(PhyloUtils.aminoAcidToFitness(startingMixture.second[i], CoreUtils.seq(0, 19)));
                args.fitnesses.put(i + 1, f);
            }

            if (!args.fix.contains("weights")) weights.set(startingMixture.first);
        }

        SubstitutionModel model = new KFMutSel0(kappa, omega, pi, weights, args.fitnesses);
        return model;
    }

    public SubstitutionModel getFMutSel0Model() {
        final Runner runner = getRunner();


        runner.setRunnerTree(args.tree);

        runner.setSites(args.sites);
        runner.setPatterns(args.patternSiteMap, args.patternWeight);

        final TsTvRatio kappa = new TsTvRatio(args.kappa);
        final Omega omega = new Omega(args.omega);
        final BaseFrequencies pi = new BaseFrequencies(args.pi);
        final Fitness f = args.fitnesses.get(1);

        final SubstitutionModel model = new FMutSel0(kappa, omega, pi, f);
        return model;
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

        return runner;
    }

}
