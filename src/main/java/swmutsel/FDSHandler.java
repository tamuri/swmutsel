package swmutsel;

import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import org.apache.commons.math3.optim.*;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.util.FastMath;
import swmutsel.cli.FDSArguments;
import swmutsel.model.*;
import swmutsel.model.parameters.*;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;
import swmutsel.utils.Pair;
import swmutsel.utils.PhyloUtils;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;

public class FDSHandler implements Handler {
    private final FDSArguments args;
    public FDSHandler(FDSArguments args) {
        this.args = args;
    }

    @Override
    public void invoke() {
        long start = System.currentTimeMillis();

        String out = "";

        TsTvRatio kappa = new TsTvRatio(args.kappa);
        BranchScaling c = new BranchScaling(args.scaling);
        BaseFrequencies pi = new BaseFrequencies(args.pi);
        SwMut mutation = new SwMut(kappa, c, pi);

        // Fix the fitness of the most common residue to 0
        int[] aminoAcidOrder = Ints.toArray(PhyloUtils.getOrderedAminoAcids(args.sites.column(args.site).values()));
        CoreUtils.msg("Fitness of residue %d (%s) fixed to 0\n", aminoAcidOrder[0] + 1, GeneticCode.getInstance().getAminoAcidCharByIndex(aminoAcidOrder[0]));

        // If user supplied fitnesses, normalise so that the fitness of the most common residue is 0 (if necessary)
        if (args.fitness.size() == GeneticCode.AMINO_ACID_STATES && args.fitness.get(aminoAcidOrder[0]) != 0) {
            double offset = args.fitness.get(aminoAcidOrder[0]);
            for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                args.fitness.set(i, args.fitness.get(i) - offset);
            }
            CoreUtils.msg("Normalised fitnesses so amino acid %s = 0\n", GeneticCode.getInstance().getAminoAcidCharByIndex(aminoAcidOrder[0]));
        }

        if (args.Zpenalty == 0) {
            CoreUtils.msg("No penalty on FDS parameter Z\n");
        } else if (args.Zpenalty > 0) {
            CoreUtils.msg("Using penalty exponential,%2.5f on FDS parameter Z\n", args.Zpenalty);
        } else {
            CoreUtils.msg("ERROR: Lambda parameter (-Zpenalty) must be positive!");
            System.exit(1);
        }

        // Output the synonymous and non-synonymous proportions for pure mutation model and no-selection model (all F=0)
        {
            Fitness fitness = new Fitness(aminoAcidOrder);
            SwMutSel mutSel = new SwMutSel(mutation, fitness);
            LikelihoodCalculator calculator = getCalculator(mutSel, fitness);
            SwMutSelFunction function = new SwMutSelFunction(calculator, mutSel);
            double lnl = function.value(function.getCurrentParameters());
            Pair<Double, Double> nonSynAndSyn = mutSel.getRhoNonSynAndSyn();
            Pair<Double, Double> neutralDnDs = mutation.getRhoNonSynAndSyn();
            String result = String.format("%d %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f ", args.site, neutralDnDs.first, neutralDnDs.second, lnl, nonSynAndSyn.first, nonSynAndSyn.second, mutSel.getExpectedSubsPerSite());
            System.out.print(result);
            out = out + result;
        }

        // Optimise site under swMutSel+F (frequency dependent selection model)
        {
            Fitness fitness;
            if (args.fitness.size() == GeneticCode.AMINO_ACID_STATES) {
                fitness = new Fitness(aminoAcidOrder, Doubles.toArray(args.fitness));
            } else {
                fitness = new Fitness(aminoAcidOrder);
            }

            FitnessFDS fitnessFDS = new FitnessFDS(args.Z, 1);
            fitnessFDS.setOptimisable(!args.noOpt);

            SwMutSel mutSel = new SwMutSel(mutation, fitness, fitnessFDS);
            // SwMutSel mutSel = new SwMutSel(mutation, fitness);

            if (args.Zpenalty == 0) {
                mutSel.setDoPenalty(false);
            } else {
                mutSel.setDoPenalty(true);
                mutSel.setPenaltyShape(1.0, args.Zpenalty);
            }

            LikelihoodCalculator calculator = getCalculator(mutSel, fitness);

            SwMutSelFunction function = new SwMutSelFunction(calculator, mutSel);
            double[] initialGuess = function.getCurrentParameters();

            double likelihood;

            if (args.noOpt) {
                likelihood = function.value(initialGuess);
            } else {
                MultivariateOptimizer optimiser = new SimplexOptimizer(new SimpleValueChecker(-1, Constants.VALUE_CONVERGENCE_TOL));
                PointValuePair optima;

                optima = optimiser.optimize(
                        new ObjectiveFunction(function),
                        GoalType.MAXIMIZE,
                        new MaxEval(Constants.MAX_EVALUATIONS), // actual number of objective function evaluations
                        new MaxIter(Constants.MAX_ITERATIONS), // number of checks for convergence
                        new InitialGuess(initialGuess),
                        new NelderMeadSimplex(initialGuess.length)
                );

                likelihood = function.value(optima.getPoint());
            }
            Pair<Double, Double> nonSynAndSyn = mutSel.getRhoNonSynAndSyn();

            String result = String.format("%8.6f %s %8.6f %8.6f %8.6f %8.6f\n", likelihood, CoreUtils.join("%8.6f", " ", fitness.get()), fitnessFDS.get(), nonSynAndSyn.first, nonSynAndSyn.second, mutSel.getExpectedSubsPerSite());
            System.out.print(result);
            out = out + result;

            long end = System.currentTimeMillis();
            CoreUtils.msg("Time taken: %s\n", CoreUtils.getMinSecString(end - start));
        }

        try {
            Files.write(out, new File("swMutSel_FDS." + args.identifier), Charset.defaultCharset());
        } catch (IOException e) {
            System.out.println("Could not write results file!");
            e.printStackTrace();
        }
    }

    private LikelihoodCalculator getCalculator(SubstitutionModel model, Fitness fitness) {
        LikelihoodCalculator calculator;
        if (args.penalty == null) {
            calculator = new LikelihoodCalculator(args.tree, args.sites.column(args.site), model);
        } else {
            calculator = new PenalisedLikelihoodCalculator(args.tree, args.sites.column(args.site), model, args.penalty, fitness);
        }
        return calculator;
    }

}
