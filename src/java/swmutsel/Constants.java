package swmutsel;

/**
 * Various constants used throughout the program but not expected to be changed by end-users.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class Constants {
    /**
     * One of the residues (usually the most observed) is fixed to this value; other fitnesses are relative to this one.
     */
    public static final double FITNESS_FIXED_FOR_RELATIVE = 0.0;

    /**
     * Random initial value for fitness parameter from range -RANDOM_INITIAL_FITNESS_RANGE to RANDOM_INITIAL_FITNESS_RANGE.
     */
    public static final int RANDOM_INITIAL_FITNESS_RANGE = 3;

    /**
     * When to assume convergence of MLE of fitness parameters (log-likelihood).
     */
    public static double CONVERGENCE_TOL = 5.0e-6;
    public static final double SMALL_CHANGE = 1e-3;
    public static final double MINIMUM_TOL = 1.0e-7;

    /**
     * Brent optimiser settings, used for one-at-a-time branch length optimisation
     */
    public static final double BRENT_REL = 1.0e-10;
    public static final double BRENT_ABS = 1.0e-14;

    /**
     * Effectively -20 = -Infinity and 20 = Infinity for fitness parameter.
     */
    public static final double FITNESS_BOUND = 20;

    /**
     * Terminate the optimisation routine after this number of evaluations.
     */
    public static final int MAX_EVALUATIONS = 12000;
    // For SimplexOptimiser, iterations ~ 0.75 * evaluations
    public static final int MAX_ITERATIONS = 6000;

    /**
     * To constrain MLE of fitness within FITNESS_BOUND, return this log-likelihood if we step outside of bounds.
     */
    public static final double VERY_BAD_LIKELIHOOD = Double.NEGATIVE_INFINITY;

    /**
     * Scaling conditional probabilities to avoid underflow. See
     * Yang, Z. (2000). Journal of Molecular Evolution, 51(5), 423–432.
     */
    public static final double SCALING_THRESHOLD = (1.0 / Math.pow(2, 128));
    public static final int SCALING_NODE_STEP = 1;

    /**
     * How to split the branch connecting heterogeneous models (e.g. 0.5 = half-way).
     */
    public static final double CLADE_BRANCH_SPLIT = 0.5;

    /**
     * For branch length optimisation constraints
     */
    public static final double INITIAL_BRANCH_LENGTH = 0.1;
    public static final double MIN_BRANCH_LENGTH = 1e-12;
    public static final double MAX_BRANCH_LENGTH = 10;
    public static final int MAX_BRANCH_LENGTH_SMALL_CHANGES = 5;

    /**
     * The filenames for parsed results files.
     */
    public static final String S_FILENAME = "S.txt";
    public static final String Q0_FILENAME = "Q0.txt";
    public static final String QS_FILENAME = "QS.txt";
    public static final String PI_FILENAME = "PiS.txt";
    public static final String PIAA_FILENAME = "PiAA.txt";
    public static final String F_FILENAME = "F.txt";

    public static final String RESOURCE_PATH = "/resources/";
    public static final String HELP_DESCRIPTION = "description.txt";
    public static final String HELP_COMMAND_FMUTSEL0 = "command_fmutsel0.txt";
    public static final String HELP_COMMAND_SIMULATE = "command_simulate.txt";
    public static final String HELP_COMMAND_SWMUTSEL = "command_swmutsel.txt";
    public static final String HELP_FULL = "help.txt";

}
