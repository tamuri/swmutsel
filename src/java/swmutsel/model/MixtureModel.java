package swmutsel.model;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 02/04/2014 15:14
 */
public abstract class MixtureModel extends SubstitutionModel {

    public abstract double getWeight(int mixture);

    public abstract SubstitutionModel[] getModels();

    @Override
    public TransProbCalculator getPtCalculator() {
        return null;
    }

    @Override
    public double[] getCodonFrequencies() {
        return new double[0];
    }

}
