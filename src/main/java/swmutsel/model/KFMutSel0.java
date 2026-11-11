package swmutsel.model;

import swmutsel.model.parameters.*;
import swmutsel.utils.CoreUtils;

import java.io.Serializable;

/**
 * java -Djava.library.path=lib -jar dist/swmutsel.jar -fmutsel0 -t mtCDNApri_unrooted.trees -gc vertebrate_mit -n test -s ../data/DatingSoftBound/mtCDNApri.txt -classes 3 -T 3
 * -kappa 20.3900008 -omega 0.2214758 -pi 0.1198626,0.3778540,0.4409021,0.0613813 -weights 0.4701069,0.3163964,0.2134967 -fitness 0.0000000,-5.2304959,-3.3770733,-5.0074518,2.1597377,-2.6939260,-3.0085947,4.4494463,0.5961836,2.9330161,2.9894690,0.2446834,0.9894206,1.9327684,-2.7000194,-0.6112096,-0.6218997,3.8835224,-1.6739847,3.4467583 0.0000000,-9.3377893,-1.7745342,-3.2675760,-4.1672616,-2.0401938,-0.2889405,-3.8835205,-4.8068900,-2.9149990,-5.3402425,-7.4313300,-0.5130639,1.1547381,-1.7395473,-4.1009040,-5.2974313,-6.1146997,-0.5793016,-5.4065138 0.0000000,0.1125513,-3.1932192,0.5988166,-2.7294258,-4.4594118,-13.8024252,-1.2143821,-3.6716985,-1.8868200,-2.3354136,-4.4590942,-1.6727464,-5.3611361,-5.1967788,0.4392429,-0.3959858,-2.9086462,-3.6215930,-1.5059550
 * -28828.08886
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 15:24
 */
public class KFMutSel0 extends MixtureModel implements Serializable {

    private static final long serialVersionUID = 7703558036907403610L;

    // The parameters of the model
    private final TsTvRatio kappa;
    private final Omega omega;
    private final BaseFrequencies pi;
    private final Probabilities weights;
    private final FitnessStore fitnesses;

    private final FMutSel0[] mixtures;

    public KFMutSel0(TsTvRatio kappa, Omega omega, BaseFrequencies pi, Probabilities weights, FitnessStore fitnesses) {
        if (weights.get().length != fitnesses.values().size())
            throw new RuntimeException("ERROR: KFMutSel0 num. of weights (" + weights.get().length + ") != num. of fitnesses (" + fitnesses.values().size() +").");

        this.kappa = kappa;
        this.omega = omega;
        this.pi = pi;
        this.weights = weights;
        this.fitnesses = fitnesses;
        this.mixtures = new FMutSel0[this.fitnesses.values().size()];
        for (int i = 0; i < this.mixtures.length; i++) {
            this.mixtures[i] = new FMutSel0(this.kappa, this.omega, this.pi, this.fitnesses.get(i + 1));
        }

        super.clearParameters();
        super.addParameters(pi, kappa, omega, weights);
        super.addParameters(this.fitnesses.values());

        build();
    }

    @Override
    public void build() {
        double totalRate = 0;
        for (int i = 0; i < this.mixtures.length; i++) {
            this.mixtures[i].setScale(0); // do not scale
            this.mixtures[i].setSenseCodonFrequencies();
            totalRate += this.mixtures[i].getAverageRate() * this.weights.get()[i];
        }

        // normalise each Q in the mixture with scaling so average across all
        // Q classes is 1
        for (FMutSel0 mixture : this.mixtures) {
            mixture.setScale(totalRate);
            mixture.build();
        }
    }


    @Override
    public boolean parametersValid() {
        for (FMutSel0 m : this.mixtures) {
            if (!m.parametersValid()) return false;
        }
        return true;
    }

    @Override
    public double getWeight(int mixture) {
        return this.weights.get()[mixture];
    }

    @Override
    public SubstitutionModel[] getModels() {
        return this.mixtures;
    }

    @Override
    public String toString() {
        String s = String.format("KFMutSel0 { -kappa %.7f -omega %.7f -pi %s -weights %s ",
                kappa.get(),
                omega.get(),
                CoreUtils.join("%.7f", ",", pi.get()),
                CoreUtils.join("%.7f", ",", weights.get()));

        for (Fitness f : fitnesses.values()) {
            s += "-fitness " + CoreUtils.join("%.7f",",",f.get()) + " ";
        }

        s += "}";

        return s;
    }
}
