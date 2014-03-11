package swmutsel.options;

import com.beust.jcommander.IStringConverter;
import com.beust.jcommander.ParameterException;
import swmutsel.model.DirichletPenalty;
import swmutsel.model.MVNPenalty;
import swmutsel.model.Penalty;
import swmutsel.utils.CoreUtils;

public class PenaltyConverter implements IStringConverter<Penalty> {
    @Override
    public Penalty convert(String s) {
        try {
            Penalty p;

            String[] definition = s.split(",");
            definition[0] = definition[0].toLowerCase();

            if (MVNPenalty.LABEL.equals(definition[0])) {
                p = new MVNPenalty(Double.parseDouble(definition[1]));
            } else if (DirichletPenalty.LABEL.equals(definition[0])) {
                p = new DirichletPenalty(Double.parseDouble(definition[1]));
            }  else {
                throw new ParameterException("Could not create penalty '" + s + "'.\n");
            }

            CoreUtils.msg("swmutsel.options.PenaltyConverter - Using %s.\n", p.toString());
            return p;

        } catch (Exception e) {
            throw new ParameterException("Could not create penalty '" + s + "'.\n");
        }
    }
}
