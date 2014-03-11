package swmutsel.options;

import com.beust.jcommander.IStringConverter;
import pal.alignment.Alignment;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.PhyloUtils;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 13:53
 */
public class AlignmentConverter implements IStringConverter<Alignment> {
    @Override
    public Alignment convert(String value) {
        CoreUtils.msg("Reading alignment: %s\n", value);
        return PhyloUtils.readAlignment(value);
    }
}
