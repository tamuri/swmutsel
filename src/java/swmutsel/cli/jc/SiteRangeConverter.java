package swmutsel.cli.jc;

import com.beust.jcommander.IStringConverter;
import com.google.common.collect.ImmutableRangeSet;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 19/11/2013 13:57
 */
public class SiteRangeConverter implements IStringConverter<RangeSet<Integer>> {

    private static final String SEP = "-";
    @Override
    public RangeSet<Integer> convert(String s) {

        RangeSet<Integer> ranges = TreeRangeSet.create();

        String[] parts = s.split(",");

        for (String part : parts) {
            int pos = part.indexOf(SEP);

            // if this is a range
            if (pos >= 0) {
                String[] split = part.split(SEP);

                ranges.add(Range.closed(Integer.parseInt(split[0]), Integer.parseInt(split[1])));


            } else {
                // it's a single site
                ranges.add(Range.singleton(Integer.parseInt(part)));
            }
        }

        return ImmutableRangeSet.copyOf(ranges);

    }
}
