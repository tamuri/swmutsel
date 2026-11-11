package swmutsel.results;

import com.google.common.base.Charsets;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.io.Files;
import com.google.common.primitives.Doubles;
import swmutsel.Constants;
import swmutsel.utils.GeneticCode;

import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 * Given the full output of the TdG12 swMutSel0 analysis, writes a file of fitnesses, ordered by site
 * <p/>
 * TODO: Handle heterogeneous fitness output (e.g. Fitness_C1, Fitness_C2 etc.) Currently, you have to do those by hand
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class FitnessExtractor {

    public void extract(List<String> _fitness) throws Exception {

        Iterator<String> items = _fitness.iterator();

        Map<Integer, double[]> fitnesses = Maps.newHashMap();

        while (items.hasNext()) {

            int site = Integer.parseInt(items.next());

            double[] f = new double[GeneticCode.AMINO_ACID_STATES];
            for (int i = 0; i < GeneticCode.AMINO_ACID_STATES; i++) {
                f[i] = Double.parseDouble(items.next());
            }

            fitnesses.put(site, f);

        }


        ArrayList<Integer> orderedKeys = Lists.newArrayList(fitnesses.keySet());
        Collections.sort(orderedKeys);

        BufferedWriter writer = Files.newWriter(new File(Constants.F_FILENAME), Charsets.US_ASCII);

        for (int site : orderedKeys) {
            writer.write(Doubles.join(" ", fitnesses.get(site)));
            writer.write("\n");
        }

        writer.close();

    }

}
