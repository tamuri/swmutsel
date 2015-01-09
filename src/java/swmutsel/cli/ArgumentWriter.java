package swmutsel.cli;

import com.google.common.collect.Lists;
import com.google.common.io.Files;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 01/11/2013 15:38
 */
public class ArgumentWriter {
    public static void saveStringToFile(String string, String file) {
        try {
            BufferedWriter writer = Files.newWriter(new File(file), Charset.defaultCharset());
            writer.write(string);
            writer.close();
        } catch (IOException e) {
            CoreUtils.msg("WARNING: Could not write checkpoint file.");
        }
    }

    public static void stringFitnesses(StringBuilder s, FitnessStore fitnesses, boolean withIndex) {
        s.append("-fitness");
        s.append("\n");

        List<Integer> sites = Lists.newArrayList(fitnesses.keySet());
        Collections.sort(sites);

        for (int site : sites) {
            if (withIndex) {
                s.append(site);
                s.append(",");
            }
            double[] f = fitnesses.get(site).get();
            s.append(CoreUtils.join("%.7f", ",", f));
            s.append("\n");
        }

        s.append("\n");
    }

    public static void stringField(StringBuilder s, String field, String value) {
        s.append(field);
        s.append("\n");

        if (value != null) {
            s.append(value);
            s.append("\n");
        }

        s.append("\n");
    }

    public static void writeSiteInformation(SwMutSelArguments args, Map<Integer, Double> siteLogLikelihood, Map<Integer, Pair<Double, Double>> siteDnDs, String file) {
        try {
            // append to the file
            BufferedWriter writer = new BufferedWriter(new FileWriter(new File(file), true));
            writer.write("site\tlnL\trel_w\trel_g");
            writer.newLine();

            List<Integer> sites = Lists.newArrayList(args.patternSiteMap.keySet());
            Collections.sort(sites);

            for (int site : sites) {
                double logLikelihood = siteLogLikelihood.get(args.patternSiteMap.get(site));
                Pair<Double, Double> dnDs = siteDnDs.get(args.patternSiteMap.get(site));
                writer.write(String.format("%s\t%.7f\t%.7f\t%.7f", site, logLikelihood, dnDs.first, dnDs.second));
                writer.newLine();
            }

            writer.newLine();
            writer.close();

        } catch (IOException e) {
            throw new RuntimeException("ERROR: Could not write to " + file);
        }
    }
}
