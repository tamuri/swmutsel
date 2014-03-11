package swmutsel;

import com.beust.jcommander.internal.Lists;
import com.google.common.base.Joiner;
import com.google.common.io.Files;
import swmutsel.model.parameters.FitnessStore;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.Collections;
import java.util.List;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 01/11/2013 15:38
 */
public class ArgumentWriter {
    public static void writeFile(Arguments arguments, String filename) {
        saveStringToFile(getOutputString(arguments), filename);
    }

    private static String getOutputString(Arguments a) {
        StringBuilder s = new StringBuilder();
        stringField(s, "-sequences", a._alignmentPath);
        stringField(s, "-geneticcode", GeneticCode.getInstance().getCurrentCodeName());


        if (a.fix.contains("branches") || a.fix.contains("ALL")) {
            stringField(s, "-tree", a.treePath);
        } else {
            stringField(s, "-tree", a.tree.toString());
        }

        stringField(s, "-kappa", String.format("%.7f", a.kappa));
        stringField(s, "-pi", String.format("%.7f,%.7f,%.7f,%.7f", a.pi[0],a.pi[1],a.pi[2],a.pi[3]));
        stringField(s, "-scaling", String.format("%.7f", a.scaling));

        if (a.penalty != null)
            stringField(s, "-penalty", a.penalty.toString());

        if (a.fix.size() > 0)
            stringField(s, "-fix", Joiner.on(',').join(a.fix));

        stringFitnesses(s, a.fitnesses);
        s.append("\n\n");
        s.append("//\n\n");
        return s.toString();
    }

    private static void saveStringToFile(String string, String file) {
        try {
            BufferedWriter writer = Files.newWriter(new File(file), Charset.defaultCharset());
            writer.write(string);
            writer.close();
        } catch (IOException e) {
            CoreUtils.msg("WARNING: Could not write checkpoint file.");
        }
    }

    private static void stringFitnesses(StringBuilder s, FitnessStore fitnesses) {
        s.append("-fitness");
        s.append("\n");

        List<Integer> sites = Lists.newArrayList(fitnesses.keySet());
        Collections.sort(sites);

        for (int site : sites) {
            s.append(site);
            s.append(",");
            double[] f = fitnesses.get(site).get();
            s.append(String.format("%.7f", f[0]));
            for (int i = 1; i < f.length; i++) s.append(String.format(",%.7f", f[i]));
            s.append("\n");
        }

        s.append("\n");
    }

    private static void stringField(StringBuilder s, String field, String value) {
        s.append(field);
        s.append("\n");
        s.append(value);
        s.append("\n\n");
    }
}
