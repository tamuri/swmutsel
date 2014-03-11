package swmutsel.model.parameters;

import java.util.Arrays;
import java.util.List;

/**
 * Maps a list of parameters to an array of [optimisable] doubles
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 31/10/2013 10:40
 */
public class Mapper {

    private static int getArraySize(List<Parameter> parameters) {
        int arraySize = 0;
        for (Parameter p : parameters) {
            if (p.isOptimisable()) arraySize += p.getOptimisableCount();
        }
        return arraySize;
    }

    public static double[] getOptimisable(List<Parameter> parameters) {
        double[] optimisable = new double[getArraySize(parameters)];

        int pos = 0;
        for (Parameter p : parameters) {
            if (p.isOptimisable()) {
                double[] opt = p.getOptimisable();
                for (double d : opt) optimisable[pos++] = d;
            }
        }

        return optimisable;
    }

    public static void setOptimisable(List<Parameter> parameters, double[] point) {
        int pos = 0;
        for (Parameter p : parameters) {
            if (p.isOptimisable()) {
                int size = p.getOptimisableCount();
                p.setOptimisable(Arrays.copyOfRange(point, pos, pos + size));
                pos += size;
            }
        }
    }
}
