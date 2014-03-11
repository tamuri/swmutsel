package swmutsel.model.parameters;

import java.util.TreeMap;

/**
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 13:58
 */
public class FitnessStore extends TreeMap<Integer, Fitness> {

    private static final long serialVersionUID = 770459384773364802L;

    public void set(int site, Fitness f) {
        super.put(site, f);
    }

}
