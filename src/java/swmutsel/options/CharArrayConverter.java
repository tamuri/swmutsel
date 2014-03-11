package swmutsel.options;

import com.beust.jcommander.IStringConverter;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class CharArrayConverter implements IStringConverter<char[]> {

    @Override
    public char[] convert(String value) {
        String[] s = value.split(",");
        char[] c = new char[s.length];
        for (int i = 0; i < s.length; i++) {
            c[i] = s[i].toCharArray()[0];
        }
        return c;
    }

}
