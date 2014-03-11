package swmutsel.options;

import com.beust.jcommander.IStringConverter;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class DoubleArrayConverter implements IStringConverter<double[]> {

  @Override
  public double[] convert(String value) {
      String[] s = value.split(",");
      double[] d = new double[s.length];
      for (int i = 0; i < s.length; i++) {
          d[i] = Double.parseDouble(s[i]);
      }
      return d;
  }

}
