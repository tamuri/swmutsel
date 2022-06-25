package swmutsel.utils;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.primitives.Doubles;

import java.io.IOException;
import java.net.ServerSocket;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class CoreUtils {

    public static void msg(String string, Object... args) {
        System.out.printf("%s - ", new SimpleDateFormat("yy-MM-dd HH:mm:ss:SSS").format(new Date()));
        System.out.printf(string, args);
    }

    public static void print2D(double[][] array, boolean flat) {
        if (flat) {
            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[i].length; j++) {
                    System.out.printf("%s,%s = %s\n", i, j, array[i][j]);
                }
            }
        } else {
            for (int i = 0; i < array.length; i++) {
                for (int j = 0; j < array[i].length; j++) {
                    System.out.printf("%s\t", array[i][j]);
                }
                System.out.println();
            }
        }
    }

    public static double sum(double[] list) {
        double t = 0.0;
        for (double d : list) {
            t += d;
        }
        return t;
    }

    public static double sum(Collection<Double> list) {
        double s = 0.0;
        for (double d : list) s += d;
        return s;
    }

    public static double[] alr(double[] p) {
        // additive log-ratio transformation
        // y = log(p1 / pn), log(p2 / pn), log(p3 / pn) ... log(pn-1 / pn)

        double sum = 0;
        for (double d : p) {
            sum += d;
        }

        double[] y = new double[p.length];

        for (int i = 0; i < p.length; i++) {
            y[i] = Math.log(p[i] / (1 - sum));
        }

        return y;
    }

    public static double[] alr_inv(double[] y) {
        // inverse additive log-ratio transformation
        // p_i = exp(y_i) / (1 + sum(exp(y_i)))
        double[] p = new double[y.length];

        double sum = 0;
        for (double d : y) {
            sum += Math.exp(d);
        }

        for (int i = 0; i < y.length; i++) {
            p[i] = Math.exp(y[i]) / (1 + sum);
        }

        return p;
    }

    public static int findFreePort() {
        ServerSocket socket = null;
        try {
            socket = new ServerSocket(0);
            socket.setReuseAddress(true);
            int port = socket.getLocalPort();
            try {
                socket.close();
            } catch (IOException e) {
// Ignore IOException on close()
            }
            return port;
        } catch (IOException e) {
        } finally {
            if (socket != null) {
                try {
                    socket.close();
                } catch (IOException e) {
                }
            }
        }
        throw new IllegalStateException("Could not find a free TCP/IP port.");
    }

   public static double[] rep(double x, int n) {
       double[] o = new double[n];
       for (int i = 0; i < n; i++) {
           o[i] = x;
       }
       return o;
   }

   public static int[] rep(int x, int n) {
       int[] o = new int[n];
       for (int i = 0; i < n; i++) {
           o[i] = x;
       }

       return o;
   }

   public static double[] seqd(double start, double end, double step) {
       List<Double> out = new ArrayList<Double>();

       out.add(start);

       double n = start;
       n += step;

       while (n <= end) {
           out.add(n);
           n += step;
       }

       double[] o = Doubles.toArray(out);

       return o;
   }

   public static int[] seq(int start, int end) {
       int[] o = new int[end - start + 1];
       for (int i = 0; i < o.length; i++) {
           o[i] = start + i;
       }
       return o;
   }

   public static <T> List<T> getFutureResults(List<Future<T>> futures) {
       List<T> results = Lists.newArrayList();

       for (Future<T> f : futures) {
           try {
               results.add(f.get());
           } catch (Exception e) {
               e.printStackTrace();
               throw new RuntimeException(e);
           }
       }

       return results;
   }

   public static <K, V> LinkedHashMap<K, V> getLinkedHashMap(Pair<K, V>... pairs) {
       LinkedHashMap<K, V> map = Maps.newLinkedHashMap();

       for (Pair<K, V> pair : pairs) {
           map.put(pair.first, pair.second);
       }

       return map;
   }

   public static double calculateMachineEpsilonDouble() {
       double machEps = 1.0d;

       do
           machEps /= 2.0d;
       while ((1.0 + (machEps / 2.0)) != 1.0);

       return machEps;
   }

   public static String convertStreamToString(java.io.InputStream is) {
        java.util.Scanner s = new java.util.Scanner(is).useDelimiter("\\A");
        return s.hasNext() ? s.next() : "";
    }

    public static String join(String format, String join, double[] array) {
        if (array.length == 0) return "";
        StringBuilder sb = new StringBuilder();
        sb.append(String.format(format, array[0]));
        for (int i = 1; i < array.length; i++) {
            sb.append(join);
            sb.append(String.format(format, array[i]));
        }
        return sb.toString();
    }

    public static Throwable getRootCause(Throwable throwable) {
        if (throwable.getCause() != null)
            return getRootCause(throwable.getCause());

        return throwable;
    }

    public static String getMinSecString(long millis) {
        return String.format("%01d:%02d",
                TimeUnit.MILLISECONDS.toMinutes(millis),
                TimeUnit.MILLISECONDS.toSeconds(millis) -
                        TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis))
        );
    }

}