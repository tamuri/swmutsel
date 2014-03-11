package swmutsel.utils;

import com.google.common.collect.Maps;

import java.util.Map;
import java.util.concurrent.atomic.AtomicLong;

/**
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public class CodeTimer {
    private static final Map<String, AtomicLong> times = Maps.newConcurrentMap();

    public static void store(String key, long startTime) {
        if (!times.containsKey(key)) {
            times.put(key, new AtomicLong());
        }
        times.get(key).getAndAdd(System.currentTimeMillis() - startTime);
    }

    public static long start() {
        return System.currentTimeMillis();
    }

    public static void printAll() {
        System.out.println("\n-------------------\nswmutsel.utils.CodeTimer\n-------------------");
        for (Map.Entry e : times.entrySet()) {
            System.out.printf("%s = %sms\n", e.getKey(), e.getValue());
        }
    }
}
