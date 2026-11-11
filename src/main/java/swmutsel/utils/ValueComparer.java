package swmutsel.utils;

import java.util.Comparator;
import java.util.Map;

public class ValueComparer<K, V extends Comparable<V>> implements Comparator<K>{
    /**
     * Little utility class that can be used to construct a SortedMap that sorts
     * keys by values, descending. It is type-safe.
     * <p/>
     * e.g.
     * <p/>
     * Map<K, V> unsorted = new HashMap<K,V>();
     * SortedMap<K, V> sorted = new TreeMap<K, V>(new ValueComparer<K, V>(unsorted));
     * sorted.putAll(unsorted);
     *
     * @param <K>
     * @param <V>
     */
        private final Map<K, V> map;

    public ValueComparer(Map<K, V> map) {
            this.map = map;
        }

        // Compare two values in a map (in descending order)
        public int compare(K key1, K key2) {
            V value1 = this.map.get(key1);
            V value2 = this.map.get(key2);
            int c = value2.compareTo(value1);
            if (c != 0) {
                return c;
            }
            Integer hashCode1 = key1.hashCode();
            Integer hashCode2 = key2.hashCode();
            return hashCode1.compareTo(hashCode2);
        }
    }