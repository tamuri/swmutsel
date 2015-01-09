package swmutsel.model.parameters;

import com.google.common.collect.Maps;

import java.util.Collection;
import java.util.Map;
import java.util.Set;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 16/10/2013 13:58
 */
public class FitnessStore extends Parameter implements Map<Integer,Fitness> {

    private static final long serialVersionUID = 770459384773364802L;
    private final Map<Integer, Fitness> store = Maps.newTreeMap();

    // Parameter abstract methods
    @Override
    public double[] getOptimisable() {
        return new double[0];
    }

    @Override
    public void setOptimisable(double[] params) {

    }

    @Override
    public int getOptimisableCount() {
        return 0;
    }

    // Map<K, V> interface
    @Override
    public int size() {
        return store.size();
    }

    @Override
    public boolean isEmpty() {
        return store.isEmpty();
    }

    @Override
    public boolean containsKey(Object key) {
        return store.containsKey(key);
    }

    @Override
    public boolean containsValue(Object value) {
        return store.containsKey(value);
    }

    @Override
    public Fitness get(Object key) {
        if (!this.store.containsKey(key)) {
            throw new RuntimeException("FitnessStore does not have entry for " + key);
        }
        return store.get(key);
    }

    @Override
    public Fitness put(Integer key, Fitness value) {
        return store.put(key, value);
    }

    @Override
    public Fitness remove(Object key) {
        return store.remove(key);
    }

    @Override
    public void putAll(Map<? extends Integer, ? extends Fitness> m) {
        store.putAll(m);
    }

    @Override
    public void clear() {
        store.clear();
    }

    @Override
    public Set<Integer> keySet() {
        return store.keySet();
    }

    @Override
    public Collection<Fitness> values() {
        return store.values();
    }

    @Override
    public Set<Entry<Integer, Fitness>> entrySet() {
        return store.entrySet();
    }
}
