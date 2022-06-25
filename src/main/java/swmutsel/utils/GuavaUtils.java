package swmutsel.utils;

import com.google.common.collect.ArrayTable;
import com.google.common.collect.Table;

import java.util.Collection;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 07/11/2013 10:50
 */
public class GuavaUtils {
    public static <R, C, V> Table<R, C, V> getColumnSubset(Table<R, C, V> table, Collection<C> columns) {

        Collection<R> rows = table.rowKeySet();

        Table<R, C, V> subTable = ArrayTable.create(rows, columns);

        for (R row : rows) {
            for (C column : columns) {
                subTable.put(row, column, table.get(row, column));
            }
        }

        return subTable;

    }
}
