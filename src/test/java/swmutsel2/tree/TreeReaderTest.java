package swmutsel2.tree;

import org.junit.Test;
import org.junit.Before;
import org.junit.After;

import java.io.File;
import java.net.URL;

import static org.junit.Assert.assertEquals;

/**
 * TreeReader Tester.
 *
 * @author Asif Tamuri
 * @version 1.0
 */
public class TreeReaderTest {
    private static final String TREE_STRING = "((Homo_sapiens: 0.461856, Macaca_mulatta: 0.762943): 0.806366, Mus_musculus: 1.262486, Bos_taurus: 0.744734);";
    private static final String TREE_FILE = "mit4s.tree";
    private static final int TREE_TIP_COUNT = 4;


    private static final String EXPECTED_TOSTRING = "Nodes: 6\n" +
            "Branches: [\n" +
            "\t(Homo_sapiens(0), (0)):0.461856\n" +
            "\t(Macaca_mulatta(1), (0)):0.762943\n" +
            "\t((0), (1)):0.806366\n" +
            "\t(Mus_musculus(2), (1)):1.262486\n" +
            "\t(Bos_taurus(3), (1)):0.744734\n" +
            "]\n" +
            "Post-order traversal (internal nodes): [(0), (1)]\n" +
            "Root: (1)\n" +
            "Incoming nodes to root: 3 [(0), Mus_musculus(2), Bos_taurus(3)]\n";

    @Before
    public void before() throws Exception {
    }

    @After
    public void after() throws Exception {
    }

    /**
     *
     * Method: read(String stringOrFile)
     *
     */
    @Test
    public void testReadString() throws Exception {
        Tree stringOut = TreeReader.read(TREE_STRING);
        assertEquals("Tree from string returned different toString()", EXPECTED_TOSTRING, stringOut.toString());
        assertEquals("Tree from string returned incorrect number of taxa.", TREE_TIP_COUNT, stringOut.getTipCount());
    }

    @Test
    public void testReadFile() throws Exception {
        URL tree = this.getClass().getResource("/" + TREE_FILE);
        String path  = new File(tree.getFile()).getAbsolutePath();

        Tree fileOut = TreeReader.read(path);
        assertEquals("Tree from file returned different toString()", EXPECTED_TOSTRING, fileOut.toString());
        assertEquals("Tree from file returned incorrect number of taxa.", TREE_TIP_COUNT, fileOut.getTipCount());
    }

}

