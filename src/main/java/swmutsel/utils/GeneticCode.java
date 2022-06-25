package swmutsel.utils;

import com.google.common.collect.Lists;
import com.google.common.primitives.Chars;
import com.google.common.primitives.Ints;

import java.util.List;

/**
 * Copied some of the Nucleotide and Codon functions we're using from PAL. Added some extra utility functions and
 * changed nucleotide order from ACGT to TCAG.
 *
 * @author Asif Tamuri (atamuri@nimr.mrc.ac.uk)
 */
public final class GeneticCode {

    private GeneticCode() {
        if (INSTANCE != null) {
            throw new IllegalStateException("Can't instantiate singleton.");
        }
    }

    public static GeneticCode getInstance() {
        return INSTANCE;
    }

    public static final int UNKNOWN_STATE = -1;
    public static final char[] AMINO_ACIDS = "ARNDCQEGHILKMFPSTWYV".toCharArray();

    // From http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    public static final String STANDARD_CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    public static final String VERTEBRATE_MITOCHONDRIAL_CODE = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG";
    public static final String PLASTID_CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    public static final String PSEUDOGENE_CODE = "FFLLSSSSYYFFCCFWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

    public static final int CODON_STATES = 64;
    public static final int AMINO_ACID_STATES = 20;
    public static final int NUCLEOTIDE_STATES = 4;

    private static final char UNKNOWN_CHARACTER = '?';
    private static final String UNKNOWN_TLA = "???";
    private static final char[] NUCLEOTIDES = "TCAG".toCharArray();
    private static final int UT_STATE = 0;
    private static final int C_STATE = 1;
    private static final int A_STATE = 2;
    private static final int G_STATE = 3;
    private static final char[][] CODONS_TLA_CHAR_ARRAY = new char[CODON_STATES][CODON_STATES];

    static {
        int i = 0;
        for (char n1 : NUCLEOTIDES) {
            for (char n2 : NUCLEOTIDES) {
                for (char n3 : NUCLEOTIDES) {
                    CODONS_TLA_CHAR_ARRAY[i++] = new char[]{n1, n2, n3};
                }
            }
        }
    }

    private static int[] SENSE_CODONS;

    private static final String[] CODONS_TLA = new String[CODON_STATES];

    static {
        int i = 0;
        for (char[] c : CODONS_TLA_CHAR_ARRAY) {
            CODONS_TLA[i++] = new String(c);
        }
    }

    // TODO: Stop codons should map to a 'stop' amino acid, 20, rather than -1! with fitness for analyses (MdR)
    private static final int[] CODONS_TO_AMINO_ACIDS = new int[CODON_STATES];
    private static final int[][] AMINO_ACIDS_TO_CODONS = new int[AMINO_ACID_STATES][];

    private static final GeneticCode INSTANCE = new GeneticCode();

    static {
        // By default, this instance uses the standard genetic code
        // Set your own code by calling setCode(String code)
        setCode(GeneticCode.STANDARD_CODE);
    }

    public String getCurrentCodeName() {
        StringBuilder aminoAcids = new StringBuilder();
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) {
            if (isSenseCodon(i)) {
                aminoAcids.append(getAminoAcidCharByIndex(getAminoAcidIndexFromCodonIndex(i)));
            } else {
                aminoAcids.append('*');
            }
        }

        String s = aminoAcids.toString();

        if (s.equals(STANDARD_CODE)) {
            return "standard";
        } else if (s.equals(VERTEBRATE_MITOCHONDRIAL_CODE)) {
            return "vertebrate_mit";
        } else if (s.equals(PLASTID_CODE)) {
            return "plastid";
        } else if (s.equals(PSEUDOGENE_CODE)) {
            return "pseudogene";
        } else {
            throw new RuntimeException("GeneticCode.getCurrentCodeName - ERROR: Unknown code!!");
        }
    }

    public synchronized static void setCode(String code) {
        // For each codon state, set the amino acid state
        for (int i = 0; i < CODON_STATES; i++) {
            if (code.charAt(i) == '*') {
                CODONS_TO_AMINO_ACIDS[i] = UNKNOWN_STATE; // stop codon
            } else {
                CODONS_TO_AMINO_ACIDS[i] = Chars.indexOf(AMINO_ACIDS, code.charAt(i)); // INSTANCE.getAminoAcidIndexByChar(code.charAt(i));
            }
        }

        // For each amino acid state, set all possible translating codon states
        for (int i = 0; i < AMINO_ACID_STATES; i++) {
            List<Integer> codonStates = Lists.newArrayList();

            // Search for this amino acid in the genetic code string
            int pos = 0;
            while (true) {
                pos = code.indexOf(AMINO_ACIDS[i], pos);
                if (pos == -1) {
                    break;
                } else {
                    codonStates.add(pos);
                }
                pos++;
            }
            AMINO_ACIDS_TO_CODONS[i] = Ints.toArray(codonStates);
        }

        SENSE_CODONS = getInstance().setSenseCodons();

        /*
        System.out.printf("swmutsel.utils.GeneticCode - Stop codons: ");
        for (int i = 0; i < CODON_STATES; i++) if (INSTANCE.isUnknownAminoAcidState(INSTANCE.getAminoAcidIndexFromCodonIndex(i))) System.out.printf("%s ", INSTANCE.getCodonTLA(i));
        System.out.println();
        System.out.printf("swmutsel.utils.GeneticCode - Methionine codon(s): ");
        for (int i : INSTANCE.getCodonIndexFromAminoAcidIndex(INSTANCE.getAminoAcidIndexByChar('M'))) System.out.printf("%s ", INSTANCE.getCodonTLA(i));
        System.out.println();
        */
    }

    public int getNucleotideIndexByChar(char c) {
        return Chars.indexOf(NUCLEOTIDES, c);
    }

    public int getAminoAcidIndexByChar(char c) {
        return Chars.indexOf(AMINO_ACIDS, c);
    }

    public char getNucleotideCharByIndex(int i) {
        return NUCLEOTIDES[i];
    }

    public char getAminoAcidCharByIndex(int i) {
        if (isUnknownAminoAcidState(i)) return UNKNOWN_CHARACTER;
        return AMINO_ACIDS[i];
    }

    public boolean isUnknownCodonState(int i) {
        return (i < 0 || i >= CODON_STATES);
    }

    public String getCodonTLA(int i) {
        if (isUnknownCodonState(i)) return UNKNOWN_TLA;
        return CODONS_TLA[i];
    }

    public char[] getNucleotidesFromCodonIndex(int i) {
        if (isUnknownCodonState(i)) return UNKNOWN_TLA.toCharArray();
        return CODONS_TLA_CHAR_ARRAY[i];
    }

    public int getCodonIndexFromNucleotides(char[] nuc) {
        String n = new String(nuc);
        for (int i = 0; i < CODONS_TLA.length; i++) {
            if (CODONS_TLA[i].equals(n)) {
                return i;
            }
        }
        return UNKNOWN_STATE;
    }

    public boolean isStopCodon(int codon) {
        return isUnknownAminoAcidState(getAminoAcidIndexFromCodonIndex(codon));
    }

    public boolean isSenseCodon(int codon) {
        return !isStopCodon(codon);
    }

    public int[] setSenseCodons() {
        // This can be set once the genetic code is set.
        List<Integer> senseCodons = Lists.newArrayList();
        for (int i = 0; i < GeneticCode.CODON_STATES; i++) if (!GeneticCode.getInstance().isStopCodon(i)) senseCodons.add(i);
        return Ints.toArray(senseCodons);
    }

    public int[] getSenseCodons() {
        return SENSE_CODONS;
    }

    public int getAminoAcidIndexFromCodonIndex(int i) {
        if (isUnknownCodonState(i)) return UNKNOWN_STATE;
        return CODONS_TO_AMINO_ACIDS[i];
    }

    public int[] getCodonIndexFromAminoAcidIndex(int i) {
        if (isUnknownAminoAcidState(i)) return new int[]{UNKNOWN_STATE};
        return AMINO_ACIDS_TO_CODONS[i];
    }

    public boolean isUnknownAminoAcidState(int i) {
        return (i < 0 || i >= AMINO_ACID_STATES);
    }

    public boolean isTransitionByChar(char firstChar, char secondChar) {
        return isTransitionByState(getNucleotideIndexByChar(firstChar), getNucleotideIndexByChar(secondChar));
    }

    public boolean isTransitionByState(int firstState, int secondState) {
        switch (firstState) {
            case A_STATE: {
                return secondState == G_STATE;
            }
            case C_STATE: {
                return secondState == UT_STATE;
            }
            case G_STATE: {
                return secondState == A_STATE;
            }
            case UT_STATE: {
                return secondState == C_STATE;
            }
            default:
                return false;
        }
    }
}
