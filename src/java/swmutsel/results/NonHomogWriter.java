package swmutsel.results;

import swmutsel.utils.GeneticCode;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Random;

/**
 * Author: Asif Tamuri (atamuri@ebi.ac.uk)
 * Date: 26/03/2013 15:00
 */
public class NonHomogWriter {
    boolean doMuts;
    boolean doNonSyn;

    public static void main(String[] args) throws Exception {
        NonHomogWriter ads = new NonHomogWriter();
        ads.doMuts = Boolean.parseBoolean(args[0]);
        ads.doNonSyn = Boolean.parseBoolean(args[1]);
        ads.run();
    }

    private void run() throws Exception {

        //boolean doMuts = true; // otherwise, substitutions
        //boolean doNonSyn = true; // otherwise, synonymous and non-synonymous changes

        // mitochondrial genome
        //char[] aaCode = {	'F', 'F', 'L', 'L', 'S', 'S', 'S', 'S', 'Y', 'Y', '*', '*', 'C', 'C', 'W', 'W', 'L', 'L', 'L', 'L', 'P', 'P', 'P', 'P', 'H', 'H', 'Q', 'Q', 'R', 'R', 'R', 'R', 'I', 'I', 'M', 'M', 'T', 'T', 'T', 'T', 'N', 'N', 'K', 'K', 'S', 'S', '*', '*', 'V', 'V', 'V', 'V', 'A', 'A', 'A', 'A', 'D', 'D', 'E', 'E', 'G', 'G', 'G', 'G'};
        // char[] aaCode = GeneticCode.VERTEBRATE_MITOCHONDRIAL_CODE.toCharArray();
        // TODO don't forget about sites below!!
        char[] aaCode = GeneticCode.STANDARD_CODE.toCharArray();
        //char[] aaCode = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG".toCharArray();


        //String dir = "/Users/atamuri/Documents/work/mitochondria/110318_TdG_Mit_ApproxVsFull/64x64/";
        //String dir = "/Users/atamuri/Documents/work/mitochondria/110328_Mit_PenalisedLnL_MaxEnt/64x64/";
        //String dir = "/Users/atamuri/Documents/work/mitochondria/110329_TdG_PB2_FullMaxEnt/pb2.exact.maxent/avian.equi/";
        //String dir = "/Users/atamuri/Documents/work/mitochondria/110329_TdG_PB2_FullMaxEnt/pb2.exact.maxent/human.adap/";
        //String dir = "/Users/atamuri/Documents/work/mitochondria/paper/response/flu.sim/9.distS/3/avian.equi/";
        // String dir = "/Users/atamuri/Documents/work/mitochondria/paper/response/flu.sim/9.distS/3/human.adap/";
        //String dir = "/Users/atamuri/Documents/work/mitochondria/paper/response/consistency/20.AA.equal.F/";
        // String dir = "./";
        String dir = "./";

        double hi = 21D;
        double low = -21D;
        int numbins = doMuts ? 168 : 167; // 40 / 0.25 ....168 for mutations, 167 for subs (to keep symmetric)

        int sites = 759; // TODO: change this when you run it!!!

        double[] bins = new double[numbins];
        double[][] kbins = new double[sites][numbins];
        double[] nonsigbins = new double[numbins];
        double[] sigbins = new double[numbins];


        BufferedReader buffS = new BufferedReader(new FileReader(dir + "S.txt"));
        BufferedReader buffQ = new BufferedReader(new FileReader(dir + "QS.txt"));
        BufferedReader buffPi = new BufferedReader(new FileReader(dir + "PiS.txt"));
        BufferedReader buffQ0 = new BufferedReader(new FileReader(dir + "Q0.txt"));

        String[] stringQ0 = buffQ0.readLine().split("\\s");
        double[] q0 = new double[4096];
        for (int iPair = 0; iPair < 4096; iPair++) {
            q0[iPair] = Double.parseDouble(stringQ0[iPair]);
        }
        buffQ0.close();

        double denom = 0D;
        for (int i = 0; i < sites; i++) {
            String lineQ = buffQ.readLine();
            String linePi = buffPi.readLine();
            String[] stringQ = lineQ.split("\\s");
            String[] stringPi = linePi.split("\\s");

            double siteSum = 0D;
            for (int iPair = 0; iPair < 4096; iPair++) {
                int iCodon = (iPair - iPair%64)/64;
                char iRes = aaCode[iCodon];
                char jRes = aaCode[iPair%64];
                double piValue = Double.parseDouble(stringPi[iCodon]);
                double qValue = Double.parseDouble(stringQ[iPair]);

                boolean isCodonChange = iCodon != (iPair % 64);
                boolean isNonSynChange = iRes != jRes;

                //if (iCodon != (iPair % 64)) {
                if (doNonSyn ? isNonSynChange : isCodonChange) {
                    // if ( (iRes != jRes) ) {
                    siteSum += piValue * (doMuts ? q0[iPair] : qValue);
                    //siteSum += piValue * q0[iPair];
                }
            }
            denom += siteSum;
        }
        //System.out.printf("denom = %s\n", denom);
        buffQ.close();
        buffPi.close();

        buffQ = new BufferedReader(new FileReader(dir + "QS.txt"));
        buffPi = new BufferedReader(new FileReader(dir + "PiS.txt"));

        double propsig = 0D;
        double propnonsig = 0D;
        int lines = 0;
        for (int i = 0; i < sites; i++) {
            lines++;
            String lineS = buffS.readLine();
            String lineQ = buffQ.readLine();
            String linePi = buffPi.readLine();
            String[] stringS = lineS.split("\\s");
            String[] stringQ = lineQ.split("\\s");
            String[] stringPi = linePi.split("\\s");

            for (int iPair = 0; iPair < 4096; iPair++) {
                int iCodon = (iPair - iPair%64)/64;
                char iRes = aaCode[iCodon];
                char jRes = aaCode[iPair%64];
                double piValue = Double.parseDouble(stringPi[iCodon]);
                double qValue = Double.parseDouble(stringQ[iPair]);
                //if ( (iRes != jRes) ) {
                //if ( (iCodon != (iPair%64)) ) {
                boolean isCodonChange = iCodon != (iPair % 64);
                boolean isNonSynChange = iRes != jRes;

                if (doNonSyn ? isNonSynChange : isCodonChange) {
                    double deltaS = Double.parseDouble(stringS[iPair]);
                    deltaS = Math.max(deltaS, -10);
                    deltaS = Math.min(deltaS, 10);
                    double val = deltaS - low;
                    int iBin = (int) (numbins * (val / (hi - low)));
                    try {
                        if (i < 734) { propnonsig += (piValue * (doMuts ? q0[iPair] : qValue)) / denom; } else { propsig += (piValue * (doMuts ? q0[iPair] : qValue)) / denom;}
                        bins[iBin] += (piValue * (doMuts ? q0[iPair] : qValue)) / denom;
                        //kbins[lines - 1][iBin] += (piValue * (doMuts ? q0[iPair] : qValue)) / denom;

                        if (i < 734) {
                            nonsigbins[iBin] += (piValue * (doMuts ? q0[iPair] : qValue)) / denom;
                        } else {
                            sigbins[iBin] += (piValue * (doMuts ? q0[iPair] : qValue)) / denom;
                        }
                        //bins[iBin] += (piValue*q0[iPair]) / denom;
                    } catch (Exception e) {
                        /*System.out.printf("deltaS = %s; low = %s; val = %s; iBin = %s; bins.length = %s\n",
                                deltaS, low, val, iBin, bins.length);*/
                    }

                }
            }
        }

        for (int i = 0 ; i < bins.length; i++) {
            System.out.printf("%s\t%s\t%s\n", bins[i], nonsigbins[i], sigbins[i]);
            // System.out.printf("%s\n", bins[i]);
        }


        System.exit(1);
        // bootstrap
        // Construct various bootstraps
        double[][] bStrapHisto = new double[1000][numbins];

        Random random = new Random();
        for (int iBoot = 0; iBoot < 1000; iBoot++) {
            for (int iLocation = 0; iLocation < sites; iLocation++) {
                int iPoint = random.nextInt(sites);
                for (int iS = 0; iS < numbins; iS++) {
                    bStrapHisto[iBoot][iS] += kbins[iPoint][iS];
                }
            }
        }

        double[] avg = new double[numbins];
        double[] avg2 = new double[numbins];

        // Normalise all distributions
        double summ = 0.0;
        for (int iBoot = 0; iBoot < 1000; iBoot++) {
            summ = 0.0;
            for (int iS = 0; iS < numbins; iS++) {
                summ += 0.25 * bStrapHisto[iBoot][iS];
            }
            for (int iS = 0; iS < numbins; iS++) {
                bStrapHisto[iBoot][iS]/= summ;
                avg[iS] += bStrapHisto[iBoot][iS];
                avg2[iS] += bStrapHisto[iBoot][iS] * bStrapHisto[iBoot][iS];
            }
        }
        summ = 0.0;
        for (int iS = 0; iS < numbins; iS++) {
            summ += 0.25 * bins[iS];
            avg[iS] /= 1000.0;
            avg2[iS] /= 1000.0;
            avg2[iS] = Math.sqrt(avg2[iS]-avg[iS]*avg[iS]);
        }

        double[][] results = new double[5][numbins];

        for (int iS = 0; iS < numbins; iS++) {
            bins[iS]/= summ;
            results[0][iS] = iS;
            results[1][iS] = bins[iS];
            results[2][iS] = avg[iS];
            results[3][iS] = avg2[iS];

            System.out.printf("%s,%s,%s,%s\n", results[0][iS], results[1][iS], results[2][iS], results[3][iS]);
        }





        //System.out.printf("lines = %s\n", lines);

        //System.out.printf("prop non-sig: %s; prop sig: %s\n", propnonsig, propsig);

    }


}
