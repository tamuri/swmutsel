package swmutsel.runner.distributed;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.collect.Lists;
import com.google.common.collect.Table;
import com.googlecode.hessianserver.HessianServer;
import com.googlecode.hessianserver.HessianServiceDefinition;
import pal.alignment.Alignment;
import swmutsel.ArgumentsProcessor;
import swmutsel.options.AlignmentConverter;
import swmutsel.options.GeneticCodeConverter;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.GeneticCode;

import java.net.InetAddress;
import java.util.Enumeration;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.LogManager;

/**
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 30/10/2013 11:34
 */
public class Slave {

    public static final String PATH = "/swmutsel";
    public static HessianServer HESSIAN_SERVER;

    public Slave(String[] args) {
        ServerOptions options = new ServerOptions();
        new JCommander(options).parse(args);

        Table<String, Integer, Byte> allSites = new ArgumentsProcessor().getCleanedSitesTable(options.alignment);

        SlaveImpl service = new SlaveImpl(allSites, options.threads);

        HessianServiceDefinition hessianService = new HessianServiceDefinition(PATH, service, SlaveAPI.class);

        int port = CoreUtils.findFreePort();
        HESSIAN_SERVER = new HessianServer(port, hessianService);

        disableLoggers();

        List<String> addresses = getHostAddress();

        CoreUtils.msg("Host address(es): %s\n", addresses);
        CoreUtils.msg("Started slave (path '%s') on port %s\n", hessianService.getPath(), port);
        CoreUtils.msg("Use -hosts %s:%s\n", addresses.get(0), port);

        HESSIAN_SERVER.start();
    }

    public static void shutdown() {
        HESSIAN_SERVER.stop();
    }

    public static void main(String[] args) {
        new Slave(args);
    }

    private void disableLoggers() {
        for (Enumeration<String> name = LogManager.getLogManager().getLoggerNames(); name.hasMoreElements(); ) {
            LogManager.getLogManager().getLogger(name.nextElement()).setLevel(Level.OFF);
        }
    }

    private static List<String> getHostAddress() {
        List<String> addresses = Lists.newArrayList();

        try {
            InetAddress in = InetAddress.getLocalHost();
            InetAddress[] all = InetAddress.getAllByName(in.getHostName());

            for (InetAddress ia : all) {
                addresses.add(ia.getHostAddress());
            }

        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        return addresses;

    }
/*


    private static void writeHostnameFile(int port) throws Exception {
        InetAddress in = InetAddress.getLocalHost();
        InetAddress[] all = InetAddress.getAllByName(in.getHostName());

        FileWriter writer = new FileWriter("host_" + in.getHostName() + "_" + new Random().nextInt(100000) +".txt");
        BufferedWriter out = new BufferedWriter(writer);

        for (InetAddress anAll : all) {
            out.write(anAll + ":" + port + "\n");
        }

        out.close();
    }
*/

    private class ServerOptions {
        @Parameter(names = {"-s", "-sequences"}, converter = AlignmentConverter.class, required = true)
        public Alignment alignment;

        @Parameter(names = "-threads")
        public int threads = 1;

        @Parameter(names = {"-gc", "-geneticcode"}, converter = GeneticCodeConverter.class, required = true)
        public GeneticCode geneticCode;
    }
}
