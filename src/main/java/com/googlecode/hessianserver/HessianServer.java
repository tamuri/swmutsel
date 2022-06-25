/*
 * $Id$
 * $URL$
 */
package com.googlecode.hessianserver;

import com.caucho.hessian.io.SerializerFactory;
import com.caucho.hessian.server.HessianSkeleton;
import com.caucho.services.server.ServiceContext;
import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;

import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

/**
 * Listens for requests from hessian clients on a given port and invokes the
 * proper methods on the corresponding api implementation.
 * <p>
 * The available services are differentiated by the path the client sends in
 * the url and should match one supplied in the {@link HessianServiceDefinition}
 * passed into the constructor.
 *
 * @author rayvanderborght
 */
public class HessianServer
{
    private static final Log log = LogFactory.getLog(HessianServer.class);

    private static int DEFAULT_PORT = 4444;

    private int port;
    private Lock serverStateTransitionLock = new ReentrantLock();
    private ServerState currentState;
    private ServerSocket listener;

    private Map<String, HessianServiceDefinition> services = new HashMap<String, HessianServiceDefinition>();

    /** */
    public HessianServer(HessianServiceDefinition... services)
    {
        this(DEFAULT_PORT, services);
    }

    /** */
    public HessianServer(int port, HessianServiceDefinition... services)
    {
        if (services == null)
            throw new NullPointerException("Must provide services to HessianServer");

        this.currentState = ServerState.STOPPED;
        this.port = port;
        for (HessianServiceDefinition service : services)
        {
            this.services.put(service.getPath(), service);
        }
    }

    /**
     * Start the server
     */
    public void start()
    {
        this.serverStateTransitionLock.lock();
        try
        {
            if (this.currentState != ServerState.STOPPED)
            {
                log.warn("Aborting server start, HessianServer is already " + this.currentState);
                return;
            }

            this.currentState = ServerState.STARTING;
            log.info("HessianServer is " + this.currentState);

            Thread t = new Thread(new Listener());
            t.start();

            this.currentState = ServerState.RUNNING;
            log.info("HessianServer is " + this.currentState);
        }
        finally
        {
            this.serverStateTransitionLock.unlock();
        }
    }

    /**
     * Stop the server
     */
    public void stop()
    {
        this.serverStateTransitionLock.lock();
        try
        {
            this.currentState = ServerState.STOPPING;
            log.info("HessianServer is " + this.currentState);

            this.listener.close();

            this.currentState = ServerState.STOPPED;
            log.info("HessianServer is " + this.currentState);
        }
        catch (IOException e)
        {
            if (this.currentState == ServerState.STOPPED)
                log.warn("Unexpected error during server stop, HessianServer is already " + this.currentState, e);
        }
        finally
        {
            this.serverStateTransitionLock.unlock();
        }
    }

    /**
     * Open a server socket listening on the defined port.  Accept and process
     * connections until the listener is closed.
     *
     * @author rayvanderborght
     */
    private class Listener implements Runnable
    {
        @Override
        public void run()
        {
            if (log.isInfoEnabled())
                log.info("Starting HessianServer on port " + HessianServer.this.port);

            try
            {
                HessianServer.this.listener = new ServerSocket(HessianServer.this.port);

                while (true)
                {
                    Socket server = HessianServer.this.listener.accept();
                    Worker worker = new Worker(server, HessianServer.this.services);
                    Thread t = new Thread(worker);
                    t.start();
                }
            }
            catch (IOException e)
            {
                if (HessianServer.this.currentState == ServerState.RUNNING)
                    log.error("IOException on socket listen", e);
            }
        }
    }

    /**
     * Worker thread for invoking an api call on a given implementation object using hessian.
     *
     * @author rayvanderborght
     */
    private static class Worker implements Runnable
    {
        /*
         * We reuse the serializer factory, which is safe in recent versions of hessian:
         * http://stackoverflow.com/questions/1474038/is-the-hessian-class-serializerfactory-thread-safe
         */
        private static final SerializerFactory serializerFactory = new SerializerFactory();

        /** */
        private Socket server;

        /** */
        private Map<String, HessianServiceDefinition> services;

        /** */
        private Worker(Socket server, Map<String, HessianServiceDefinition> services)
        {
            this.server = server;
            this.services = services;
        }

        /** */
        public void run()
        {
            try
            {
                InputStream in = this.server.getInputStream();
                PrintStream out = new PrintStream(this.server.getOutputStream());

                String requestLine = readHttpRequestLine(in);
                String path = this.getRequestPath(requestLine);

                HessianServiceDefinition service = this.services.get(path);

                if (!requestLine.startsWith("POST") || "".equals(path) || service == null)
                {
                    out.print("HTTP/1.1 500 Internal Server Error\r\n");
                    out.print("Content-Type: text/html; charset=UTF-8\r\n");
                    out.print("Content-Length: 30\r\n");
                    out.print("\r\n");
                    out.print("<h1>Hessian Requires POST</h1>");
                    return;
                }

                while (!"".equals(readHttpHeaderLine(in)))
                {
                    // skip straight to the content; we don't care about the headers
                }

                // older hessian versions made me skip the next line of input
                // for some unknown reason.  they must've fixed a bug or
                // something so thankfully we don't have to do this anymore...
                // NOTE: if you're not using a recent hessian jar, you may need
                // put this back in:
                // readHttpHeaderLine(in);

                out.write("HTTP/1.1 200\r\n".getBytes());
                out.write("Content-Type: application/x-hessian\r\n\r\n".getBytes());

                // AUT: begin() has changed in Hessian-4.0.37 - fix
                ServiceContext.begin(null, null, path, null);

                HessianSkeleton hessianObjectSkeleton = new HessianSkeleton(service.getImplementation(), service.getApi());
                hessianObjectSkeleton.invoke(in, out, serializerFactory);

                out.write("\r\n".getBytes());

                out.close();
                this.server.close();
            }
            catch (Exception e)
            {
                log.error("Problem handling remote api call: ", e);
            }
        }

        /* */
        private String getRequestPath(String requestLine)
        {
            int beginIndex = requestLine.indexOf('/');
            int endIndex = requestLine.indexOf(' ', beginIndex);

            return requestLine.substring(beginIndex, endIndex);
        }

        /* */
        private static String readHttpRequestLine(InputStream in) throws IOException
        {
            return readHttpHeaderLine(in);
        }

        /* */
        private static String readHttpHeaderLine(InputStream in) throws IOException
        {
            byte[] buffer = new byte[8 * 1024];	// no defined limit, but apache's is 8k so let's go with that
            int next;
            int count = 0;
            while (true)
            {
                next = in.read();
                if (next < 0 || next == '\n')
                    break;

                if (next != '\r')
                    buffer[count++] = (byte) next;

                if (count >= buffer.length)
                    throw new IOException("HTTP Header too long");
            }
            return new String(buffer, 0, count);
        }
    }

    /** */
    private enum ServerState
    {
        STARTING,
        RUNNING,
        STOPPING,
        STOPPED;
    }
}
