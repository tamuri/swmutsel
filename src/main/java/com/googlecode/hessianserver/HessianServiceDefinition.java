package com.googlecode.hessianserver;

import java.util.Arrays;

/**
 * Defines a service for {@link HessianServer} to make remotable
 * <p>
 * This is a combination of the path the client will provide on the url when
 * calling the HessianServer, the api to use (which will be an interface)
 * and finally the implementation (which implements the given interface).
 *
 * @author rayvanderborght
 */
public class HessianServiceDefinition
{
    private String path;
    public String getPath() { return this.path; }

    private Object implementation;
    public Object getImplementation() { return this.implementation; }

    private Class<?> api;
    public Class<?> getApi() { return this.api; }

    /**
     * Defines a service to be exposed by {@link HessianServer}
     *
     * @param path The path associated with this service
     * @param implementation The implementation of the api interface
     * @param api The interface implemented by the implementation object
     * @throws NullPointerException If any of the arguments provided are null
     * @throws IllegalArgumentException If the implementation object does not implement the given api
     */
    public HessianServiceDefinition(String path, Object implementation, Class<?> api)
    {
        if (path == null || implementation == null || api == null)
            throw new NullPointerException("Null input is not allowed");

        Class<?>[] interfaces = implementation.getClass().getInterfaces();
        if (interfaces == null || !Arrays.asList(interfaces).contains(api))
        {
            throw new IllegalArgumentException(
                    "The implementation class: " + implementation.getClass().getCanonicalName() +
                            " does not implement the given api: " + api.getCanonicalName());
        }
        this.path = path;
        this.implementation = implementation;
        this.api = api;
    }
}
