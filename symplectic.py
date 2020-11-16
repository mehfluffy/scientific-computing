def odeNDLeapfrog(x0, v0, start, end, h, spacing, dvdt, *args):
    """
    Finds x_i(t), v_i(t) for i = 1, 2, ..., ndim
    from a 2nd order ODE using Leapfrog.
    All quantities must have compatible units.
        args: additional arguments for dvdt()
        x0, array of initial positions in kpc
        v0: array of initial velocities in kpc/Gyr
        start, end: interval [t_a, t_b] of integration in Gyr
        h: timestep between integration points, (t_b - t_a) must be
            divisible by h.
        spacing: time intervals at which to store synchronized x and v,
            must be multiple of h.
        dvdt: 2nd derivative of x(t), callable function of vector x.
        args: additional arguments for dvdt (unpacked).
    Returns arrays of x and v values, at output intervals
    """
    N = int((end-start)/h)
    every = spacing/h
    # store initial values
    x = np.array(x0)
    v = np.array(v0)
    # current value arrays
    x_val = x
    v_val = v
    # jumpstart
    v_val = v_val + h/2 * dvdt(x_val, *args)
    # remaining values
    j = 1
    while j <= N:
        x_val = x_val + h * v_val #/ kpc2km   # velocity becomes kpc/s
        if j%every != 0:
            v_val = v_val + h * dvdt(x_val, *args)
        else:
            # resync and store
            x = np.vstack((x, x_val))
            v_store = v_val + h/2 * dvdt(x_val, *args)
            v = np.vstack((v, v_store))
            v_val = v_val + h * dvdt(x_val, *args)
        j += 1
    return x.transpose(), v.transpose()
