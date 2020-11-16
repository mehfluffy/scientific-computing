def trapint(f, a, b, N):
    """
    Definite integration using the trapezoidal rule.
    f: callable function to be integrated along 1st axis
    a,b: integration from a to b
    N: number of steps
    """
    h = (b-a)/N
    x = np.linspace(a, b, N)
    fs = np.array(f(x))
    for i in [0, N-1]:
        fs[i] = fs[i]/2
    integral = h*fs.sum()
    return integral


def odeRK4(y0, dydx, start, end, N):
    """
    Finds y(x) from y(0) and dy/dx using Runge-Kutta
    y0: initial value of y(x)
    dydx: derivative dy/dx, callable function of (x,y)
    start, end: interval of integration
    N: number of steps
    """
    h = (end-start)/N
    xvals = np.linspace(start, end, N)
    yvals = np.empty((N))
    yvals[0] = y0
    yn = y0
    for i in range(N-1):
        xn = xvals[i+1]
        f1 = dydx(xn, yn)
        f2 = dydx(xn+h/2, yn+h/2)
        f3 = dydx(xn+h/2, yn+(h/2)*f2)
        f4 = dydx(xn+h, yn+h*f3)
        yn = yn + (h/6)*f1 + (h/3)*f2 + (h/3)*f3 + (h/6)*f4
        yvals[i+1] = yn
    return yvals


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
