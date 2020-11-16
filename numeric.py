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
