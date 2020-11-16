def pearsoncoef(x,y):
    xm = x.mean(); ym = y.mean()
    Sxy = ((x-xm)*(y-ym)).sum()
    Sxx = ((x-xm)*(x-xm)).sum()
    Syy = ((y-ym)*(y-ym)).sum()
    r = Sxy/(Sxx**0.5*Syy**0.5)
    return r
    
def weightedmoments(x,w):
    gam = w.sum()
    XY = (x*w).sum()
    mean = XY/gam
    var = (1/gam)*(w*(x-mean)**2).sum()
    std = np.sqrt(var)
    skew = (1/gam)*(w*(x-mean)**3).sum()/std**3
    kurt = (1/gam)*(w*(x-mean)**4).sum()/std**4
    return mean, var, std, skew, kurt
    
def denoise(s, n=10):
    """
    Denoises 1D signal s with the procedure - 
    apply FFT, filter by taking only n terms, apply inverse FFT.
    Returns denoised signal and complex terms
    (use attributes .real and .imag).
    """
    C = np.fft.fft(s)
    N = len(s)
    start = n; end = N-n
    C[start:end] = 0.0
    D = np.fft.ifft(C)
    return D, C/N
