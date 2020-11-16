def binary_to_decimal(b):
    
    neg = False
    if b[0] == "-":                     # if number is negative
        b = b[1:len(b)]                   # slice the sign away
        neg = True
    
    df = 0
    dw = 0
    
    w = b
    for i in range(0,len(b)):
        if b[i] == "." :                # if number is fraction
            f = b[i+1:len(b)]             # fractional part
            w = b[0:i]                    # whole part
            for j in range(0,len(f)):
                df = df + int(f[j])* 2**(-(j+1))
    for k in range(0,len(w)):
        dw = dw + int(w[k])* 2**(len(w)-k-1)
    
    d = dw+df
    if neg:
        d = -d
    return d
