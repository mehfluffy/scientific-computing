#!/usr/bin/python
"""
A script to find the root of a 1D function using bisection.
"""

from numpy import *

def f(x, func):
    return(eval(func))

def bisection(a, b, func, tol=1e-4, N=100):
    """
    Only takes functions with roots, and only output one of the roots.
    Starting conditions: f is continuous, f(a) and f(b) are opposite sign.
    """
    i = 1
    if f(a, func) == 0:
        print("f({}) = 0".format(a))
        return a
    if f(b, func) == 0:
        print("f({}) = 0".format(b))
        return b
    while i <= N:
        p = (a+b)/2
        if f(p,func) == 0 or abs(b-a) < tol:
            print("f(x) = 0 at x = {}".format(p))
            print("found after {} iterations.".format(i))
            return p
        if f(a, func) * f(p, func) > 0:
            a = p
        else: b = p
        i+=1
    print("Method failed after {} iterations".format(N))
    return None


tol = 1e-8
N = 1000
while True:
    func = input("\nEnter f(x): ")
    a,b = eval(input("Enter x interval as a,b: "))
    x = linspace(a, b, N*10)
    y = f(x, func)
    if y.max() * y.min() > 0 or isinstance(y, ndarray) == False:   # false means y is constant
        print("f(x) does not have a root (in the given interval).")
    elif f(a, func) * f(b, func) > 0:   # this does not imply lack of root, but only user error
        print("f(a) and f(b) must be of opposite sign.")
    else:
        root = bisection(a, b, func, tol, N)
        break
