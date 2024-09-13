import numpy as np
from scipy.optimize import curve_fit

def _linfit(x, a, b):
    return a + b*x

def _potfit(x, a, b):
    return a*(x**b)

def _invfit(x, a, b):
    return a/x + b

def linfit(x, y):
    return curve_fit(_linfit, x, y)

def potfit(x, y):
    return curve_fit(_potfit, x, y)

def logfit(x, y):
    (a, b), pcov = curve_fit(_linfit, np.log(x), np.log(y))

    return a, b

def invfit(x, y):    
    (a, b), pcov = curve_fit(_invfit, x, y)

    return a, b