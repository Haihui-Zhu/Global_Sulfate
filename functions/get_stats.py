
import numpy as np
from scipy.stats import linregress

def get_stats(data1,data2):
    # remove nan
    x = data1.flatten()
    y = data2.flatten()
    valid = ~np.isnan(x) & ~np.isnan(y)  # ~ is the bitwise NOT operator, used here to invert the boolean array
    x = x[valid]
    y = y[valid]

    correlation = np.corrcoef(x, y)[0, 1]
    r2 = correlation**2

    # Perform linear regression
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    intercept = intercept + np.mean(x)*slope*(1-1/correlation); 
    slope = slope/correlation
    
    it = '-+'
    itv = abs(intercept)
    N = len(x)
    statstrs = f"y = {slope:.2f}x {it[int(intercept>0)]} {itv:.2f}\n $r^2$ = {r2:.2f}\nN = {N}"
    
    return statstrs, slope, intercept, N, x, y