import numpy as np
import random
from scipy.optimize import minimize

def chi2(p,ydata,xdata,yerr,f):
    #As long as user defines a function they pass this should for any function
    #hopefully
    chisq = 0
    for i in range(len(ydata)):
        chisq = chisq + pow((1/(yerr[i]))*(ydata[i] - f(xdata[i],p)),2)
        #where p is a list or vector of parameters
      
    return chisq

def exponential(x,p):
    return p[0]*np.exp(p[1]*x)

def log(x,p):
    return p[0]*np.log(p[1]*x)

def difference(p,ydata,xdata,yerr,f):
    #As long as user defines a function they pass this should for any function
    #hopefully
    difference = []
    for i in range(len(ydata)):
        difference.append((1/(yerr[i]))*(ydata[i] - f(xdata[i],p)))
        #where p is a list or vector of parameters
      
    return difference


        
