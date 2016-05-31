import nonlinearfit as nlf
import numpy as np
import random
import scipy.optimize

def quadratic(x,p):
    return p[0] + p[1]*x + p[2]*x*x

def cosinefunc(x,p):
    return p[0]*np.cos(p[1]*x+p[2])

x = range(100)
y = []
for i in range(len(x)):
    #y.append(2.02*np.exp(0.49*i))
    y.append(4*np.exp(-3*x[i]))



yerr = []
for i in range(len(y)):
    yerr.append(0.01*y[i]) # yerr[i]*y[i]/100.0

p0 = [3,-2]

p = [4.1,-3.001]

chisq = nlf.chi2(p,y,x,yerr,nlf.exponential)

print chisq/100
minchi = scipy.optimize.fmin(nlf.chi2,p0,args=(y,x,yerr,nlf.exponential,))
print minchi

print scipy.optimize.leastsq(nlf.difference,p0, args = (y,x, yerr, nlf.exponential),full_output=False)

#print p, cov, infodict, mesg, ier



#print minimize(quadratic,1,args=([-5,2,1],))

#print minimize(quadratic,p0,args=(y,x,))



