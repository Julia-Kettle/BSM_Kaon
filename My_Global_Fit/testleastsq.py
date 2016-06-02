from scipy.optimize import leastsq
import numpy as np


xdata = np.array(range(2000))
ydata = np.array(range(0,4000,2))


def funcquad(p,x):
    return p[2]*x*x+p[1]*x+p[0]

def diff(p,y,x):
    return  y - funcquad(p,x)

pinit = np.array([1.0,1.0,1.0])

print np.shape(pinit)
print np.shape(xdata)


p,cov,infodict,mesg,ier=leastsq(diff, pinit,args=(ydata,xdata),full_output=True)
#print pfinal

print p
print cov
print infodict
print mesg
print ier

#print success
    