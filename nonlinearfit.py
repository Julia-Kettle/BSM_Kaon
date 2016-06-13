import numpy as np
import random
from scipy.optimize import leastsq
from scipy import stats
from plotting import *

def chi2(p,ydata,yerr,xdata,x2data,f):
    #As long as user defines a function they pass this should for any function
    #hopefully
    chisq = 0
    for i in range(len(ydata)):
        #chisq = chisq + pow((1/(yerr[i]))*(ydata[i] - f(xdata[i],p)),2)
        #where p is a list or vector of parameters
        chisq = chisq +pow(residuals(p,ydata[i],xdata[i],x2data[i]),2)
      
    return chisq

def exponential(x,p):
    return p[0]*np.exp(p[1]*x)

def log(x,p):
    return p[0]*np.log(p[1]*x)

def difference(p,ydata,yerr,xdata,f):
    #As long as user defines a function they pass this should for any function
    #hopefully
    difference = []
    for i in range(len(ydata)):
        difference.append((1/(yerr[i]))*(ydata[i] - f(xdata[i],p)))
        #where p is a list or vector of parameters
      
    return difference


def residuals(p,y,x1,x2):
    return (y - globalfunction(p,x1,x2))


def globalfit(ydata,yerr,x1data,x2data,pinit):
    ydata = np.array(ydata)
    x1data=np.array(x1data)
    x2data=np.array(x2data)
    p, cov,infodict,mesg,ier = leastsq(residuals,pinit,args=(ydata,x1data,x2data),full_output=True)

    chisq = chi2(p,ydata,yerr,x1data,x2data,globalfunction)
    
    return p, cov, chisq
    
def globalfunction(p,asq,msq):
    return p[0] + p[1]*msq + p[2]*asq


def bootglobalfit(y1,y2,a1,a2,m1,m2,nboot,p0,mp,xlabel,ylabel,name,y3=[],a3=[],m3=[]):
    #have y1 for different strange masses - also have
    #data passed will have form [ms,nboot+1]
    err1 = np.std(y1[:][0:nboot],axis=1)
    err2 = np.std(y2[:][0:nboot],axis=1)
    yerr =  np.hstack((err1,err2))
    asq = np.zeros([len(m1)+len(m2)])
    asq[0:len(m1)] = a1
    for i in range(len(m2)):
	asq[len(m1)+i] = a2
    #msq=np.hstack((m1,m2))
    bp=[]
    byp = []
    bcov = []
    bchisq = []
    if type(y3)==np.ndarray and type(a3)==np.ndarray and type(m3)==np.ndarray:
	np.hstack((asq,a3))
	err3 = np.std(y3[:][0:nboot],axis=1)
    for j in range(nboot+1):
	y = []
	msq = []
	for i in range(len(y1)):
	    y.append(y1[i,j])
	    msq.append(m1[i,j])
	for i in range(len(y2)):
	    y.append(y2[i,j])
	    msq.append(m2[i,j])
	if type(y3)==np.ndarray and type(a3)==np.ndarray and type(m3)==np.ndarray:
	    np.hstack((y,y3[:,j]))
	    np.hstack((msq,m3[:,j]))
	y = np.array(y)
	msq = np.array(msq)

	p, cov, chisq = globalfit(y,yerr,msq,asq,p0)
	bp.append(p)
	bcov.append(cov)
	bchisq.append(chisq)

	yp = globalfunction(p,0,mp)
	byp.append(yp)

    perr = np.std(bp[0:nboot],axis=0)
    pcent = bp[-1]
    covcent = bcov[-1]
    if type(y3)==np.ndarray and type(a3)==np.ndarray and type(m3)==np.ndarray:
	plot_Global(pcent,globalfunction,y1,y2,err1,err2,m1,m2,mp,a1,a2,xlabel,ylabel,name,y3,err3,m3,a3)
    else:
	plot_Global(pcent,globalfunction,y1,y2,err1,err2,m1,m2,mp,a1,a2,xlabel,ylabel,name)
    
     
    
    return pcent, covcent, bp, bcov, bchisq, byp

def newbootglobalfit(y,a,m,nboot,p0,mp):
    #are passed equal length lists of data pnts, a, m.
    err = np.std(y[:][0:nboot],axis=1)
    bp = np.array([nboot+1])
    bcov = np.array([nboot+1])
    bchisq = np.array([nboot+1])

    for i in range(nboot+1):
        p, cov, chisq = globalfit(y,yerr,msq,asq,p0)
        bp[i]=p
        bcov[i]=cov
        bchisq[i]=chisq

    perr = np.std(bp[0:nboot],axis=0)
    pC = bp[-1]
    covC = bcov[-1]
    
    ynew=[0,0]
    y1=[]
    y2=[]
    for i in range(len(y)):
        if i < 3:
            y1.append(y[i])
        elif i < 5:
            y2.append(y[i])
        else:
            y3.append(y[i])
    
    new_plot_Global(pC,globalfunction,y1,y2,y3,err,m,a,p0)


        
