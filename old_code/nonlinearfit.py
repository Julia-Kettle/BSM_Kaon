import numpy as np
import random
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy import stats
from plotting import *


##################################################################################
def chi2(p,ydata,yerr,xdata,x2data):
    #As long as user defines a function they pass this should for any function
    #hopefully
    chisq = 0
    for i in range(len(ydata)):
        #chisq = chisq + pow((1/(yerr[i]))*(ydata[i] - f(xdata[i],p)),2)
        #where p is a list or vector of parameters
        chisq = chisq +pow(residuals(p,ydata[i],yerr[i], xdata[i],x2data[i]),2)
      
    return chisq
###################################################################################
def exponential(x,p):
    return p[0]*np.exp(p[1]*x)
###################################################################################
def log(x,p):
    return p[0]*np.log(p[1]*x)
###################################################################################
def difference(p,ydata,yerr,xdata,f):
    #As long as user defines a function they pass this should for any function
    #hopefully
    difference = []
    for i in range(len(ydata)):
        difference.append((1/(yerr[i]))*(ydata[i] - f(xdata[i],p)))
        #where p is a list or vector of parameters
      
    return difference
###################################################################################


###################################################################################
def residuals(p,y,ey,x1,x2):
    res =  ((y - globalfunction(p,x1,x2))/ey)
    return res
###################################################################################


###################################################################################
def globalfit(ydata,yerr,x1data,x2data,pinit):
    ydata = np.array(ydata)
    yerr = np.array(yerr)
    x1data=np.array(x1data)
    x2data=np.array(x2data)
    p, cov,infodict,mesg,ier = leastsq(residuals,pinit,args=(ydata,yerr,x1data,x2data),full_output=True,maxfev=10000)
    chisq = chi2(p,ydata,yerr,x1data,x2data) 
    return p, cov, chisq
###################################################################################
    
###################################################################################
def globalfunction(p,msq,asq):
    return p[0] + p[1]*msq + p[2]*asq
###################################################################################


###################################################################################
def bootglobalfit(y1,y2,a1,a2,m1,m2,nboot,p0,mp,xlabel,ylabel,name,y3=[],a3=[],m3=[]):
    #have y1 for different strange masses - also have
    #data passed will have form [ms,nboot+1]
    err1 = np.std(y1[:][0:nboot],axis=1)
    err2 = np.std(y2[:][0:nboot],axis=1)
    yerr =  np.hstack((err1,err2))
    #yerr = err1
    #set the array with the errors on y
    bp=[]
    byp = []
    bcov = []
    bchisq = []
    if type(y3)==np.ndarray and a3 and type(m3)==np.ndarray:
        err3 = np.std(y3[:][0:nboot],axis=1)
        yerr = np.hstack((yerr,err3[:]))
        #loop through boots and for each append all ys of that boot
    for j in range(nboot+1):
        y = []
        msq = []
        asq = []
        for i in range(len(y1)):
            y.append(y1[i,j])
            msq.append(m1[i,j])
            asq.append(1.0/(a1[j]*a1[j]))
        for i in range(len(y2)):
            y.append(y2[i,j])
            msq.append(m2[i,j])
            asq.append(1.0/(a2[j]*a2[j]))
        if type(y3)==np.ndarray and a3 and type(m3)==np.ndarray:
            y = np.hstack((y,y3[:,j]))
            msq = np.hstack((msq,m3[:,j]))
            asq.append(1.0/(a3[0][j]*a3[0][j]))
            asq.append(1.0/(a3[1][j]*a3[1][j]))
        y = np.array(y)
        msq = np.array(msq)
        p, cov, chisq = globalfit(y,yerr,msq,asq,p0)
        bp.append(p)
        bcov.append(cov)
        bchisq.append(chisq)

        yp = globalfunction(p,0,mp)
        byp.append(yp)
    eyp = np.std(byp)
    perr = np.std(bp[0:nboot],axis=0)
    pcent = bp[-1]
    #pcent = [0,0,0]
    #for iboot in range(500):
    #    for j in range(3):
    #        pcent[j] = pcent[j] + bp[-1][j]
    #for j in range(3):
    #    pcent[j] = pcent[j]/500
    covcent = bcov[-1]
    if type(y3)==np.ndarray and a3 and type(m3)==np.ndarray:
        plot_Global(pcent,globalfunction,y1,y2,err1,err2,eyp,m1,m2,mp,a1,a2,xlabel,ylabel,name,y3,err3,m3,a3)
    else:
        plot_Global(pcent,globalfunction,y1,y2,err1,err2,eyp,m1,m2,mp,a1,a2,xlabel,ylabel,name)
    return pcent, covcent, bp, bcov, bchisq, byp
###################################################################################

