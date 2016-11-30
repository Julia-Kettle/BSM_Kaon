from scipy.optimize import leastsq
from scipy.optimize import minimize
import numpy as np
class Fit:

    def __init__(self):
        pass
    def residuals(self,params,y,erry,x1,x2,coeff):
        #returns difference between data and model over err
        res =  ((y - self.globalfunction(params,x1,x2,coeff))/erry)
        return res

    def globalfit(self,ydata,yerr,x1data,x2data,par_init,coeff):
        #do global fit through leastsq, and calculate chi2
        params, cov,infodict,mesg,ier = leastsq(self.residuals,par_init,args=(ydata,yerr,x1data,x2data,coeff),full_output=True,maxfev=100000000)
        #print minimize(self.chi2,pinit,args=(ydata,yerr,x1data,x2data),method='Nelder-Mead') #alternative minimisaiton method
        chisq, chidof = self.chi2(params,ydata,yerr,x1data,x2data,coeff)
        return params, cov, chisq, chidof

    def globalfunction(self,params,msq,asq,coeff):
        #we define the function, either the linear ansatz or chiral with the extra log term
        func = params[0] + params[1]*msq + params[2]*asq + params[0]*coeff*np.log(msq)*msq
        return func

    def chi2(self,params,ydata,yerr,xdata,x2data,coeff):
        #calculate and retun the chi2 and chi2/dof
        chisq = 0
        for i in range(len(ydata)):
            chisq = chisq +pow(self.residuals(params,ydata[i],yerr[i],xdata[i],x2data[i],coeff),2)
            chidof = chisq/len(ydata)
        return chisq, chidof
