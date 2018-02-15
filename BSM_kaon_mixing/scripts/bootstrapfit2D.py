import numpy as np 
import matplotlib.pyplot as plt
import fit
import globaldefs
from datasets import *

# -----------   Start of class  -----------------
class BootstrapFit2D(object):
    
    ##########################################################################################
    #class to perform non linear fit with bootstraps for function of form y=f(x1,x2)
    # loops through the bootstraps and performs nonlinear fit, returning nboots sets of parameters
    def __init__(self,datasets,fitfunction,p_init,coeff=None):
        
        #initialise object with list of datasets
        self.dslist=datasets

        #set number of points in datasets + no of bootstraps
        self.ndp = datasets[0].no_points
        self.nboots=len(self.dslist)-1 # as last value should be central value
        
        #fitfunction of form y=f(x1,x2,coeff)
        self.fitfunction=fitfunction
        #allows for coeff to be passed to fitfunction
        self.coeff=coeff    

        #number of parameters in fit taken from initial par list provided
        self.nparams=len(p_init)
        self.p_init=p_init
        #set up empty parameter results and chi/dof arrays to store results
        self.bparams=np.zeros([self.nboots+1,self.nparams])
        self.bchidof=np.zeros([self.nboots+1])

        #weighting factor default set to 1 for all points.
        self.wf=np.ones(self.ndp)
        print("fit object created for data set - \n"+str(self.dslist[-1]))

    ##########################################################################################
    #performs the non-linear fit
    def dofit(self):
        #set up fit results arrays
        bcovmat=np.zeros([self.nboots+1,self.nparams,self.nparams])
        byphys=np.zeros([self.nboots+1])
        bchisq=np.zeros([self.nboots+1])
        bqvals=np.zeros([self.nboots+1])

        #loop through each bootstrap
        for i in range(self.nboots+1):
            #perform fit
            #print(self.dslist[i].y,self.dslist[i].ey,self.dslist[i].x1,self.dslist[i].x2)
            self.bparams[i],bcovmat[i],bchisq[i],self.bchidof[i],bqvals[i] = fit.nonlinearfit(self.fitfunction,self.p_init,
            self.dslist[i].y,self.dslist[i].ey,self.dslist[i].x1,self.dslist[i].x2,self.coeff,wf=self.wf)

    ##########################################################################################
    #takes values of x1 and x2 to return a value of y given by the fitted function for each bootstrap of the parameters
    def bootphysvalue(self,x1_phys,x2_phys):
        byphys=np.zeros(self.nboots+1)
        for i in range(self.nboots+1):
            byphys[i] = self.fitfunction(self.bparams[i],x1_phys,x2_phys,self.coeff)
        eyphys=np.std(byphys)
        return byphys

    def appendPhysVal(self,x1_phys,x2_phys):
        byphys=self.bootphysvalue(x1_phys,x2_phys)
        eyphys=np.std(byphys[:-1])
        group=self.dslist[-1].group[-1]+1
        name='phys'

        for i in range(self.nboots+1):
            dpPhys=DataPoint(group,name,x1_phys,x2_phys,byphys,eyphys)
            self.dslist[i].append(dpPhys)

    ##########################################################################################
    #returns a value of y for given x1,x2 from fitfunction with the central value of the parameters from the fit
    def centralphysvalue(self,x1_phys,x2_phys):
        byphys = self.bootphysvalue(x1_phys,x2_phys)
        return byphys[-1],np.std(byphys[:-1])

    ##########################################################################################
    #returns arrays [xa,xb,xc,...,xN] and [ya,yb,yc,...,yN] which we call "fit lines" where other variable held constant
    # intended for use in plotting
    # constval gives the values which x1/2 held constant to, xrange range over which variable x is sampled
    # axis sets which x is varied
    def return_fit_line(self,constval,xrange,axis=0):
            #set up heavily populated list of x variables 
            xlist=np.linspace(xrange[0],xrange[1],1000)
            #user defines along which axis fit line goes, 0=x1,1=x2
            if axis == 0:
                #x1 is variable x2 held constant
                x1=xlist
                x2=constval
            elif axis ==1:
                x1=constval
                x2=xlist
            else:
                print("Error: axis must be 0 or 1")
            #returns the list and the y=f(xlist,const,coeff)
            return xlist, self.fitfunction(self.bparams[-1],x1,x2,self.coeff)

    ##########################################################################################
    # plots fit lines to a provided axis ax
    #user choose which axis (0=x1,1=x2) is varied. other is taken as datapoints values
    def plot_fitlines(self,ax,axis):
        
        #set up xrange
        #xrange=[self.dslist[-1].min(axis),self.dslist[-1].max(axis)]
        #span=xrange[1]-xrange[0]
        #xrange[0]-=0.1*span
        #xrange[1]+=0.1*span

        xrange=ax.get_xlim()


        #constant is x2 if axis = 0 and x1 is variable and vice vers
        xconst = self.dslist[-1].x1 if axis==1 else self.dslist[-1].x2

        #store the values of x1(2) which are taken from data points values
        constvalues=[]
        for i in range(len(xconst)):
            #check we're not replicating fitline if 2 datapoints share an x value
            if xconst[i] not in constvalues:
                #store value for above check
                constvalues.append(xconst[i])
                x,y=self.return_fit_line(xconst[i],xrange,axis=axis)
                ax.plot(x,y,self.dslist[-1].linestyle[i],color=self.dslist[-1].colour[i])#,linestyle=self.dslist[-1].linestyle[i])
        ax.set_xlim(xrange)
        return ax
        #get fit line
        
        

    ##########################################################################################

