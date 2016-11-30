import numpy as np
import random
import math
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy import stats
import matplotlib.pyplot as plt

class GlobalFit:
    def __init__(self,betas,data,m2f2,a2,val,ischeme,ikin,ibasis,dolog):
        #where ibeta is now a list [1,1,2,2,2,3,4] of the indices
        self.scheme_name = ['MOM','ms'][ischeme]
        self.basis_name = ['SUSY','Lattice'][ibasis]
        self.kin_name = ['E','qq','qg','gq','gg'][ikin]
        self.lattice_name = ['24','32','48','64']
        self.nchan = 5
        self.betas = betas
        self.data = data
        self.x1 = m2f2
        self.x2 = a2
        self.nboots = 500
        self.mp= pow(0.13957018/(4*math.pi*0.13),2)
        self.val=val
        self.dolog = dolog
        print "do log value is - ", self.dolog
        if self.val == "R":
            val_no = 0
        elif self.val == "B":
            val_no = 1
        self.coeff_set = [[[1.5,2.5,2.5,1.5,1.5],[1.5,1.5,1.5,2.5,2.5]],[[-0.5,0.5,0.5,-0.5,-0.5],[-0.5,-0.5,-0.5,0.5,0.5]]][val_no][ibasis]
    def SetUp(self):
        print "---Doing Global Fitting---"
        for ic in range(self.nchan):
            print "Scheme - ",self.scheme_name,"  Basis - ",self.basis_name,"  Kinematics - ",self.kin_name
            print "Doing ",self.val
            print "Channel - ",ic
            self.coeff = self.coeff_set[ic]
            bp=[]
            bchisq=[]
            bchidof=[]
            byp=[]
            for iboot in range(self.nboots+1):
                datapoints = np.zeros(len(self.betas))
                dataerror = np.zeros(len(self.betas))
                x1points = np.zeros(len(self.betas))
                x2points = np.zeros(len(self.betas))
                j=0
                for i in range(len(self.betas)):
                    if self.betas[i] == self.betas[i-1]:
                        j=j+1
                    else:
                        j=0
                    datapoints[i] = self.data[i][ic,0,iboot]
                    dataerror[i] = np.std(self.data[i][ic,0,:])
                    x1points[i] = self.x1[i][0,0,iboot]
                    x2points[i] = self.x2[i][iboot]
                #Here we need to carry out the global fi
                
                p, cov, chisq, chidof = self.globalfit(datapoints,dataerror,x1points,x2points,[1,1,1])
                bp.append(p)
                bchisq.append(chisq)
                bchidof.append(chidof)
                byp.append(self.globalfunction(p,self.mp,0))
            yperr = np.std(byp)
            pcent = bp[-1]
            print "---Fitting Results---"
            print "parameters: "
            print "p0 - ", pcent[0], " +/- ", np.std(bp[0])
            print "p1 - ", pcent[1], " +/- ", np.std(bp[1])
            print "p2 - ", pcent[2], " +/- ", np.std(bp[2])
            print "\n"
            print "chisq - ", bchisq[-1], " +/- ", np.std(bchisq)
            print "chisq/dof - ", bchidof[-1], "+/-", np.std(bchidof)
            name="./plots/" +self.val + str(ic+1) + "_"+ self.basis_name+"_"+self.scheme_name+"_"+self.kin_name
            #print "---Plotting for channel ", ic, "---"
            self.plot_Global(pcent,datapoints,dataerror,x1points,x2points,yperr,r"$m^2_{\pi}/4\pi f^2_{\pi}$",r"$"+self.val+"_{"+str(ic+1)+"}$",name)
            #We only call plotting after loop so last value (ie central values) plotted

    ###################################################################################
    def residuals(self,p,y,ey,x1,x2):
        res =  ((y - self.globalfunction(p,x1,x2))/ey)
        return res
    ###################################################################################


    ###################################################################################
    def residuals_chiral(self,p,y,ey,x1,x2):
        res =  ((y - self.globalfunction_chiral(p,x1,x2))/ey)
        return res
    ###################################################################################

    def function_to_fit(self,p,y,ey,x1,x2):
        val = 0
        for i in range(len(y)):
            val += pow(self.residuals(p,y[i],ey[i],x1[i],x2[i]),2)
        return val

    ###################################################################################
    def globalfit(self,ydata,yerr,x1data,x2data,pinit):
        ydata = np.array(ydata)
        yerr = np.array(yerr)
        x1data=np.array(x1data)
        x2data=np.array(x2data)
        p, cov,infodict,mesg,ier = leastsq(self.residuals,pinit,args=(ydata,yerr,x1data,x2data),full_output=True,maxfev=100000000)
        #print minimize(self.function_to_fit,pinit,args=(ydata,yerr,x1data,x2data),method='Nelder-Mead')
        diffs = infodict['fvec']
        sum1 = 0
        for i in range(len(diffs)):
            sum1 += diffs[i]*diffs[i]
        #print sum1/len(diffs)
        chisq, chidof = self.chi2(p,ydata,yerr,x1data,x2data)
        return p, cov, chisq, chidof
    ###################################################################################

    ###################################################################################
    def globalfunction(self,p,msq,asq):
        if self.dolog == 1:
            func = p[0] + p[1]*msq + p[2]*asq + p[0]*self.coeff*np.log(msq)*msq
        else:
            func = p[0] + p[1]*msq + p[2]*asq
        return func
    ###################################################################################

    def globalfunction_chiral(self,p,msq,asq):
        return p[0] + p[1]*msq + (p[2]*asq + p[3]*np.log(msq))*msq

    def chi2(self,p,ydata,yerr,xdata,x2data):
        #As long as user defines a function they pass this should for any function
        #hopefully
        chisq = 0
        for i in range(len(ydata)):
            #chisq = chisq + pow((1/(yerr[i]))*(ydata[i] - f(xdata[i],p)),2)
            #where p is a list or vector of parameters
            chisq = chisq +pow(self.residuals(p,ydata[i],yerr[i], xdata[i],x2data[i]),2)
            chidof = chisq/len(ydata)
        return chisq, chidof


    def plot_Global(self,p,y,yerr,m,a,yperr,xlabel,ylabel,name):
        '''
        call this from boot strap global
        plot R or B against (m/f)^2
        plot function fitted
        plot physical mass + continuum'''
        """
        mpsq_fpsq_phys = pow(0.13957018/(4*pi*0.13),2)
        mp_phys = 0.13957018
        fp_phys = 0.13
        """
        ydata = []
        mdata = []
        adata = []
        dataerr=[]
        tempy=[]
        tempm=[]
        tempa=[]
        temperr=[]
        for i in range(len(self.betas)):
            if i > 0 and self.betas[i]==self.betas[i-1]:
                dummy = 1
            else:
                if i>0:
                    ydata.append(tempy)
                    mdata.append(tempm)
                    adata.append(tempa)
                    dataerr.append(temperr)
                tempy=[]
                tempm=[]
                tempa=[]
                temperr=[]
            tempy.append(y[i])
            tempm.append(m[i])
            tempa.append(a[i])
            temperr.append(yerr[i])
            if i==len(self.betas)-1:
                ydata.append(tempy)
                mdata.append(tempm)
                adata.append(tempa)
                dataerr.append(temperr)

        plt.figure()
        lab = ['24cubed','32cubed','48cubed','64cubed']
        col = ['r','b','g','y']
        for i in range(4):
            plt.errorbar(mdata[i],ydata[i],yerr=dataerr[i],ecolor=col[i],fmt='v',color=col[i],markersize=5,label=lab[i])

        xmax=plt.xlim()[1]
        xmin=plt.xlim()[0]
        #save limits 

        x = range(0,int(math.ceil(xmax*1.1*100000)))
        #set range over which function plotted.
        f1=[]
        f2=[]
        f3=[]
        f4=[]
        for i in range(len(x)):
            x[i] = x[i]/100000.0
            f1.append(self.globalfunction(p,x[i],adata[0][0]))
            f2.append(self.globalfunction(p,x[i],adata[1][0]))
            f3.append(self.globalfunction(p,x[i],adata[2][0]))
            f4.append(self.globalfunction(p,x[i],adata[3][0]))
        #find the continuum results
        yp = self.globalfunction(p,self.mp,0)
        #p1 = self.globalfunction(p,self.mp,adata[0][0])
        #p2 = self.globalfunction(p,self.mp,adata[1][0])

        #plot the fit lines
        plt.plot(x,f1,'r-')
        plt.plot(x,f2,'b-') 
        plt.plot(x,f3,'g-')
        plt.plot(x,f4,'y-')
        plt.errorbar(self.mp,yp,yerr=yperr,fmt='mo',ecolor='m',markersize=6.5,label='continuum limit phys pt')
        #plt.plot(mp,p1,'rx',markersize=5)
        #plt.plot(mp,p2,'bx',markersize=5)


        ylims=plt.ylim()
        plt.vlines(self.mp,ylims[0],ylims[1],colors='k', linestyles='dashed',label='physical pion mass')
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=3, mode="expand", borderaxespad=0.,prop={'size':7})
        plt.xlabel(xlabel,fontsize=16)
        plt.ylabel(ylabel,fontsize=20)
        plt.ylim(ylims)
        plt.xlim([0,1.05*xmax])
        plt.savefig(name+'.eps',format='eps')
        print "Saving ", name, ".eps"
        plt.close()
