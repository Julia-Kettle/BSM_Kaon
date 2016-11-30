import numpy as np
import random
import math
from scipy import stats
import matplotlib.pyplot as plt
from .fit import Fit

class GlobalFit(Fit):


    def __init__(self,betas,data,m2f2,a2,val,ischeme,ikin,ibasis,dolog):
        """
        set up the data needed for the fitting process
        """
        #Initialises the class with a list of beta indeices, data and masses and ainvs.
        #Can also set the scheme, basis, kinematics and whether lin or chiral fit
        Fit.__init__(self)
        self.x1_phys= pow(0.13957018/(4*math.pi*0.13),2) #physical m/4pif
        self.x2_phys=0
        self.scheme_name = ['MOM','ms'][ischeme]
        self.basis_name = ['SUSY','Lattice'][ibasis]
        self.kin_name = ['E','qq','qg','gq','gg'][ikin]
        self.lattice_name = ['24','32','48','64']
        self.betas = betas
        self.data = data #Either the R or B data
        self.x1 = m2f2 #m^2/4pi f^2
        self.x2 = a2 #a^2
        self.val=val #string either R or B
        self.dolog = dolog #chiral fit if 1
        self.fittype = ["linear","chiral"][dolog]
        if self.val == "R":
            val_no = 0
        elif self.val == "B":
            val_no = 1
        if dolog == 1:
            self.coeff_set = [[[1.5,2.5,2.5,1.5,1.5],[1.5,1.5,1.5,2.5,2.5]],[[-0.5,0.5,0.5,-0.5,-0.5],[-0.5,-0.5,-0.5,0.5,0.5]]][val_no][ibasis]
        else:
            self.coeff_set = np.zeros([5])
        #the coefficients of the log term in the chiral fit


    def SetUp(self):
        """
        Set up the data arrays such that order is [channels,bootstraps,datapoints(ml)]
        """

        print "---Doing Global Fitting---"
        print "Scheme - ",self.scheme_name,"  Basis - ",self.basis_name,"  Kinematics - ",self.kin_name,"  Fit type - ",self.fittype
        self.data=np.swapaxes(self.data,0,2)
        self.data=np.swapaxes(self.data,0,1)
        self.x1=np.squeeze(self.x1)
        self.x1=np.swapaxes(self.x1,0,1)
        self.x2=np.swapaxes(self.x2,0,1)
        self.dataerror=np.std(self.data,axis=1)


    def Run(self):
        """
        Run the global fit bootstrap by bootstrap for all  channels
        """

        #set up arrays to store results from globalfit
        store_params=np.zeros(np.shape(self.data)[:-1]+(3,))
        store_cov=np.zeros(np.shape(self.data)[:-1]+(3,3))
        store_chisq=np.zeros(np.shape(self.data)[:-1])
        store_chidof=np.zeros(np.shape(self.data)[:-1])
        store_y_phys=np.zeros(np.shape(self.data)[:-1])

        #loop over the channels and bootrstraps
        for index in np.ndindex(np.shape(self.data)[:-1]):

            #run the global fit
            store_params[index][:], store_cov[index], store_chisq[index], store_chidof[index] = Fit().globalfit(self.data[index],self.dataerror[index[0]],self.x1[index[-1]],self.x2[index[-1]],[1,1,1],self.coeff_set[index[0]])

            #find the physical result using the params from the fit
            store_y_phys[index] = Fit().globalfunction(store_params[index],self.x1_phys,self.x2_phys,self.coeff_set[index[0]])

        #store central and error values of params and phys point
        params_cent=store_params[:,-1]
        params_err=np.std(store_params,axis=1)
        yperr =np.std(store_y_phys,axis=1)
        #loop over the channels
        for index in np.ndindex(np.shape(self.data)[0]):

            #print out the fitting results
            print "---Fitting Results---"
            print "chisq - ", store_chisq[index][-1], " +/- ", np.std(store_chisq[index])
            print "chisq/dof - ", store_chidof[index][-1], "+/-", np.std(store_chidof[index])
            print "parameters: "
            print "p0 - ", params_cent[0], " +/- ", np.std(store_params[index][:,0])
            print "p1 - ", params_cent[1], " +/- ", np.std(store_params[index][:,1])
            print "p2 - ", params_cent[2], " +/- ", np.std(store_params[index][:,2])

            #set folder based on fit ansatz and filename
            if self.dolog == 0:
                subfolder = "linear/"
            elif self.dolog == 1:
                subfolder = "chiral/"
            name="./plots/" + subfolder +self.val + str(index[0]+1) + "_"+ self.basis_name+"_"+self.scheme_name+"_"+self.kin_name
            self.coeff = self.coeff_set[index[0]]
            #plot the fit results & data points
            self.plot_Global(params_cent[index],self.data[index][-1],self.dataerror[index],self.x1[-1],self.x2[-1],yperr[index],r"$m^2_{\pi}/4\pi f^2_{\pi}$",r"$"+self.val+"_{"+str(index[0]+1)+"}$",name)

        return store_y_phys, store_y_phys[:,-1], np.std(store_y_phys,axis=1)




    def plot_Global(self,params,y,yerr,m,a,yperr,xlabel,ylabel,name):
        """
        Plots data and it's fitted function against x1
        Calculates and plots the physical point
        Saves the plot to a pdf
        """
        #lists to store data
        ydata = []
        mdata = []
        adata = []
        dataerr=[]

        #restructure data into 1 list per data set        
        beta_dict={i:self.betas.count(i) for i in self.betas}
        indx=0
        for k,v in beta_dict.items():
            ydata.append(y[indx:indx+v])
            dataerr.append(yerr[indx:indx+v])
            mdata.append(m[indx:indx+v])
            adata.append(a[indx:indx+v])
            indx=indx+v

        #labelling and colours for each set
        lab = ['24cubed','32cubed','48cubed','64cubed']
        col = ['r','b','g','y']
        plt.figure()

        #plot data with errorbars
        for i in range(4):
            plt.errorbar(mdata[i],ydata[i],yerr=dataerr[i],ecolor=col[i],fmt='v',color=col[i],markersize=5,label=lab[i])

        #save graph limits so not changed by generation of fit lines
        xmax=plt.xlim()[1]
        xmin=plt.xlim()[0]

        #calculate the fit lines over a range
        x = range(0,int(math.ceil(xmax*1.1*100000)))
        f1=[]
        f2=[]
        f3=[]
        f4=[]
        for i in range(len(x)):
            x[i] = x[i]/100000.0
            f1.append(Fit().globalfunction(params,x[i],adata[0][0],self.coeff))
            f2.append(Fit().globalfunction(params,x[i],adata[1][0],self.coeff))
            f3.append(Fit().globalfunction(params,x[i],adata[2][0],self.coeff))
            f4.append(Fit().globalfunction(params,x[i],adata[3][0],self.coeff))

        #find the physical point
        yp = Fit().globalfunction(params,self.x1_phys,self.x2_phys,self.coeff)
        print "---continuum point is: yp = ", str(yp), " +/- ", str(yperr)


        #p1 = self.globalfunction(p,self.mp,adata[0][0])
        #p2 = self.globalfunction(p,self.mp,adata[1][0])

        #plot the fit lines
        plt.plot(x,f1,'r-')
        plt.plot(x,f2,'b-') 
        plt.plot(x,f3,'g-')
        plt.plot(x,f4,'y-')
        plt.errorbar(self.x1_phys,yp,yerr=yperr,fmt='mo',ecolor='m',markersize=6.5,label='continuum limit phys pt')

        ylims=plt.ylim()
        #plt.plot(mp,p1,'rx',markersize=5)
        #plt.plot(mp,p2,'bx',markersize=5)

        #plot a vertical line at the physical mass
        plt.vlines(self.x1_phys,ylims[0],ylims[1],colors='k', linestyles='dashed',label='physical pion mass')
        #plot a key
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=3, mode="expand", borderaxespad=0.,prop={'size':7})

        #set labels
        plt.xlabel(xlabel,fontsize=16)
        plt.ylabel(ylabel,fontsize=20)
        #reset limits with old value
        plt.ylim(ylims)
        plt.xlim([0,1.05*xmax])
        plt.savefig(name+'.pdf',format='pdf')
        print "Saving ", name, ".pdf"
        plt.close()
