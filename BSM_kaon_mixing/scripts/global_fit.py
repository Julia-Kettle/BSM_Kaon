import numpy as np
import random
import math
from scipy import stats
import matplotlib.pyplot as plt
import fit

def str_with_err(value, error):
    sign=np.sign(value)
    value=abs(value)
    digits = -int(math.floor(math.log10(error)))+1
    prepnt = len(str(int(value//1)))
    return "{0:{pre}.0f}.{1:{dig}.0f}({2:.0f})".format(sign*int(value//1),value%1*10**digits, error*10**digits,pre=prepnt,dig=digits)

class GlobalFit():


    def __init__(self,betas,data,m2f2,a2,val,ischeme,ikin,ibasis,dolog):
        """
        set up the data needed for the fitting process
        """
        #Initialises the class with a list of beta indeices, data and masses and ainvs.
        #Can also set the scheme, basis, kinematics and whether lin or chiral fit
        self.x1_phys= pow(0.13957018/(4*math.pi*0.13),2) #physical m/4pif
        self.x2_phys=0
        self.scheme_name = ['MOM','ms'][ischeme]
        self.basis_name = ['SUSY','Lattice'][ibasis]
        self.kin_name = ['E','qq','gg'][ikin]
        self.lattice_name = ['24','32','48','64','48fine']
        self.betas = betas
        self.data = data #Either the R or B data
        self.x1 = m2f2 #m^2/4pi f^2
        self.x2 = a2 #a^2
        self.val=val #string either R or B
        self.dolog = dolog #chiral fit if 1 both if 2
        self.fittype = ["linear","chiral","both"][dolog]
        if self.val == "R":
            val_no = 0
        elif self.val == "B":
            val_no = 1
        if dolog == 1:
            self.coeff_set = [[[1.5,2.5,2.5,1.5,1.5],[1.5,1.5,1.5,2.5,2.5]],[[-0.5,0.5,0.5,-0.5,-0.5],[-0.5,-0.5,-0.5,0.5,0.5]]][val_no][ibasis]
        elif dolog == 0:
            self.coeff_set = np.zeros([5])
        else: self.coeff_set = [ np.zeros([5]),[[[1.5,2.5,2.5,1.5,1.5],[1.5,1.5,1.5,2.5,2.5]],[[-0.5,0.5,0.5,-0.5,-0.5],[-0.5,-0.5,-0.5,0.5,0.5]]][val_no][ibasis]
 ]

        #the coefficients of the log term in the chiral fit


    def SetUp(self):
        """
        Set up the data arrays such that order is [channels,bootstraps,datapoints(ml)] by swapping axes. 
        Currently  [datapoints,channels,bootstraps]
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
        if self.dolog != 2:
            store_params=np.zeros(np.shape(self.data)[:-1]+(3,))
            store_cov=np.zeros(np.shape(self.data)[:-1]+(3,3))
            store_chisq=np.zeros(np.shape(self.data)[:-1])
            store_chidof=np.zeros(np.shape(self.data)[:-1])
            store_y_phys=np.zeros(np.shape(self.data)[:-1])
        else:
            store_params_lin=np.zeros(np.shape(self.data)[:-1]+(3,))
            store_cov_lin=np.zeros(np.shape(self.data)[:-1]+(3,3))
            store_y_phys_lin=np.zeros(np.shape(self.data)[:-1])
            store_chisq_lin=np.zeros(np.shape(self.data)[:-1])
            store_chidof_lin=np.zeros(np.shape(self.data)[:-1])

            store_params_chir=np.zeros(np.shape(self.data)[:-1]+(3,))
            store_cov_chir=np.zeros(np.shape(self.data)[:-1]+(3,3))
            store_y_phys_chir=np.zeros(np.shape(self.data)[:-1])
            store_chisq_chir=np.zeros(np.shape(self.data)[:-1])
            store_chidof_chir=np.zeros(np.shape(self.data)[:-1])

        #loop over the channels and bootrstraps
        for index in np.ndindex(np.shape(self.data)[:-1]):

            #run the global fit
            if self.dolog != 2:
                store_params[index][:], store_cov[index], store_chisq[index], store_chidof[index] = fit.nonlinearfit(self.data[index],self.dataerror[index[0]],self.x1[index[-1]],self.x2[index[-1]],[1,1,1],self.coeff_set[index[0]])
                store_y_phys[index] = fit.globalfunction(store_params[index],self.x1_phys,self.x2_phys,self.coeff_set[index[0]])
            else:
                #linear fit
                store_params_lin[index][:], store_cov_lin[index], store_chisq_lin[index], store_chidof_lin[index] = fit.nonlinearfit(self.data[index],self.dataerror[index[0]],self.x1[index[-1]],self.x2[index[-1]],[1,1,1],self.coeff_set[0][index[0]])
                #chiralfit
                store_params_chir[index][:], store_cov_chir[index], store_chisq_chir[index], store_chidof_chir[index] = fit.nonlinearfit(self.data[index],self.dataerror[index[0]],self.x1[index[-1]],self.x2[index[-1]],[1,1,1],self.coeff_set[1][index[0]])
                #get yphys from fit results
                store_y_phys_lin[index] = fit.globalfunction(store_params_lin[index],self.x1_phys,self.x2_phys,self.coeff_set[0][index[0]])
                store_y_phys_chir[index] = fit.globalfunction(store_params_chir[index],self.x1_phys,self.x2_phys,self.coeff_set[1][index[0]])

        #store central and error values of params and phys point
            if self.dolog !=2:
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
            else:
                params_cent_lin=store_params_lin[:,-1]
                params_err_lin=np.std(store_params_lin,axis=1)
                yperr_lin =np.std(store_y_phys_lin,axis=1)
                params_cent_chir=store_params_chir[:,-1]
                params_err_chir=np.std(store_params_chir,axis=1)
                yperr_chir =np.std(store_y_phys_chir,axis=1)
                #loop over the channels
        print "& phys pt (lin) & $\chi^2$ /dof & phys pt (chir) & $\chi^2$ /dof \\"+"\\"
        print "\hline"
        for index in np.ndindex(np.shape(self.data)[0]):
            print "$"+str(self.val)+"_"+str(index[0]+1)+"$ & "+ str_with_err(store_y_phys_lin[index][-1],yperr_lin[index]) + " & " + str_with_err(store_chidof_lin[index][-1],np.std(store_chidof_lin[index])) + "&"+ str_with_err(store_y_phys_chir[index][-1],yperr_chir[index]) + " & " + str_with_err(store_chidof_chir[index][-1],np.std(store_chidof_chir[index])) + "\\" + "\\"
        for index in np.ndindex(np.shape(self.data)[0]):
            #print out the fitting results
            print "------Fitting Results------"
            print "---linear fit---"
            print "chisq - ", store_chisq_lin[index][-1], " +/- ", np.std(store_chisq_lin[index])
            print "chisq/dof - ", store_chidof_lin[index][-1], "+/-", np.std(store_chidof_lin[index])
            print "parameters: "
            print "p0 - ", params_cent_lin[0], " +/- ", np.std(store_params_lin[index][:,0])
            print "p1 - ", params_cent_lin[1], " +/- ", np.std(store_params_lin[index][:,1])
            print "p2 - ", params_cent_lin[2], " +/- ", np.std(store_params_lin[index][:,2])
            print "---chiral fit---"
            print "chisq - ", store_chisq_chir[index][-1], " +/- ", np.std(store_chisq_chir[index])
            print "chisq/dof - ", store_chidof_chir[index][-1], "+/-", np.std(store_chidof_chir[index])
            print "parameters: "
            print "p0 - ", params_cent_chir[0], " +/- ", np.std(store_params_chir[index][:,0])
            print "p1 - ", params_cent_chir[1], " +/- ", np.std(store_params_chir[index][:,1])
            print "p2 - ", params_cent_chir[2], " +/- ", np.std(store_params_chir[index][:,2])
            #set folder based on fit ansatz and filename
            if self.dolog == 0:
                subfolder = "linear/"
            elif self.dolog == 1:
                subfolder = "chiral/"
            elif self.dolog == 2:
                subfolder = "full/"
            name="../plots/" + subfolder +self.val + str(index[0]+1) + "_"+ self.basis_name+"_"+self.scheme_name+"_"+self.kin_name
            self.coeff_lin = self.coeff_set[0][index[0]]
            self.coeff_chir = self.coeff_set[1][index[0]]
            #plot the fit results & data points
            if self.dolog != 2:
                self.plot_Global(params_cent[index],self.data[index][-1],self.dataerror[index],self.x1[-1],self.x2[-1],yperr[index],r"$m^2_{\pi}/(4\pi f_{\pi})^2$",r"$"+self.val+"_{"+str(index[0]+1)+"}$",name)
            else:
                self.plot_Global_full(params_cent_lin[index],params_cent_chir[index],self.data[index][-1],self.dataerror[index],self.x1[-1],self.x2[-1],yperr_lin[index],yperr_chir[index],r"$m^2_{\pi}/(4\pi f_{\pi})^2$",r"$"+self.val+"_{"+str(index[0]+1)+"}$",name)
 
        return 1


    def plot_Global_full(self,params_lin,params_chir,y,yerr,m,a,yperr_lin,yperr_chir,xlabel,ylabel,name):
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
        lab = [r'$24^3$',r'$32^3$',r'$48^3$',r'$64^3$',r'$48^3 fine$']
        col = ['r','b','g','y','orange']
        plt.figure()
        plt.subplot(1,1,1)


        #find the physical point
        yp_lin = fit.globalfunction(params_lin,self.x1_phys,self.x2_phys,self.coeff_lin)
        yp_chir = fit.globalfunction(params_chir,self.x1_phys,self.x2_phys,self.coeff_chir)
        print "---linear continuum point is: yp = ", str(yp_lin), " +/- ", str(yperr_lin)
        print "---chiral continuum point is: yp = ", str(yp_chir), " +/- ", str(yperr_chir)

        #plot data with errorbars
        for i in range(5):
            try:
                plt.errorbar(mdata[i],ydata[i],yerr=dataerr[i],ecolor=col[i],fmt='v',color=col[i],markersize=5,label=lab[i])
                print dataerr[i]
            except:
                print "beta ", i, " plot failed"
                pass
        #save graph limits so not changed by generation of fit lines
        xmax=plt.xlim()[1]
        xmin=plt.xlim()[0]

        #calculate the fit lines over a range
        x = range(0,int(math.ceil(xmax*1.1*1000)))
        f1_lin=[]
        f2_lin=[]
        f3_lin=[]
        f4_lin=[]
        f5_lin=[]
        f1_chir=[]
        f2_chir=[]
        f3_chir=[]
        f4_chir=[]
        f5_chir=[]
        f_chir=[]
        f_lin=[]
        f_chir_min=[]
        f_chir_max=[]
        for i in range(len(x)):
            x[i] = x[i]/1000.0
            f_chir.append(fit.globalfunction(params_chir,x[i],0,self.coeff_chir))
            f_lin.append(fit.globalfunction(params_lin,x[i],0,self.coeff_lin))
            try:
                f1_lin.append(fit.globalfunction(params_lin,x[i],adata[0][0],self.coeff_lin))
                f1_chir.append(fit.globalfunction(params_chir,x[i],adata[0][0],self.coeff_chir))
            except:
                pass
            try:
                f2_lin.append(fit.globalfunction(params_lin,x[i],adata[1][0],self.coeff_lin))
                f2_chir.append(fit.globalfunction(params_chir,x[i],adata[1][0],self.coeff_chir))
            except:
                pass
            try:
                f3_lin.append(fit.globalfunction(params_lin,x[i],adata[2][0],self.coeff_lin))
                f3_chir.append(fit.globalfunction(params_chir,x[i],adata[2][0],self.coeff_chir))
            except:
                pass
            try:
                f4_lin.append(fit.globalfunction(params_lin,x[i],adata[3][0],self.coeff_lin))
                f4_chir.append(fit.globalfunction(params_chir,x[i],adata[3][0],self.coeff_chir))
            except:
                pass
            try:
                f5_lin.append(fit.globalfunction(params_lin,x[i],adata[4][0],self.coeff_lin))
                f5_chir.append(fit.globalfunction(params_chir,x[i],adata[4][0],self.coeff_chir))
            except:
                pass
        print "fits calculalated"

        #p1 = self.globalfunction(p,self.mp,adata[0][0])
        #p2 = self.globalfunction(p,self.mp,adata[1][0])

        #plot the fit lines
        plt.plot(x,f_chir,'c')
        plt.plot(x,f_lin,'m')
        print "physical fit lines plotted"
        try:
            plt.plot(x,f1_lin,'r--')
            #plt.plot(x,f1_chir,'r--')
        except:
            pass
        try:
            plt.plot(x,f2_lin,'b--') 
            #plt.plot(x,f2_chir,'b--') 
        except:
            pass
        try:
            plt.plot(x,f3_lin,'g--')
            #plt.plot(x,f3_chir,'g--')
        except:
            pass
        try:
            plt.plot(x,f4_lin,'y--')
            #plt.plot(x,f4_chir,'y--')
        except:
            pass
        try:
            plt.plot(x,f5_lin,'k--')
            #plt.plot(x,f5_chir,'orange--')
        except:
            pass
        print "linear fits plotted"
        plt.errorbar(self.x1_phys,yp_lin,yerr=yperr_lin,fmt='mo',ecolor='m',markersize=6.5,label=r'phys pt extrapolation (linear)')
        plt.errorbar(self.x1_phys,yp_chir,yerr=yperr_chir,fmt='co',ecolor='c',markersize=6.5,label=r'phys pt extrapolation (chiral)')

        ylims=plt.ylim()
        #plt.plot(mp,p1,'rx',markersize=5)
        #plt.plot(mp,p2,'bx',markersize=5)
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.,prop={'size':15})
        #plot a vertical line at the physical mass
        plt.vlines(self.x1_phys,ylims[0],ylims[1],colors='k', linestyles='dashed',label='physical pion mass')
        #plot a key
        

        #set labels
        plt.xlabel(xlabel,fontsize=16)
        plt.ylabel(ylabel,fontsize=20)
        #reset limits with old value
        plt.ylim(ylims)
        plt.xlim([0,1.05*xmax])
        plt.subplots_adjust(top=0.80)
        plt.savefig(name+'.pdf',format='pdf')
        print "Saving ", name, ".pdf"
        plt.close()
        print "plotted and closed"


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
            try:
                plt.errorbar(mdata[i],ydata[i],yerr=dataerr[i],ecolor=col[i],fmt='v',color=col[i],markersize=5,label=lab[i])
            except:
                pass
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
            try:
                f1.append(fit.globalfunction(params,x[i],adata[0][0],self.coeff))
            except:
                pass
            try:
                f2.append(fit.globalfunction(params,x[i],adata[1][0],self.coeff))
            except:
                pass
            try:
                f3.append(fit.globalfunction(params,x[i],adata[2][0],self.coeff))
            except:
                pass
            try:
                f4.append(fit.globalfunction(params,x[i],adata[3][0],self.coeff))
            except:
                pass

        #find the physical point
        yp = fit.globalfunction(params,self.x1_phys,self.x2_phys,self.coeff)
        print "---continuum point is: yp = ", str(yp), " +/- ", str(yperr)


        #p1 = self.globalfunction(p,self.mp,adata[0][0])
        #p2 = self.globalfunction(p,self.mp,adata[1][0])

        #plot the fit lines
        try:
            plt.plot(x,f1,'r-')
        except:
            pass
        try:
            plt.plot(x,f2,'b-') 
        except:
            pass
        try:
            plt.plot(x,f3,'g-')
        except:
            pass
        try:
            plt.plot(x,f4,'y-')
        except:
            pass
        plt.errorbar(self.x1_phys,yp,yerr=yperr,fmt='mo',ecolor='m',markersize=6.5,label='continuum limit phys pt')

        ylims=plt.ylim()
        #plt.plot(mp,p1,'rx',markersize=5)
        #plt.plot(mp,p2,'bx',markersize=5)

        #plot a vertical line at the physical mass
        plt.vlines(self.x1_phys,ylims[0],ylims[1],colors='k', linestyles='dashed',label='physical pion mass')
        #plot a key
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=3, mode="expand", borderaxespad=0.,prop={'size':9})

        #set labels
        plt.xlabel(xlabel,fontsize=16)
        plt.ylabel(ylabel,fontsize=20)
        #reset limits with old value
        plt.ylim(ylims)
        plt.xlim([0,1.05*xmax])
        plt.savefig(name+'.pdf',format='pdf')
        print "Saving ", name, ".pdf"
        plt.close()
