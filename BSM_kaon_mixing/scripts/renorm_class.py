from scipy import stats
import config as cnfg
from struct import *
import numpy as np
from copy import copy
from copy import deepcopy
from decimal import *
import fit
from math import pi
import utils
import file_io
import plot


class Renormalisation:

    def __init__(self,lattice,ischeme,ikin,ibasis):
        "Initialise the global variables, and those specific to this lattice"
        print "Renormalising"
        self.ikin = ikin
        self.ischeme = ischeme
        self.ibasis = ibasis
        self.scheme_name = ['MOM','ms'][ischeme]
        self.basis_name = ['SUSY','Lattice'][ibasis]
        self.kin_name = ['E','qq','gg'][ikin]
        self.lattice_name = lattice.name
        self.lattice_dim = lattice.name if 'smeared' not in lattice.name else lattice.name[0:2]
        cnfg.nboots = 500
        self.N_ren = [8.0/3,4.0/3,-2.0,5.0/3,1.0] #This is the factor in the bag parameter
        self.N = [8.0/3,-5.0/3,1.0/3,2.0,2.0/3]
        cnfg.nchan=5
        self.T = [ [np.array([[1,0,0,0,0],[0,0,0,1,0],[0,0,0,-0.5,0.5],[0,0,1,0,0],[0,-0.5,0,0,0]])
], [np.identity(5) ] ][ibasis]
        self.N=np.squeeze(np.dot(self.T,self.N_ren))
        self.Mkphys = 0.5*(493.677+497.614)*pow(10,-3)
        self.Fkphys = 156.2*pow(10,-3)
        cnfg.nboots = 500
        lattice_name_Zs = {'24':'24cubed','24smeared':'24cubed','32':'32cubed','48':'24cubed','48smeared':'24cubed','64':'32cubed','48fine':'48fine'}
        Zfilename = "../data/2GeV/Z" + self.scheme_name + "_boot_mu_match_" + lattice_name_Zs[self.lattice_name] + "_" + self.kin_name + "_block.out"
        self.bZ_scheme=np.zeros([cnfg.nchan,cnfg.nchan,101])
        file_io.read_file_array(Zfilename,self.bZ_scheme)
        dictZA = {'24': 0.71273,'32':0.74404 ,'48':0.71273, '64':0.74404,'48fine':0.76104 }
        #dictZP = {'24MOME':0.57764,'32MOME':0.57773,'24msE':0.69512,'32msE':0.69464,'24MOMqq':0.65635,'32MOMqq':0.65832,'24msqq':0.6563*1.0499,'32msqq':0.65832*1.04993,'24MOMgg':0.69513,'32MOMgg':0.69471,'24msgg':0.69513*1.0157,'32msgg':0.69471*1.0157,'48MOME':0.57764,'64MOME':0.57773,'48msE':0.69512,'64msE':0.69464,'48MOMqq':0.65635,'64MOMqq':0.65864,'48msqq':0.6563*1.0499,'64msqq':0.65864*1.04993,'48MOMgg':0.69513,'64MOMgg':0.69471,'48msgg':0.69513*1.0157,'64msgg':0.69471*1.0157}
        #Below is Zp at 2Gev - there so I can test with the 48fine at 2GeV BUT I have no values for 48 fine. 
        dictZP = {'24MOME':0.57764,'32MOME':0.57773,'24msE':0.69512,'32msE':0.69464,'24MOMqq':0.6423,'32MOMqq':0.6940,'24msqq':0.6372,'32msqq':0.6506,'24MOMgg':0.5974,'32MOMgg':00.6585,'24msgg':0.5924,'32msgg':0.6320,'48MOME':0.57764,'64MOME':0.57773,'48fineMOME':1,'48msE':0.69512,'64msE':0.69464,'48finemsE':1,'48MOMqq':0.6423,'64MOMqq':0.6940,'48fineMOMqq':1,'48msqq':0.6372,'64msqq':0.6506,'48finemsqq':1,'48MOMgg':0.5974,'64MOMgg':00.6585,'48fineMOMgg':1,'48msgg':0.5924,'64msgg':0.6320,'48finemsgg':1,}
        self.ZP = dictZP[self.lattice_dim+self.scheme_name+self.kin_name]
        self.ZA = dictZA[self.lattice_dim]

    def Renorm_R(self,R):
        #(chan,ml,ms,nboot)
        #Scale the ratios R by Z factor
        R = np.swapaxes(R,0,-1)
        dimR = np.shape(R)
        ren_R = np.zeros(dimR)
        ren_phys_R = np.zeros(dimR)
        for index in np.ndindex(np.shape(R)[:-1]):
            ren_R[index] = np.dot(self.bZ_scheme[:,:,-1],R[index])/self.bZ_scheme[0][0][-1]
            ren_phys_R[index] = np.dot(self.T,ren_R[index])
        ren_R = np.swapaxes(ren_R,0,-1)
        ren_phys_R = np.swapaxes(ren_phys_R,0,-1)
        R = np.swapaxes(R,0,-1)
        return ren_phys_R
        #print "Renormalised R is- \n",self.ren_phys_R[:,0,0,-1]

    def Renorm_B(self,B):
        #Scale the bag B by the Z factors
        dimB = np.shape(B)
        B = np.swapaxes(B,0,-1)
        ren_B = np.swapaxes(np.zeros(dimB),0,-1)
        ren_phys_B = np.swapaxes(np.zeros(dimB),0,-1)
        B_phys = np.swapaxes(np.zeros(dimB),0,-1)
        for index in np.ndindex(np.shape(B)[:-1]):
            ren_B[index] = (self.ZA*self.ZA)*(np.dot(self.bZ_scheme[:,:,-1],B[index]))
            ren_phys_B[index] = np.dot(self.T,ren_B[index])
            for ich in range(1,5):
                B_phys[index][ich] = ren_phys_B[index][ich]/(self.ZP*self.ZP*self.N[ich])
            B_phys[index][0] = ren_phys_B[index][0]/(self.ZA*self.ZA*self.N[0])

        B = np.swapaxes(B,0,-1)
        ren_B = np.swapaxes(ren_B,0,-1)
        ren_phys_B = np.swapaxes(ren_phys_B,0,-1)
        B_phys = np.swapaxes(B_phys,0,-1)
        return B_phys

    #this is our main function
    #we want to pass the class/function (do we need a class?) Incorporate into the renormalisation??
    #pass R/B - loop through the channels,lights,mval2, to get a matrix [strange,nboot]
    #then loop through the boots interpolating


class StrangeAdjustment:

    def __init__(self,lattice):
        self.lattice = lattice
        self.m_s = copy(lattice.m_val_s)
        self.m_s_phys = lattice.m_s_phys
        self.m_l = copy(lattice.m_sea_l)
        #(channels,ml,ms,nboot)


    def interpolation(self,ydata):
        #put in form (channels,ml,nboot,ms)
        ydata=np.swapaxes(ydata,2,3)
        y_int = np.zeros(np.shape(ydata)[:-1])
        dydms = np.zeros(np.shape(ydata)[:-1])
        for index in np.ndindex(np.shape(ydata)[:-2]):
            for iboot in range(np.shape(ydata)[-2]):
                y = ydata[index+(iboot,)]
                yerr = np.std(ydata[index],axis=0)
                y_int[index+(iboot,)] , dydms[index+(iboot,)] = fit.interp(self.m_s,y,yerr,self.m_s_phys[0],1)
        return y_int, dydms

    def slope_extrap(self,dydms,ml_phys):
        print self.m_l, ml_phys[0]
        #put in form (channels,nboot,ml)
        dydms = np.swapaxes(dydms,1,2)
        dydms_int = np.zeros(np.shape(dydms)[:-1])
        #print dydms[:,-1,:]
        for index in np.ndindex(np.shape(dydms)[:-2]):
            for iboot in range(np.shape(dydms)[-2]):
                y = dydms[index+(iboot,)]
                yerr = np.std(dydms[index],axis=0)
                #print "yerr", yerr
                dydms_int[index+(iboot,)], dummy = fit.interp(self.m_l,y,yerr,ml_phys[0],1)
                #print "dydms_interr", np.std(dydms_int[index])
            xlabel="aml"
            ylabel="slope"
            name="slope"+str(index+(1,))+"_"+self.lattice.name
            #plot.plot_interp(dydms[index,-1],np.std(dydms[index]),self.m_l,dydms_int[index,-1],np.std(dydms_int[index]),ml_phys[0],dummy,name,xlabel,ylabel)

        return dydms_int

    def extrap_1pt(self,ydata,dydms):
        #put in form (channels, ml,nboot,ms)
        ydata=np.swapaxes(ydata,2,3)
        
        y_int = np.zeros(np.shape(ydata)[:-1])
        for index in np.ndindex(np.shape(ydata)[:-1]):
            #print index, index[0], index[-1], index+(0,)
            y_int[index] = ydata[index+(0,)] + dydms[index[0],index[-1]]*(+self.m_s_phys[0]-self.m_s[0])
        return y_int




