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


# I want to restructure this a bit. 
# I have the same renormalisation factors for B and R so a class makes sense. 
# Shared values: ZBK, ZS 
# let's have static methods to read renorm factors. 

class Renormalisation:

    def __init__(self,zbk,zs,zv,basis):
        self.zbk    =   zbk
        self.zs     =   zs
        self.zv     =   zv
        self.nboots =   len(zs)
        self.basisSetUp(basis)

    def basisSetUp(self,basis):
        # set up the the normalisation factors and transformation matrix depending on basis
        if(basis):
            self.T=np.array([[1,0,0,0,0],[0,0,0,1,0],[0,0,0,-0.5,0.5],[0,0,1,0,0],[0,-0.5,0,0,0]])
        else:
            self.T=np.identity(5)
        self.N=np.squeeze(np.dot(self.T,np.array([8.0/3,4.0/3,-2.0,5.0/3,1.0]) ))

    def renorm_B(self,B):

        # check dimensions, if 4th dimension due to multiple strange, renormalise each strange mass.
        if len(np.shape(B)) == 4:
            renB=self.renorm_multiStrange(B,self.renorm_B)
        
        elif len(np.shape(B)) == 3:
            renB=np.empty_like(B)  #empty renorm matrix to fill
            for il in range(np.size(B,1)):
                for iboot in range(np.size(B,2)):
                    renB[:,il,iboot] = np.dot(self.T, np.dot(self.zbk[:,:,iboot],B[:,il,iboot]) )
                    renB[:,il,iboot] /= self.N   # divide by normalisation 
                    # BSM elements multiplied by (ZV/ZS)^2 : Zij/ZV^2 -> Zij/ZS^2 due to different defs of matrix elems
                    for ich in range(1,5):
                        renB[ich,il,iboot] *= pow(self.zv[iboot],2)/pow(self.zs[iboot],2)
        else:
            print "Error: array should have dimensions (channels,n_lightmass,boots) or (channels,n_lightmass,n_strangemass,boots)"
            sys.exit(-2)
        return renB

    def renorm_R(self,R):
        #renormalizes R by multiplying with renormalisation Z_Bk.R[ich]
        if len(np.shape(R)) == 4:
            renR=self.renorm_multiStrange(R,self.renorm_R)
        elif len(np.shape(R)) == 3:
            renR=np.empty_like(R)
            for il in range(np.size(R,1)):
                for iboot in range(np.size(R,2)):
                    renR[:,il,iboot] = np.dot(self.T, np.dot(self.zbk[:,:,iboot],R[:,il,iboot])/self.zbk[0,0,iboot])
        else:
            print "Error: array should have dimensions (channels,n_lightmass,boots) or (channels,n_lightmass,n_strangemass,boots)"
            sys.exit(-1)
        return renR

    def renorm_multiStrange(self,mat,renFunc):
        renmat = np.empty_like(mat)
        for ims in range(np.size(mat,2)):
            renmat[:,:,ims,:] = renFunc(mat[:,:,ims,:])
        return renmat

    @staticmethod
    def readZBK(renlattice_name,scheme,kin,directory):
        finame  =   directory + "/Z" + scheme + "_boot_mu_match_" + renlattice_name + "_" + kin + "_block.out"
        fi      =   open(finame,'r')
        lines   =   fi.readlines()
        nboot   =   len(lines)/25 - 1  # where we have 5 by 5 by (nboots + 1)
        zbk     =   np.zeros([5,5,nboot+1])  
        i,j     =   0,0 
        iline   =   0
        while iline < len(lines):
            for iboot in range(nboot+1):
                zbk[i,j,iboot]  =   float(lines[iline])
                iline+=1
            # increase j by one, if j=5 increase i by 1 and j -> 0
            j+=1
            i+=(j)//5
            j=j%5
        return zbk
    
    @staticmethod
    def readZS(renlattice_name,scheme,kin,directory):
        finame  =   directory + "/Zs_" + scheme + "_boot_" + renlattice_name + "_" + kin + ".out"
        fi      =   open(finame,'r')
        lines   =   fi.readlines()
        zs      =   np.zeros(len(lines))
        for i in range(len(lines)):
            zs[i]   =   float(lines[i])
        return zs

    @staticmethod
    def readZV(renlattice_name,directory):
        finame  =   directory + "/Zv_boot_" + renlattice_name + ".out"
        fi      =   open(finame,'r')
        lines   =   fi.readlines()
        zv      =   np.zeros(len(lines))
        for i in range(len(lines)):
            zv[i]   =   float(lines[i])
        return zv

"""
class Renormalisation:

    def __init__(self,lattice,ischeme,ikin,ibasis,do2GeV=False):
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
        self.do2GeV=do2GeV
        cnfg.nboots = 500
        cnfg.nchan=5

        #define the normalisations N, dependent on basis
        self.N_ren = [8.0/3,4.0/3,-2.0,5.0/3,1.0] #This is the factor in the bag parameter
        #self.N = [8.0/3,-5.0/3,1.0/3,2.0,2.0/3]
        self.T = [ [np.array([[1,0,0,0,0],[0,0,0,1,0],[0,0,0,-0.5,0.5],[0,0,1,0,0],[0,-0.5,0,0,0]])], [np.identity(5) ] ][ibasis]
        self.N=np.squeeze(np.dot(self.T,self.N_ren))

        self.Mkphys = 0.5*(493.677+497.614)*pow(10,-3)
        self.Fkphys = 156.2*pow(10,-3)

        self.read_renorm_factors()


        print "Z_bk - "
        for i in range(5):
            print self.bZ_scheme[i,0,-1],self.bZ_scheme[i,1,-1],self.bZ_scheme[i,2,-1],self.bZ_scheme[i,3,-1],self.bZ_scheme[i,4,-1]

        print "Z_bk/(Z[0,0]*ZP^2) - "
        for i in range(5):
            print self.bZ_scheme[i,0,-1]/pow(self.ZP[-1],2),self.bZ_scheme[i,1,-1]/pow(self.ZP[-1],2),self.bZ_scheme[i,2,-1]/pow(self.ZP[-1],2),self.bZ_scheme[i,3,-1]/pow(self.ZP[-1],2),self.bZ_scheme[i,4,-1]/pow(self.ZP[-1],2)
        
        try:
            print("Z_V/A - ",self.ZA[-1])
        except:
            print ""
        try:
            print("Z_S/P - ",self.ZP[-1])
        except:
            print ""


    def read_renorm_factors(self):
        lattice_name_Zs = {'24':'24cubed','24smeared':'24cubed','32':'32cubed','48':'24cubed','48smeared':'24cubed','64':'32cubed','48fine':'48fine','64smeared':'32cubed','32smeared':'32cubed','48finesmeared':'48fine'}
        Zvfilename="../data/Zv_boot_" + lattice_name_Zs[self.lattice_name]+".out"
        Zsfilename="../data/Zs_" + self.scheme_name + "_boot_" + lattice_name_Zs[self.lattice_name] + "_" + self.kin_name + ".out"
        #Zfilename = "../data/2GeV/Z" + self.scheme_name + "_boot_mu_match_" + lattice_name_Zs[self.lattice_name] + "_" + self.kin_name + "_block.out"
        self.bZ_scheme=np.zeros([cnfg.nchan,cnfg.nchan,501])
        self.ZA = np.zeros([501])
        self.ZP = np.zeros([501])

        #try: #hack to get old Z
        doOld = False
        #hack to use 1609 version
        doSUSYConv = False
        if ( '24' in Zsfilename or '32' in Zsfilename)  and (doOld) and self.kin_name=='qq':
            Zfilename = "../data/Z" + self.scheme_name + "_boot_mu_match_" + lattice_name_Zs[self.lattice_name] + "_" + self.kin_name + "_block.out"
            self.bZShort=np.zeros([cnfg.nchan,cnfg.nchan,101])
            file_io.read_file_array(Zfilename,self.bZShort)
            self.bZ_scheme  =   np.dstack((self.bZShort[:,:,:-1],self.bZShort[:,:,:-1],self.bZShort[:,:,:-1],self.bZShort[:,:,:-1],self.bZShort))
        elif ( '24' in Zsfilename or '32' in Zsfilename)  and (doSUSYConv) and self.kin_name=='qq' and self.scheme_name == 'MOM':
            print "Doing SUSY CONV 1609"
            Zfilename = "../data/Z" + self.scheme_name + "_boot_mu_match_" + lattice_name_Zs[self.lattice_name] + "_" + self.kin_name + "_block_500_Lattice.out"
            file_io.read_file_array(Zfilename,self.bZ_scheme)
        elif (self.do2GeV == True):
            Zfilename = "../data/2GeV/Z" + self.scheme_name + "_boot_mu_match_" + lattice_name_Zs[self.lattice_name] + "_" + self.kin_name + "_block.out"
            self.bZShort=np.zeros([cnfg.nchan,cnfg.nchan,101])
            file_io.read_file_array(Zfilename,self.bZShort)
            self.bZ_scheme  =   np.dstack((self.bZShort[:,:,:-1],self.bZShort[:,:,:-1],self.bZShort[:,:,:-1],self.bZShort[:,:,:-1],self.bZShort))
        else:
            print "Doing Lattice 1708"
            Zfilename = "../data/Z" + self.scheme_name + "_boot_mu_match_" + lattice_name_Zs[self.lattice_name] + "_" + self.kin_name + "_block_500.out"
            file_io.read_file_array(Zfilename,self.bZ_scheme)
        file_io.read_file_array(Zvfilename,self.ZA)
        file_io.read_file_array(Zsfilename,self.ZP)


    def renorm_R(self,R):
        #renormalizes R by multiplying with renormalisation Z_Bk.R[ich]
        if len(np.shape(R)) == 4:
            ren_R=self.renorm_R_mult_str(R)
        elif len(np.shape(R)) == 3:
            ren_R=np.empty_like(R)
            for il in range(np.size(R,1)):
                for iboot in range(np.size(R,2)):
                    #ren_R[:,il,iboot] = np.dot(self.bZ_scheme[:,:,-1],R[:,il,iboot])/self.bZ_scheme[0,0,-1]
                    #ren_R[:,il,iboot] = np.dot(self.T,ren_R[:,il,iboot])
                    ren_R[:,il,iboot] = np.dot(self.bZ_scheme[:,:,iboot],R[:,il,iboot])/self.bZ_scheme[0,0,iboot]
                    ren_R[:,il,iboot] = np.dot(self.T,ren_R[:,il,iboot])
        else:
            print "Error: array should have dimensions (channels,n_lightmass,boots) or (channels,n_lightmass,n_strangemass,boots)"
            sys.exit(-1)
        return ren_R

    def renorm_R_mult_str(self,R):
        ren_R = np.empty_like(R)
        for ims in range(np.size(R,2)):
            ren_R[:,:,ims,:] = self.renorm_R(R[:,:,ims,:])
        return ren_R

    def renorm_B_mult_str(self,B):
        ren_B = np.empty_like(B)
        for ims in range(np.size(B,2)):
            ren_B[:,:,ims,:] = self.renorm_B(B[:,:,ims,:])
        return ren_B

    def renorm_B(self,B):
        if len(np.shape(B)) == 4:
            ren_B=self.renorm_B_mult_str(B)
        elif len(np.shape(B)) == 3:
            ren_B=np.empty_like(B)
            for il in range(np.size(B,1)):
                for iboot in range(np.size(B,2)):
                    ren_B[:,il,iboot] = np.dot(self.bZ_scheme[:,:,iboot],B[:,il,iboot])
                    ren_B[:,il,iboot] = np.dot(self.T,ren_B[:,il,iboot])
                    for ich in range(5):
                        if ich == 0:
                            ren_B[ich,il,iboot] /= self.N[ich]
                        else:
                            ren_B[ich,il,iboot] *= ((pow(self.ZA[iboot],2))/(pow(self.ZP[iboot],2)*self.N[ich]))
                print "B before renorm - ", B[:,il,-1]/self.N
                print "B after renrom  - ", ren_B[:,il,-1]
        else:
            print "Error: array should have dimensions (channels,n_lightmass,boots) or (channels,n_lightmass,n_strangemass,boots)"
            sys.exit(-2)

        return ren_B
"""


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


'''
        dictZA = {'24': 0.71273,'32':0.74404 ,'48':0.71273, '64':0.74404,'48fine':0.761108 }
        #dictZP = {'24MOME':0.57764,'32MOME':0.57773,'24msE':0.69512,'32msE':0.69464,'24MOMqq':0.65635,'32MOMqq':0.65832,'24msqq':0.6563*1.0499,'32msqq':0.65832*1.04993,'24MOMgg':0.6563,'32MOMgg':0.6585,'24msgg':0.69513*1.0157,'32msgg':0.69471*1.0157,'48MOME':0.57764,'64MOME':0.57773,'48msE':0.69512,'64msE':0.69464,'48MOMqq':0.65635,'64MOMqq':0.65864,'48msqq':0.6563*1.0499,'64msqq':0.65864*1.04993,'48MOMgg':0.69513,'64MOMgg':0.69471,'48msgg':0.69513*1.0157,'64msgg':0.69471*1.0157,'48fineMOME':1,'48finemsE':1,'48fineMOMqq':1,'48fineMOMgg':0.663341,'48finemsqq':1,'48finemsgg':0.685909,}
        
        dictZP = {'24MOME':0.57764,'24MOMgg':0.6563,'24MOMqq':0.6945,'24msE':0.69512,'24msgg':0.6908,'24msqq':0.7060,'32MOME':0.57773,'32MOMgg':0.6585,'32MOMqq':0.6940,'32msE':0.69464,'32msgg':0.6931,'32msqq':0.7056,'48MOME':0.57764,'48MOMgg':0.6563,'48MOMqq':0.6945,'48msE':0.69512,'48msgg':0.6908,'48msqq':0.7060,'64MOME':0.57773,'64MOMgg':0.6585,'64MOMqq':0.6940,'64msE':0.69464,'64msgg':0.6931,'64msqq':0.7056,'48fineMOME':1,'48fineMOMgg':0.658285,'48fineMOMqq':1,'48finemsE':1,'48finemsgg':0.685909,'48finemsqq':1,}

        #Below is Zp at 2Gev - there so I can test with the 48fine at 2GeV BUT I have no values for 48 fine.
        #dictZP = {'24MOME':0.57764,'32MOME':0.57773,'24msE':0.69512,'32msE':0.69464,'24MOMqq':0.6423,'32MOMqq':0.6940,'24msqq':0.6372,'32msqq':0.6506,'24MOMgg':0.5974,'32MOMgg':00.6585,'24msgg':0.5924,'32msgg':0.6320,'48MOME':0.57764,'64MOME':0.57773,'48fineMOME':1,'48msE':0.69512,'64msE':0.69464,'48finemsE':1,'48MOMqq':0.6423,'64MOMqq':0.6940,'48fineMOMqq':1,'48msqq':0.6372,'64msqq':0.6506,'48finemsqq':1,'48MOMgg':0.5974,'64MOMgg':00.6585,'48fineMOMgg':1,'48msgg':0.5924,'64msgg':0.6320,'48finemsgg':1,}
        #self.ZP = dictZP[self.lattice_dim+self.scheme_name+self.kin_name]
        print self.ZP
        #self.ZA = dictZA[self.lattice_dim]
'''
