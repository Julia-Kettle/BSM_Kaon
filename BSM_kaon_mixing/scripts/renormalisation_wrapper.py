from scipy import stats
from struct import *
import numpy as np
from copy import deepcopy
from decimal import *
from linearfit import *
from math import pi
import utils
from .renormalisation import Renormalisation

class RenormalisationWrapper(Renormalisation):

    def __init__(self,ibeta,ischeme,ikin,ibasis):
        "Initialise the global variables, and those specific to this lattice"
        self.ibeta = ibeta
        self.ikin = ikin
        self.ischeme = ischeme
        self.ibasis = ibasis
        self.scheme_name = ['MOM','ms'][ischeme]
        self.basis_name = ['SUSY','Lattice'][ibasis]
        self.kin_name = ['E','qq','qg','gq','gg'][ikin]
        self.lattice_name = ['24','32','48','64'][ibeta]
        self.lattice_name_Zs = ['24','32','24','32'][ibeta]
        self.m_sea_l = [[0.0050, 0.0100, 0.0200],[0.0040, 0.0060, 0.0080],[0.00078],[0.000678]][ibeta]
        self.m_val_s = [[0.0300, 0.0350, 0.0400], [0.0250, 0.0300],[0.0362],[0.02661]][ibeta]
        self.m_sea_s_phys = [0.03224,0.02447,0.03580,0.02359][ibeta]
        self.m_l_phys = [0.00078,0.000678,0.00078,0.000678][ibeta]#Is this shortcut or only used to interpolate the slope
        self.n_mseal = [3,3,1,1][ibeta] #check the 2nd value of this
        self.n_mval1 = [1,1,1,1][ibeta]
        self.n_mval2 = [3,2,1,1][ibeta]
        self.nboots = 500
        self.N_ren = [8.0/3,4.0/3,-2.0,5.0/3,1.0] #This is the factor in the bag parameter
        self.N = [8.0/3,-5.0/3,1.0/3,2.0,2.0/3]
        self.nchan=5
        self.T = [ [np.array([[1,0,0,0,0],[0,0,0,1,0],[0,0,0,-0.5,0.5],[0,0,1,0,0],[0,-0.5,0,0,0]])
], [np.identity(5) ] ][ibasis]
        self.Mkphys = 0.5*(493.677+497.614)*pow(10,-3)
        self.Fkphys = 156.2*pow(10,-3)
        self.nboots = 500
    def read_ainvs(self):
        #set up the ainv. For the physical have generated fake bootstraps using a pseudorand gaussian dist with sd corresponding to error
        filename = 'common_data/boot_ainv_'+self.lattice_name+'cubed_IW_'+str(self.nboots)
        self.bootainv = utils.Read_in().read_file_list(filename)
        ainv_cent = [1.7848,2.3833,1.7295,2.3586]
        self.bootainv[-1] = ainv_cent[self.ibeta]
        self.a2 = []
        for b in range(len(self.bootainv)):
            self.a2.append(pow(self.bootainv[b],-2))




    ######################################Return the Data#########################################################
    def return_3pts(self):
        """
        Reads in the 3 point functions into suitably sized arrays
        """
        msea = self.m_sea_l
        ms = self.m_val_s
        latticedim = self.lattice_name

        filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p_over_2p.txt'
        self.matel_3p_2p = np.zeros([5,len(msea),len(ms),self.nboots+1])
        utils.Read_in().read_file_array(filename,self.matel_3p_2p) #3pt/(2pt*2pt) ~ Bag

        filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p.txt'
        self.matel_3p=np.zeros([5,len(msea),len(ms),self.nboots+1])
        utils.Read_in().read_file_array(filename,self.matel_3p) #3pt

        filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p_rat.txt'
        self.matel_3p_rat = np.zeros([4,len(msea),len(ms),self.nboots+1])
        utils.Read_in().read_file_array(filename,self.matel_3p_rat) #3pt1/3pti ~ R

        print "Unrenormalised bag is - ", self.matel_3p_2p[:,0,0,-1]
        print "Unrenormalised bare ratios is - ", self.matel_3p_rat[:,0,0-1]
        #print "Done" + "\n"


    def return_myBag(self):
        """
        Reads in the new physical point bags and bare ratios
        """
        ms  = deepcopy(self.m_val_s)
        msea = deepcopy(self.m_sea_l)
        kaon = 's' + str(ms[0]) + '-l'+ str(msea[0])
        lat = self.lattice_name+"cubed"
        Ninv = [3.0/8,3.0/4,-0.5,3.0/5,1.0] #removes scaling to be consistent with Nicolas, re-added in later.

        self.matel_3p_2p = np.zeros([self.nchan,1,1,self.nboots+1])
        self.matel_3p_rat = np.zeros([self.nchan-1,1,1,self.nboots+1])    

        if self.ibeta == 2:
            dt = 40
        elif self.ibeta == 3:
            dt=52

        for ic in range(self.nchan):

            # 3 & 4 need to switched - convention choice
            if ic == 2:
                filenameR='../Fits/' + lat  + '/3ptrat/channel' +str(2*ic+2) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
                filenameB='../Fits/' + lat  + '/bag/channel' +str(2*ic+2) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
            elif ic ==3:
                filenameR='../Fits/' + lat  + '/3ptrat/channel' +str(2*ic-2) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
                filenameB='../Fits/' + lat  + '/bag/channel' +str(2*ic-2) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
            else:
                filenameR='../Fits/' + lat  + '/3ptrat/channel' +str(2*ic) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
                filenameB='../Fits/' + lat  + '/bag/channel' +str(2*ic) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'

            utils.Read_in().read_file_array(filenameB,self.matel_3p_2p[ic])
            if ic > 0:
                utils.Read_in().read_file_array(filenameR,self.matel_3p_rat[ic-1])
            dataR = utils.Read_in().read_file_list(filenameR)

            for iboot in range(self.nboots+1):
                self.matel_3p_2p[ic,0,0,iboot] = self.matel_3p_2p[ic,0,0,iboot]*self.N_ren[ic]


    def return_myMeson(self):
        #print " Reading Physical point Mass data..."
        #setup masses & lat
        ms = deepcopy(self.m_val_s)
        msea = deepcopy(self.m_sea_l)
        lat = self.lattice_name+"cubed"
        mval=np.hstack([msea,ms]) # list of light,  strange val quarks
        self.b_M=np.zeros([1,1,2,self.nboots+1])
        self.b_F=np.zeros([1,1,2,self.nboots+1])
        self.b_msq_fsq=np.zeros([1,1,2,self.nboots+1])
        for imv2 in range(2):
            if imv2 == 1:
                meson = 'l' + str(mval[0]) + '-l' + str(msea[0])
            elif imv2 == 0:
                meson = 's' + str(mval[1]) + '-l'+ str(msea[0])
            #set up matrices
            #set up filenames & read data
            filenameF='../Fits/' + lat + '/mass/mass-' + meson + '/decay-' + meson +'_boots.dat'
            filenameM='../Fits/' + lat + '/mass/mass-' + meson + '/mass-' + meson +'_boots.dat'
            utils.Read_in().read_file_array(filenameM,self.b_M[0,0,imv2,:])
            utils.Read_in().read_file_array(filenameF,self.b_F[0,0,imv2,:])
            for iboot in range(self.nboots+1):
                self.b_msq_fsq[0,0,imv2,iboot] = pow(self.b_M[0,0,imv2,iboot]/(4*pi*self.b_F[0,0,imv2,iboot]),2)
            #print "Done" + "\n"

    def return_Jamie_data(self):
        #Start reading Jamie's data (2pts)
        #print "Reading 2pts..."
        #get quark masses + physical
        msea = deepcopy(self.m_sea_l) #mv1=msea always
        ms = deepcopy(self.m_val_s)
        ms_phys = deepcopy(self.m_sea_s_phys)
        filend_basis = ['_9.22', '_12.52']

        #set up matrices to store results
        self.b_F = np.zeros([1,len(msea),len(ms)+1,self.nboots+1])
        self.b_M = np.zeros([1,len(msea),len(ms)+1,self.nboots+1])
        self.b_msq_fsq = np.zeros([1,len(msea),len(ms)+1,self.nboots+1])

        for i in range(len(msea)):
                mv2 = np.hstack((ms,msea[i])) # 2nd valence, either s or l
                for j in range(len(mv2)):
                    filenameF = 'data_Jamie/'+self.lattice_name + '/FK_4CHAN_mu' +utils.Utilities().str5sf(msea[i])+'_mq'+utils.Utilities().str5sf(msea[i])+'_ms'+utils.Utilities().str5sf(mv2[j]) + filend_basis[self.ibeta]
                    filenameM = 'data_Jamie/'+self.lattice_name + '/MK_4CHAN_mu' +utils.Utilities().str5sf(msea[i])+'_mq'+utils.Utilities().str5sf(msea[i])+'_ms'+utils.Utilities().str5sf(mv2[j]) + filend_basis[self.ibeta]
                    n, nboot, dataF = utils.Read_in().read_results_Jamie(filenameF)
                    n, nboot, dataM = utils.Read_in().read_results_Jamie(filenameM)
                    for l in range(len(dataM)):
                            self.b_F[0,i,j,l] = dataF[l]
                            self.b_M[0,i,j,l] = dataM[l]
                            self.b_msq_fsq[0,i,j,l]= pow((self.b_M[0,i,j,l]/(4*pi*self.b_F[0,i,j,l])),2)
        #print "Done"+"\n"


    def calcR(self):
        dim = []
        dim1 = np.shape(self.matel_3p_rat)
        dim2 = np.shape(self.b_F)
        for i in range(len(dim1)):
            dim.append(min(dim1[i],dim2[i]))
        dim[0] = self.nchan
        self.R = np.zeros(dim)
        #Mkphys = 139.57018
        #Fkphys=130
        for index, val in np.ndenumerate(self.R):
            if index[0] == 0:
                self.R[index] = pow((self.Fkphys/self.Mkphys),2)*pow((self.b_M[index])/self.b_F[index],2)
            else:
                self.R[index] = pow((self.Fkphys/self.Mkphys),2)*pow((self.b_M[(0,)+index[1:]])/self.b_F[(0,)+index[1:]],2)*self.matel_3p_rat[(index[0]-1,) + index[1:]]

    def Read_Z(self,filename):
        nboot=100
        self.bZ_scheme=np.zeros([self.nchan,self.nchan,nboot+1])
        utils.Read_in().read_file_array(filename,self.bZ_scheme)



    def ZA_ZP(self):
        if self.ibeta > 1:
            iBeta = self.ibeta-2
        else:
            iBeta  = self.ibeta
        self.ZA=[0.71273,0.74404][iBeta]#,0.71076,0.74293][iBeta]
        self.ZP=[[0.57732,0.57733],[0.69512,0.69464],[0.6563,0.6585],[0.6563*1.05259,0.6585*1.05259],[0.6945,0.6940],[0.6945*1.01664,0.6940*1.01664],[0.6945,0.6940],[0.6563,0.6585],[0.6563*1.05259,0.6585*1.05259]][self.ikin][iBeta]
        print "ZA - ",self.ZA,"   ZP - ", self.ZP
    '''
    def ZA_ZP(self):
        if self.ibeta > 1:
            iBeta = self.ibeta-2
        else:
            iBeta  = self.ibeta
        self.ZA=[0.71273,0.74404][iBeta]
        self.ZP=[[0.57732,0.57733],[0.69512,0.69464],[0.6563,0.6585],[0.6563*1.05259,0.6585*1.05259],[0.6945,0.6940],[0.6945*1.01664,0.6940*1.01664],[0.6945,0.6940],[0.6563,0.6585],[0.6563*1.05259,0.6585*1.05259]][self.ikin][iBeta]
    ''' 




    #this is our main function
    #we want to pass the class/function (do we need a class?) Incorporate into the renormalisation??
    #pass R/B - loop through the channels,lights,mval2, to get a matrix [strange,nboot]
    #then loop through the boots interpolating

    def data_Interp(self, data):
        print "-----Interpolating non-physical mass data-----"
        dim = np.shape(data)
        #shape channels, sea, mval2, boot
        Iy = np.zeros([dim[0],dim[1],dim[3]])
        store_dydms=np.zeros([dim[0],dim[1],dim[3]])
        for ich in range(dim[0]): #channel
            for isea in range(dim[1]): #sea mass
                yerr=np.std(data[ich,isea],axis=1)
                store_params=[]
                for ib in range(dim[3]): #boot
                    ydata = []
                    for ivs in range(dim[2]):
                        #strange val - #last entry is light
                        ydata.append(data[ich,isea,ivs,ib])
                    Iy[ich,isea,ib], store_dydms[ich,isea,ib] = Interp(self.m_val_s,ydata,yerr,self.m_sea_s_phys,1)
                    #store_dydms[ich,isea,iv2,ib] = params
                #print "channel -",ich, "\t ml -",  m_l[isea]
                #for ivs in range(dim[2]):
                #    print "x({}):  {}\ty:   {} +/- {}".format(ivs,m_s[ivs],R[ich,isea,ivs,iv2,-1],np.std(R[ich,isea,ivs,iv2,:500]))
                #print "x(p):  {}\typ:  {} +/- {}".format(m_s_phys,IR[ich,isea,iv2,-1],np.std(IR[ich,isea,iv2,:500]))
                #print "dR/dms:  {} +/- {}\n".format(store_dRdms[ich,isea,iv2,-1],np.std(store_dRdms[ich,isea,iv2,:500]))
        return Iy, store_dydms



    def slope_Interp(self,dydms):
        print "Interpolate dR/dms to the near physical point data"
        #shape channels, sea, mval2, boot
        dimdy=np.shape(dydms)
        print self.lattice_name, " - dy/dms = ", dydms[0,:,-1] 
        Idy = np.zeros([dimdy[0],dimdy[2]])
        for ich in range(dimdy[0]): #channel
            yerr1=np.std(dydms[ich,:])
            for ib in range(dimdy[2]): #boot
                ydata = []
                yerr = []
                for il in range(dimdy[1]):	
                    #strange val - #last entry is light
                    ydata.append(dydms[ich,il,ib])
                    yerr.append(yerr1)
                Idy[ich,ib], params=Interp(self.m_sea_l,ydata,yerr,self.m_l_phys,1)
        return Idy

    def phys_Interp(self,y_sim,dydms):
        print "Interpolating physical pt data to true strange physical point"
        dimy = np.shape(y_sim) #(4/5,1,1,1,501)
        y_p = np.zeros([dimy[0],dimy[1],dimy[3]])
        print "uninterpolated value is ", y_sim[0,:,-1]
        for ich in range(dimy[0]):
            for isea in range(dimy[1]):
                    for istr in range(dimy[2]):
                        for iboot in range(dimy[3]):
                            y_p[ich,isea,iboot]=y_sim[ich,isea,istr,iboot] + dydms[ich,iboot]*(self.m_sea_s_phys - np.squeeze(self.m_val_s))
                    #print "x_sim:  {}\ty_sim:  {}\tdy/dx:  {}\txp:  {}\typ:  {} +/- {}".format(x_sim,y_sim[ich,isea,iv2,istr,-1],dRdms[ich,iv2,-1],x_p,y_p[ich,isea,iv2,-1],np.std(y_p[ich,isea,iv2,:500]))

        print "interpolated physical value is", y_p[0,:,-1]
        return y_p




    def ReadAllDataIn(self):
        if self.ibeta < 2:
            self.return_Jamie_data()
            self.return_3pts()

        else:
            self.return_myBag()
            self.return_myMeson()
        self.calcR()
        self.B = self.matel_3p_2p

    def Return_All_Renormalised_Data(self):
        return self.B_phys, self.ren_phys_R
    """
    def CalculateMatrixElem(self,IB,IR):
        #optional function to calculate the matrix elements back
        self.IB=IB
        self.IR=IR
        dimB = np.shape(IB)

        matel = np.zeros(dimB)
        matel_alt = np.zeros([dimB[0]-1,dimB[1],dimB[2],dimB[3]])
        mk=0.4956
        fk=0.1562
        md=0.003162
        ms=0.08735
        for ic in range(5):
            for isea in range(dimB[1]):
                for iv2 in range(dimB[2]):
                    for iboot in range(dimB[3]):
                        if ic == 0:
                            matel[ic,isea,iv2,iboot] = self.N[ic]*IB[ic,isea,iv2,iboot]*mk*fk*mk*fk
                        else:
                            matel[ic,isea,iv2,iboot] = self.N[ic]*IB[ic,isea,iv2,iboot]*pow(mk,4)*pow(fk,2)/pow((md +ms),2)
                            matel_alt[ic-1,isea,iv2,iboot] = IR[ic,isea,iv2,iboot]*matel[0,isea,iv2,iboot]
        return matel, matel_alt
    """
    def Run_All(self):
        self.ReadAllDataIn()
        self.Renormalise(1,1)
        self.read_ainvs()
        return self.B_phys, self.ren_phys_R, self.b_msq_fsq, self.b_M, self.b_F, self.a2

    def Set_IR(self,IR):
        self.IR = IR
    
    def Set_IB(self,IB):
        self.IB = IB



