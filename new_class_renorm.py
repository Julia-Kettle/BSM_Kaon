from scipy import stats
from struct import *
import numpy as np
from copy import deepcopy
from decimal import *
from utils import *
from math import pi
from global_settings import *
from linearfit import *

class Renormalisation:

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
        if ibeta > 1:
            self.lattice_name_m2beta = ['24','32'][ibeta-2]
        self.m_sea_l = [[0.0050, 0.0100, 0.0200],[0.0040, 0.0060, 0.0080],[0.00078],[0.000678]][ibeta]
        self.m_val_s = [[0.0300, 0.0350, 0.0400], [0.0250, 0.0300],[0.0362],[0.02661]][ibeta]
        self.m_sea_s_phys = [0.03224,0.02447,0.03580,0.02359][ibeta]
        self.m_l_phys = [0.00078,0.000678,0.00078,0.000678][ibeta]
        self.n_mseal = [3,3,1,1][ibeta] #check the 2nd value of this
        self.n_mval1 = [1,1,1,1][ibeta]
        self.n_mval2 = [3,2,1,1][ibeta]
        self.nboots = 500
        self.N_ren = [8.0/3,4.0/3,-2.0,5.0/3,1.0] #This is the factor in the bag parameter
        self.N = [8.0/3,-5.0/3,1.0/3,2.0,0/3]
        self.NrenInv = [3.0/8,3.0/4,-0.5,3.0/5,1.0]
        self.nchan=5
        self.T = [ [np.array([[1,0,0,0,0],[0,0,0,1,0],[0,0,0,-0.5,0.5],[0,0,1,0,0],[0,-0.5,0,0,0]])
], [np.identity(5) ] ][ibasis]
        self.Mkphys = 0.5*(493.677+497.614)*pow(10,-3)
        self.Fkphys = 156.2*pow(10,-3)

    def read_ainvs(self):
        #set up the ainv. For the physical have generated fake bootstraps using a pseudorand gaussian dist with sd corresponding to error
        filename = 'common_data/boot_ainv_'+self.lattice_name+'cubed_IW_'+str(self.nboots)
        self.bootainv , err = self.read_bootstraps(filename)
        ainv_cent = [1.7848,2.3833,1.7295,2.3586]
        self.bootainv[-1] = ainv_cent[self.ibeta]
        self.a2 = []
        for b in range(len(self.bootainv)):
            self.a2.append(pow(self.bootainv[b],-2))


    ###################################################FILEREADING#############################################
    def read_datacol_file(self,filename):
        "reads in data from a file in list format"
        fo = open(filename,'r')
        lines = fo.readlines()
        data = []
        for line in lines:
            if type(line)==str:
                line = float(line)
                data.append(line)
        return data

    def store_col_3D(self,x,y,data):
        "rearranges list of data into np array of dimension [x,y,boots]"
        col = len(data)/(x*y)
        array = np.zeros(x,y,col)
        counter=0
        for i in range(x):
            for j in range(y):
                for k in range(col):
                    array[i,j,k] = data[counter]
                    counter += 1
        return data

    def read_Z(fself,filename,nchan,nboot):
        "reads in the Z boots and puts into np array dim [nchan,nchan,nboots]"
        "could be replicated using read data_col_file & store_col_3d" 
        "Could be got rid of"
        Z = np.zeros([nchan,nchan,nboot+1])
        data = self.read_datacol_file(filename)
        counter=0
        for i in range(nchan):
            for j in range(nchan):
                for k in range(nboot+1):
                    Z[i,j,k] = data[counter]
                    counter += 1
        return Z

    def read_bootstraps(self,filename):
        '''
        bootstrap files are files with nboot+1 LINES.
        contain bootstraps + error for one variable
        '''
        3
        #do I need this? does barely anything
        data = self.read_datacol_file(filename)
        error = stats.sem(data)
        return data, error

    def read_results(self,filename,nc,n_mseal,n_mval1,n_mval2):
        #read data from file into 5D array
        data=self.read_datacol_file(filename)
        array = np.zeros([nc,n_mseal,n_mval1,n_mval2,self.nboots+1])
        counter = 0
        for j in range(n_mseal):
            for k in range(n_mval1):
                for l in range(n_mval2):
                    for i in range(nc):
                        for m in range(self.nboots+1):
                            array[i,j,k,l,m] = data[counter]
                            counter+=1
        return array


    def read_results_Jamie(self,filename):
        '''
        Need to read in binary files. consist of (4byte int) n, nboot, 
        (8byte double) nboot+1 times dummy data, (4byte int) n, nboot, 
        (8byte double) nboot+1 times real data.
        Use struct module to read in 4 or 8 bytes(depending on file posn)
        then convert back to int or float
        When converting use >, to indicate order, bigendian
        Values gained are bootstraps + central of either fk, mk or mk0fk(crossterms)
        '''
        fo = open(filename,"rb")
        sect = fo.read(4) #read 4bytes of data (1st n)
        count=0;
        dummy = []
        data = []
        while sect !="":
            count=count+1
            if count<=1: #first n
                n = unpack('>i',sect)[0]#convert from bin to int (bugendian)
                sect = fo.read(4) 
            elif count <=2: #first nboot
                nboot = unpack('>i',sect)[0]
                sect = fo.read(8)
            elif count <= nboot+2: #the dummy variables except last
                dummy.append((unpack('>d',sect))[0])
                sect = fo.read(8)
            elif count<= nboot+3: #last dummy variable - switch back to 4 byte
                dummy.append(unpack('>d',sect)[0])
                sect = fo.read(4)
            elif count <= nboot+4:#2nd n
                n = unpack('>i',sect)[0]
                sect = fo.read(4)
            elif count <= nboot+5:#2nd nboot
                nboot = unpack('>i',sect)[0]
                sect = fo.read(8)
            else: #data
                data.append(unpack('>d',sect)[0])
                sect = fo.read(8)

        return n, nboot, data

    ######################################Return the Data#########################################################
    def return_3pts(self):
        #Start reading 3pt functions
        print "Reading 3pts..."
        msea = self.m_sea_l
        ms = self.m_val_s
        latticedim = self.lattice_name
        filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p_over_2p.txt'
        self.matel_3p_2p = self.read_results(filename,5,len(msea),len(ms),1) #3pt/(2pt*2pt) ~ Bag
        filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p.txt'
        self.matel_3p= self.read_results(filename,5,len(msea),len(ms),1) #3pt
        filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p_rat.txt'
        self.matel_3p_rat = self.read_results(filename,4,len(msea),len(ms),1) #3pt1/3pti ~ R
        print "Done" + "\n"


    def return_myBag(self):
        print " Reading Physical point Bag data..."
        ms  = deepcopy(self.m_val_s)
        msea = deepcopy(self.m_sea_l)
        kaon = 's' + str(ms[0]) + '-l'+ str(msea[0])
        lat = self.lattice_name+"cubed"
        Ninv = [3.0/8,3.0/4,-0.5,3.0/5,1.0] #removes scaling to be consistent with Nicolas, re-added in later.
        self.matel_3p_2p = np.zeros([self.nchan,1,1,1,self.nboots+1])
        self.matel_3p_rat = np.zeros([self.nchan-1,1,1,1,self.nboots+1])    
        filename=[]
        if self.ibeta == 2:
            dt = 40
        elif self.ibeta == 3:
            dt=40
        #loop throught the channels
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
            self.matel_3p_2p[ic] = self.read_results(filenameB,1,1,1,1)
            if ic > 0:
                self.matel_3p_rat[ic-1] = self.read_results(filenameR,1,1,1,1)
            dataR = self.read_bootstraps(filenameR)
            for iboot in range(self.nboots+1):
                self.matel_3p_2p[ic,0,0,0,iboot] = self.matel_3p_2p[ic,0,0,0,iboot]*self.N_ren[ic]
        print "Done"+"\n"


    def return_myMeson(self):
        print " Reading Physical point Mass data..."
        #setup masses & lat
        ms = deepcopy(self.m_val_s)
        msea = deepcopy(self.m_sea_l)
        lat = self.lattice_name+"cubed"
        mval=np.hstack([msea,ms]) # list of light,  strange val quarks
        self.b_M=np.zeros([1,1,2,1,self.nboots+1])
        self.b_F=np.zeros([1,1,2,1,self.nboots+1])
        self.b_msq_fsq=np.zeros([1,1,2,1,self.nboots+1])
        for imv2 in range(2):
            if imv2 == 1:
                meson = 'l' + str(mval[0]) + '-l' + str(msea[0])
            elif imv2 == 0:
                meson = 's' + str(mval[1]) + '-l'+ str(msea[0])
            #set up matrices
            #set up filenames & read data
            filenameF='../Fits/' + lat + '/mass/mass-' + meson + '/decay-' + meson +'_boots.dat'
            filenameM='../Fits/' + lat + '/mass/mass-' + meson + '/mass-' + meson +'_boots.dat'
            self.b_M[0,0,imv2,0,:] = self.read_results(filenameM,1,1,1,1)
            self.b_F[0,0,imv2,0,:] = self.read_results(filenameF,1,1,1,1)
            for iboot in range(self.nboots+1):
                self.b_msq_fsq[0,0,imv2,0,iboot] = pow(self.b_M[0,0,imv2,0,iboot]/(4*pi*self.b_F[0,0,imv2,0,iboot]),2)
            print "Done" + "\n"

    def return_Jamie_data(self):
        #Start reading Jamie's data (2pts)
        print "Reading 2pts..."
        #get quark masses + physical
        msea = deepcopy(self.m_sea_l) #mv1=msea always
        ms = deepcopy(self.m_val_s)
        ms_phys = deepcopy(self.m_sea_s_phys)
        filend_basis = ['_9.22', '_12.52']

        #set up matrices to store results
        self.b_F = np.zeros([1,len(msea),len(ms)+1,1,self.nboots+1])
        self.b_M = np.zeros([1,len(msea),len(ms)+1,1,self.nboots+1])
        self.b_msq_fsq = np.zeros([1,len(msea),len(ms)+1,1,self.nboots+1])

        for i in range(len(msea)):
            mv1 = [deepcopy(msea[i])] #1st valence l + sea

            for k in range(len(mv1)):
                    mv2 = np.hstack((ms,msea[i])) # 2nd valence, either s or l
                    for j in range(len(mv2)):

                        filenameF = 'data_Jamie/'+self.lattice_name + '/FK_4CHAN_mu' +str5sf(msea[i])+'_mq'+str5sf(mv1[k])+'_ms'+str5sf(mv2[j]) + filend_basis[self.ibeta]
                        filenameM = 'data_Jamie/'+self.lattice_name + '/MK_4CHAN_mu' +str5sf(msea[i])+'_mq'+str5sf(mv1[k])+'_ms'+str5sf(mv2[j]) + filend_basis[self.ibeta]
                        n, nboot, dataF = self.read_results_Jamie(filenameF)
                        n, nboot, dataM = self.read_results_Jamie(filenameM)
                        for l in range(len(dataM)):
                                self.b_F[0,i,j,k,l] = dataF[l]
                                self.b_M[0,i,j,k,l] = dataM[l]
                                self.b_msq_fsq[0,i,j,k,l]= pow((self.b_M[0,i,j,k,l]/(4*pi*self.b_F[0,i,j,k,l])),2)
        print "Done"+"\n"



    def calcB(self):
        dim = np.shape(self.matel_3p_2p)
        self.B = np.zeros(dim)
        for i in range(dim[0]):
            for j in range(dim[1]):
                    for k in range(dim[2]):
                        for l in range(dim[3]):
                                for m in range(dim[4]):
                                    self.B[i,j,k,l,m] = self.matel_3p_2p[i,j,k,l,m]/self.N_ren[i]


    def calcR(self):
        """
        Calculates the ratio of each channel 3pt over sm 3pt.
        timesed by a factor (Fk/Mk)^2phys*(Mpi/Fpi)^2
        """
        dim = []
        dim1 = np.shape(self.matel_3p_rat)
        dim2 = np.shape(self.b_F)
        for i in range(len(dim1)):
            dim.append(min(dim1[i],dim2[i]))
        dim[0] = self.nchan
        self.R = np.zeros(dim)
        #Mkphys = 139.57018
        #Fkphys=130
        for i in range(dim[1]):
            for j in range(dim[2]):
                    for k in range(dim[3]):
                        for m in range(dim[4]):
                                for l in range(1,dim[0]):
                                    self.R[0,i,j,k,m] = pow((self.Fkphys/self.Mkphys),2)*pow((self.b_M[0,i,j,k,m]/self.b_F[0,i,j,k,m]),2)
                                    self.R[l,i,j,k,m] = pow((self.Fkphys/self.Mkphys),2)*pow((self.b_M[0,i,j,k,m]/self.b_F[0,i,j,k,m]),2)*self.matel_3p_rat[l-1,i,j,k,m]


    def Read_Z(self,filename):
        #Read the Z factors into 5by5bynbootsmatrix.
        #But bootstraps on 
        nboot=100
        fi = open(filename,'r')
        lines = fi.readlines()
        self.bZ_scheme=np.zeros([self.nchan,self.nchan,nboot+1])
        counter = 0
        #Z_cent=np.zeros([self.nchan,self.nchan])
        #Z_err=np.zeros([self.nchan,self.nchan])
        for ich in range(self.nchan):
            for jch in range(self.nchan):
                for iboot in range(nboot+1):
                    self.bZ_scheme[ich][jch][iboot] = lines[counter]
                    counter += 1
                #Z_cent[ich][jch]=Z_boots[ich][jch][-1]
                #Z_err=np.std(Z_boots,axis=2)




    def ZA_ZP(self):
        if self.ibeta > 1:
            iBeta = self.ibeta-2
        else:
            iBeta  = self.ibeta
        self.ZA=[0.71273,0.74404][iBeta]
        self.ZP=[[0.57732,0.57733],[0.69512,0.69464],[0.65635,0.65832],[0.65635*1.0499,0.65832*1.0499],[0.69513,0.69471],[0.69513*1.0157,0.69471*1.0157],[0.69513,0.69471],[0.65635,0.65832],[0.65635*1.0499,0.65832*1.0499]][self.ikin][iBeta]


    def Renormalise(self,doR,doB):
        #Pass the filename with the Zfactors, basis, kinematics and vol to renormalise the ratio and the bag
        self.ZA_ZP() 
        if self.ibeta < 2:
            Zfilename = "data/Z" + self.scheme_name + "_boot_mu_match_" + self.lattice_name + "cubed_" + self.kin_name + "_block.out"
        else:
            Zfilename = "data/Z" + self.scheme_name + "_boot_mu_match_" + self.lattice_name_m2beta + "cubed_" + self.kin_name + "_block.out"
        self.Read_Z(Zfilename)
        self.N=np.squeeze(np.dot(self.T,self.N_ren))
        if doB == 1:
            self.MultBZ()
        if doR ==1:
            self.MultRZ()


    def MultRZ(self):
        #Scale the ratios R by Z factor
        dimR = np.shape(self.R)
        self.ren_R = np.zeros(dimR)
        self.ren_phys_R = np.zeros(dimR)
        for iml in range(dimR[1]):
            for ims in range(dimR[2]):
                    for imq in range(dimR[3]):
                        for iboot in range(dimR[4]):
                                self.ren_R[:,iml,ims,imq,iboot] = (np.dot(self.bZ_scheme[:,:,-1],self.R[:,iml,ims,imq,iboot]))/self.bZ_scheme[0][0][-1]
                                self.ren_phys_R[:,iml,ims,imq,iboot] = np.dot(self.T,self.ren_R[:,iml,ims,imq,iboot])


    def MultBZ(self):
        #Scale the bag B by the Z factors
        dimB = np.shape(self.B)
        self.ren_B = np.zeros(dimB)
        self.ren_phys_B = np.zeros(dimB)
        self.B_phys = np.zeros(dimB)
        for iml in range(dimB[1]):
            for ims in range(dimB[2]):
                    for imq in range(dimB[3]):
                        for iboot in range(dimB[4]):
                                self.ren_B[:,iml,ims,imq,iboot] = (self.ZA*self.ZA)*(np.dot(self.bZ_scheme[:,:,-1],self.B[:,iml,ims,imq,iboot]))
                                self.ren_phys_B[:,iml,ims,imq,iboot] = np.dot(self.T,self.ren_B[:,iml,ims,imq,iboot])
                                for ich in range(5):
                                    self.B_phys[ich,iml,ims,imq,iboot] = self.ren_phys_B[ich,iml,ims,imq,iboot]/(self.ZP*self.ZP*self.N[ich])
                                self.B_phys[0,iml,ims,imq,iboot] = self.ren_phys_B[0,iml,ims,imq,iboot]/(self.ZA*self.ZA*self.N[0])




    #this is our main function
    #we want to pass the class/function (do we need a class?) Incorporate into the renormalisation??
    #pass R/B - loop through the channels,lights,mval2, to get a matrix [strange,nboot]
    #then loop through the boots interpolating

    def data_Interp(self, data):
        print "-----Interpolating non-physical mass data-----"
        dim = np.shape(data)
        #shape channels, sea, mval2, boot
        Iy = np.zeros([dim[0],dim[1],dim[3],dim[4]])
        store_dydms=np.zeros([dim[0],dim[1],dim[3],dim[4]])
        for ich in range(dim[0]): #channel
            for isea in range(dim[1]): #sea mass
                for iv2 in range(dim[3]): #valence 2
                    yerr=np.std(data[ich,isea],axis=2)
                    yerr = yerr[:,0]
                    store_params=[]
                    for ib in range(dim[4]): #boot
                        ydata = []
                        for ivs in range(dim[2]):
                            #strange val - #last entry is light
                            ydata.append(data[ich,isea,ivs,iv2,ib])
                        Iy[ich,isea,iv2,ib], store_dydms[ich,isea,iv2,ib] = Interp(self.m_val_s,ydata,yerr,self.m_sea_s_phys,1)
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
        Idy = np.zeros([dimdy[0],dimdy[2],dimdy[3]])
        for ich in range(dimdy[0]): #channel
            for iv2 in range(dimdy[2]): #valence 2
                yerr1=np.std(dydms[ich,:,:],axis=2)
                for ib in range(dimdy[3]): #boot
                    ydata = []
                    yerr = []
                    for il in range(dimdy[1]):	
                        #strange val - #last entry is light
                        ydata.append(dydms[ich,il,iv2,ib])
                        yerr.append(yerr1[il,0])
                    Idy[ich,iv2,ib], params=Interp(self.m_sea_l,ydata,yerr,self.m_l_phys,1)
        return Idy

    def phys_Interp(self,y_sim,dydms):
        print "Interpolating physical pt data to true strange physical point"
        dimy = np.shape(y_sim) #(4/5,1,1,1,501)
        y_p = np.zeros([dimy[0],dimy[1],dimy[2],dimy[4]])
        for ich in range(dimy[0]):
            for isea in range(dimy[1]):
                for iv2 in range(dimy[2]):
                    for istr in range(dimy[3]):
                        for iboot in range(dimy[4]):
                            y_p[ich,isea,iv2,iboot]=y_sim[ich,isea,iv2,istr,iboot] + dydms[ich,iv2,iboot]*(self.m_sea_s_phys - np.squeeze(self.m_val_s))
                    #print "x_sim:  {}\ty_sim:  {}\tdy/dx:  {}\txp:  {}\typ:  {} +/- {}".format(x_sim,y_sim[ich,isea,iv2,istr,-1],dRdms[ich,iv2,-1],x_p,y_p[ich,isea,iv2,-1],np.std(y_p[ich,isea,iv2,:500]))
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

    def Run_All(self):
        self.ReadAllDataIn()
        self.Renormalise(1,1)
        self.read_ainvs()
        #self.IM = self.data_Interp(b_M)
        #self.IF = self.data_Interp(b_F)
        return self.B_phys, self.ren_phys_R, self.b_msq_fsq[:,:,:,:,:], self.b_M, self.b_F, self.a2

    def Set_IR(self,IR):
        self.IR = IR
    
    def Set_IB(self,IB):
        self.IB = IB



