'''
        
Perform the Global Fit. 
Read Jamie's 2pt data for 24, 32
Read Nicolas' 3pt data for 24, 32
Calculate ratios
Renormalise
Interpolate to physical strange
Perform global fit against pion mass/decay sqr ratio and lattice spacing


Need to include the new physical point data too. 
'''

from random import *
import nonlinearfit as nlf
import linearfit as lf
import numpy as np
import random
from readfiles import *
import resampling as rsamp
import setup
from utils import *
from copy import deepcopy
from decimal import *
import scipy as sp
from calculations import *
from math import *
from renormalisation import *
import sys
from global_settings import *

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("logfile.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    

def Renormalisation(iBeta,iParams,nboot,scheme,kin):

        b_FK, b_MK, b_mpsq_fpsq = return_Jamie_data(iBeta,iParams,nboot) #read and return meson mass + decay
        matel_3p_2p, matel_3p, matel_3p_rat = return_3pts(iBeta,iParams,nboot) #read and return 3pt functions

        ##########################################
        #Calculate the Bag paramater and Ratio R 
        print "Calculate bag and R ..."
        B = matel_3p_2p
        R = calc_R(b_FK,b_MK,matel_3p_rat)
        ##########################################


        ##########################################################################################################################
        #Renormlise B and R
        Zfilename = "data/Z" + name_scheme[scheme] + "_boot_mu_match_" + volname[iBeta] + "_" + name_kin[kin] + "_block.out"
        print Zfilename
        Rren_phys,Bren_phys = Renormalise(ibas,kin,iBeta,Zfilename,R=R,B=B)
        ##########################################################################################################################

        ###########################################################
        #Interpolate results to physical masses

        msea = deepcopy(iParams.m_sea_l)
        ms = deepcopy(iParams.m_val_s)
        ms_phys = deepcopy(iParams.m_sea_s_phys)

        IR,IdR = lf.R_Interp(Rren_phys[:,:,:,:,:], ms[:], ms_phys,msea)
        IB,IdB = lf.R_Interp(Bren_phys[:,:,:,:,:], ms[:], ms_phys,msea)
        ###########################################################
 
        dimB = np.shape(IB)
        for ich in range(5): 
            OutfilenameB = "data/renormalisation/"+volname[iBeta]+"/B" + str(ich+1) + "_boot_" + volname[iBeta] + "_" + name_kin[kin] + ".dat"
            OutfilenameR = "data/renormalisation/"+volname[iBeta]+"/R" + str(ich+1) + "_boot_" + volname[iBeta] + "_" + name_kin[kin] + ".dat"
            foB = open(OutfilenameB,'w')
            foR = open(OutfilenameR,'w')
            for iml in range(dimB[1]):
                for iboot in range(501):
                    foB.write(str(IB[ich,iml,0,iboot]))
                    foB.write("\n")
                    foR.write(str(IR[ich,iml,0,iboot]))
                    foR.write("\n")

        return IR,IB,IdR,IdB

def main(scheme,kin,bas):



    params = [] #holds class with all data for each lattice.
    b_ainv=[]
    b_ZV=[]
    nboot = 500

    #Initilasie lists which hold all lattices data
    IRlst=[]
    IBlst=[]
    IdRlst=[]
    IdBlst=[]
    Rlst = []
    Blst = []
    mlst=[]
    mllst=[]
    flst = []
    m2_f2lst=[]
    ainv=[]

    ###################
    # set up kinematics 
    #eventually will need to loop through kinematics - put whole file within function, passing these values??
    #kin=0
    #ikin=0
    #scheme=0
    #ibas=0
    ###################

    for iBeta in range(4): #loop through each lattice
        print "----------------------setting up as----------------------------"
        #read lattice inverse - again where do we use this?? Currently not using bootstraps for ainve. Should I be?
        filename = 'common_data/boot_ainv_'+latticedim[iBeta]+'cubed_IW_'+str(nboot)
        bootainv , err = read_bootstraps(filename)
        b_ainv.append(bootainv)
        print iBeta, latticedim[iBeta],bootainv[-1]
        #read ZV of 2pt - why is this here? Do we not do this again within renormalisation??
        #filename = 'common_data/boot_ZV_'+latticedim[iBeta]+'cubed_IW_'+str(nboot)
        #bootzv, err = rf.read_bootstraps(filename)
        #b_ZV.append(bootzv)

        mpsq_fpsq_phys = pow(0.13957018/(4*pi*0.13),2)
        mp_phys = 0.13957018
        fp_phys = 0.13
    
    
    for iBeta in range(2): #loop through 24, and 32 (Non-Physical)

        #class cotaining paramaters set up    
        iParams = setup.Params(iBeta)
        params.append(iParams)
        ainv.append(deepcopy(iParams.ainv))
        ml=deepcopy(iParams.m_sea_l)

        print '----------', iParams.latticename, '----------'

        IR,IB,IdR,IdB = Renormalisation(iBeta,iParams,nboot,scheme,kin)

        #############################
        #Append matrices of results
        #to lists
        IRlst.append(IR)
        IBlst.append(IB)
        IdRlst.append(IdR)
        IdBlst.append(IdB)
        mllst.append(ml)
        mlst.append(b_MK)
        flst.append(b_FK)
        m2_f2lst.append(b_mpsq_fpsq)
        #############################

    print "Finished reading + interpolating Nicolas' data"
    #need to update the interpolation such that slope kept the same for all the different light masses?


    for iBeta in [2,3]: #loop through the first 2 phsy data lattices (reuse Z frm 24 & 32)  	

        print '----------', iParams.latticename, '----------'
        iParams = setup.Params(iBeta)
        params.append(iParams)
        ainv.append(deepcopy(iParams.ainv))
        B_new, R_new = return_myBag(iBeta,iParams,nboot)
        b_MK, b_FK, b_mksq_fksq = return_myMeson(iBeta,iParams,nboot,'kaon')
        b_MP, b_FP, b_mpsq_fpsq = return_myMeson(iBeta,iParams,nboot,'pion')

        #read all the new physical files here		
        #still need the R data + pion mass + decays
        R=calc_R(b_FK,b_MK,R_new)


        ##########################################################################################################################
        #Renormlise B and R
        Zfilename = "data/Z" + name_scheme[scheme] + "_boot_mu_match_" + volname[iBeta-2] + "_" + name_kin[kin] + "_block.out" 
        Rren_phys, Bren_phys = Renormalise(ibas,kin,iBeta,Zfilename,R=R,B=B_new)
        ##########################################################################################################################


        dR = IdRlst[iBeta-2]
        dB = IdBlst[iBeta-2]
        ml_phys_pt = deepcopy(iParams.m_sea_l)
        ml = mllst[iBeta-2]
        ms_phys_pt = deepcopy(iParams.m_val_s)
        ms_phys = deepcopy(iParams.m_sea_s_phys)

        IdR,dummy = lf.dR_Interp(dR,ml,ml_phys_pt[0])
        IRp=lf.phys_Interp(Rren_phys,IdR,ms_phys,ms_phys_pt[0])

        IdB,dummy = lf.dR_Interp(dB,ml,ml_phys_pt[0])
        IBp=lf.phys_Interp(Bren_phys,IdB,ms_phys,ms_phys_pt[0])

        
        #we have a multidimensional matrix need to save in a file how,
        #we have channel, ms, ml, mv2,and iboot, Just save in a list? As long as consistent 
        dimB = np.shape(IBp)
        for ich in range(5): 
            OutfilenameB = "data/renormalisation/"+volname[iBeta]+"/B" + str(ich+1) + "_boot_" + volname[iBeta] + "_" + name_kin[kin] + ".dat"
            OutfilenameR = "data/renormalisation/"+volname[iBeta]+"/R" + str(ich+1) + "_boot_" + volname[iBeta] + "_" + name_kin[kin] + ".dat"
            foB = open(OutfilenameB,'w')
            foR = open(OutfilenameR,'w')
            for iml in range(dimB[1]):
                for iboot in range(501):
                    foB.write(str(IBp[ich,iml,0,iboot]))
                    foB.write("\n")
                    foR.write(str(IRp[ich,iml,0,iboot]))
                    foR.write("\n")


        Blst.append(Bren_phys)
        Rlst.append(R)
        IRlst.append(IRp)
        #IBlst.append(B_new)
        IBlst.append(IBp)
        mlst.append(b_MP)
        flst.append(b_FP)
        m2_f2lst.append(b_mpsq_fpsq)

    p0 = [1,1,1]



    #p, cov, bp, bcov, bchisq = nlf.bootglobalfit(y1[0,0,0],y2[0,0,0],1/ainv[0],1/ainv[1],m1[0,0,0],m2[0,0,0],500,p0,mpsq_fpsq_phys)
    #actually do the global fit (looping through the sea masses)


    dimIR = np.shape(IR)
    dimIB = np.shape(IB)
    Rphys=[]
    Bphys=[]

    #a_store = np.zeros([2])
    #a_store[0] = 1/(ainv[2]*ainv[2])
    #a_store[1] = 1/(ainv[3]*ainv[3])
    #a_store = np.zeros([2,501])
    

    for ic in range(dimIB[0]): #loop channels
        IB3_store = np.zeros([2,nboot+1])
        m3_store = np.zeros([2,nboot+1])
        IR3_store = np.zeros([2,nboot+1])
        for iboot in range(nboot+1):
            IB3_store[0,iboot] = (IBlst[2])[ic,0,0,iboot]
            IB3_store[1,iboot] = (IBlst[3])[ic,0,0,iboot]
            m3_store[0,iboot] = (m2_f2lst[2])[0,0,0,0,iboot]
            m3_store[1,iboot] = (m2_f2lst[3])[0,0,0,0,iboot]
            IR3_store[0,iboot] = (IRlst[2])[ic,0,0,iboot]
            IR3_store[1,iboot] = (IRlst[3])[ic,0,0,iboot]
        for iv2 in range(dimIB[2]): #loop valence 2 (=1 here, unitary=sea)
            IB1_store = np.zeros([dimIB[1]-1,dimIB[3]])
            IB2_store = np.zeros([dimIB[1],dimIB[3]])
            IR1_store = np.zeros([dimIB[1]-1,dimIB[3]])
            IR2_store = np.zeros([dimIB[1],dimIB[3]])
            m1_store = np.zeros([dimIB[1]-1,dimIB[3]])
            m2_store = np.zeros([dimIB[1],dimIB[3]])
            for iboot in range(dimIB[3]): #loop boots
                for isea in range(dimIB[1]): #loops sea quark
                    if isea <= 1:
                        IB1_store[isea,iboot] = (IBlst[0])[ic,isea,iv2,iboot]
                        IR1_store[isea,iboot] = (IRlst[0])[ic,isea,iv2,iboot]
                        m1_store[isea,iboot]=(m2_f2lst[0])[0,isea,-1,iv2,iboot]
                    IB2_store[isea,iboot] = (IBlst[1])[ic,isea,iv2,iboot]
                    IR2_store[isea,iboot] = (IRlst[1])[ic,isea,iv2,iboot]
                    #m2_store[isea,iboot]=pow((mlst[1])[0,isea,-1,iv2,iboot]*ainv[1],2)
                    m2_store[isea,iboot]=(m2_f2lst[1])[0,isea,-1,iv2,iboot]
        nameR="./plots/R" + str(ic+1) + "_"+ name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]
        nameB="./plots/B"+str(ic+1)+"_"+name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]
        Rpcent, covcent, Rbp, Rbcov, Rbchisq, Ryp = nlf.bootglobalfit(IR1_store,IR2_store,b_ainv[0],b_ainv[1],m1_store,m2_store,500,p0,pow(mp_phys/(4*pi*fp_phys),2),r"$m^2_{\pi}/4\pi f^2_{\pi}$",r"$R_{"+str(ic+1)+"}$",nameR,IR3_store,b_ainv[2:],m3_store)
        Bpcent, Bcovcent, Bbp, Bbcov, Bbchisq, Byp = nlf.bootglobalfit(IB1_store,IB2_store,b_ainv[0],b_ainv[1],m1_store,m2_store,500,p0,pow(mp_phys/(4*pi*fp_phys),2),r"$m^2_{\pi}/4\pi f^2_{\pi}$",r"$B_{"+str(ic+1)+"}$",nameB,IB3_store,b_ainv[2:],m3_store)
        Rphys.append(Ryp)
        Bphys.append(Byp)

    fo = open("./results/R_B_"+name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]+".dat",'w')
    for i in range(5):
        j=i+1
        fo.write(str(j) + "\t" + str(Rphys[i][-1]) + "\t"  + str(np.std(Rphys[i][0:-1])) + "\t" +  str(Bphys[i][-1]) + "\t" + str(np.std(Bphys[i][0:-1])) + "\n")
    fo.close()

if __name__ == "__main__":
    old_stdout = sys.stdout
    log_file = open("message.log","w")
    sys.stdout = Logger()
    for ischeme in [0]:#range(2):
        for ikin in [0]:#range(5):
            for ibas in [1]:#range(2):
                main(ischeme,ikin,ibas)
    sys.stdout = old_stdout
    log_file.close()
