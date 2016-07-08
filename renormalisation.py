#renormalisation.py 
#Julia Kettle
#17/05/16

import numpy as np
import time

#####################################################
def Read_Z(filename):
    #Read the Z factors into 5by5bynbootsmatrix.
    #But bootstraps on 
    nch=5
    nboot=100
    fi = open(filename,'r')
    lines = fi.readlines()
    i=0
    Z_boots=np.zeros([nch,nch,nboot+1])
    Z_cent=np.zeros([nch,nch])
    Z_err=np.zeros([nch,nch])
    for ich in range(nch):
        for jch in range(nch):
            for iboot in range(nboot+1):
                Z_boots[ich][jch][iboot] = lines[i]
                i=i+1
            Z_cent[ich][jch]=Z_boots[ich][jch][-1]
            Z_err=np.std(Z_boots,axis=2)

    return Z_boots, Z_cent, Z_err
#####################################################            

###################################################################################################################################################################################################################        
def ZA_ZP(ikin,ibeta):
	if ibeta > 1:
		ibeta = ibeta-2
        ZA=[0.71273,0.74404]
        ZP=[[0.57732,0.57733],[0.69512,0.69464],[0.65635,0.65832],[0.65635*1.0499,0.65832*1.0499],[0.69513,0.69471],[0.69513*1.0157,0.69471*1.0157],[0.69513,0.69471],[0.65635,0.65832],[0.65635*1.0499,0.65832*1.0499]]
        return ZA[ibeta], ZP[ikin][ibeta]
###################################################################################################################################################################################################################

#########################################################################################################################################################
def MultRZ(R,ZschemeB,T):
    #Scale the ratios R by Z factor
	dimR = np.shape(R)
        ren_R = np.zeros(dimR)
        ren_phys_R = np.zeros(dimR)
        for iml in range(dimR[1]):
            for ims in range(dimR[2]):
                    for imq in range(dimR[3]):
                        for iboot in range(dimR[4]):
                                for ich in range(5):
                                    for jch in range(5):
                                            ren_R[ich,iml,ims,imq,iboot] += ZschemeB[ich,jch]*R[jch,iml,ims,imq,iboot]/ZschemeB[0][0] 
                                for ich in range(5):
                                    for jch in range(5):
                                            ren_phys_R[ich,iml,ims,imq,iboot] += T[ich,jch]*ren_R[jch,iml,ims,imq,iboot] #convert to correct basis, T=I in lattice basis
	return ren_R, ren_phys_R
##########################################################################################################################################################


###################################################################################################################
def MultBZ(B,ZschemeB,ZA,ZP,T,N):
    #Scale the bag B by the Z factors
    dimB = np.shape(B)
    ren_B = np.zeros(dimB)
    ren_phys_B = np.zeros(dimB)
    B_phys = np.zeros(dimB)
    for iml in range(dimB[1]):
        for ims in range(dimB[2]):
                for imq in range(dimB[3]):
                    for iboot in range(dimB[4]):
                            for ich in range(5):
                                for jch in range(5):
                                        ren_B[ich,iml,ims,imq,iboot] += ZschemeB[ich,jch]*B[jch,iml,ims,imq,iboot]*(ZA*ZA)
                            for ich in range(5):
                                for jch in range(5):
                                        ren_phys_B[ich,iml,ims,imq,iboot] += T[ich,jch]*ren_B[jch,iml,ims,imq,iboot]
                            B_phys[0,iml,ims,imq,iboot] = ren_phys_B[0,iml,ims,imq,iboot]/(ZA*ZA*N[0])
                            for ich in range(1,5):
                                            B_phys[ich,iml,ims,imq,iboot] = ren_phys_B[ich,iml,ims,imq,iboot]/(ZP*ZP*N[ich])
    return ren_B, ren_phys_B, B_phys
####################################################################################################################


##################################################################################################
def Renormalise(ibas,ikin,ibeta,Zfilename,R=[],B=[]):
    #Pass the filename with the Zfactors, basis, kinematics and vol to renormalise the ratio and the bag
    N_ren = [8.0/3,4.0/3,-2.0,5.0/3,1.0] #This is the factor in the bag parameter
    #N_ren = [1.0, 1.0, 1.0, 1.0, 1.0]
    ZA, ZP = ZA_ZP(ikin,ibeta)
    #Zfilename = "data/Z" + scheme_name + "_boot_mu_match_" + volname + name_kin +"_block.out"
    print Zfilename
    ZschemeB, ZschemeC, ZschemeE = Read_Z(Zfilename)
    if ibas==1:#susy
        T=np.array([[1,0,0,0,0],[0,0,0,1,0],[0,0,0,-0.5,0.5],[0,0,1,0,0],[0,-0.5,0,0,0]])
    elif ibas==0:#renorm
        T=np.identity(5)
    N=[0,0,0,0,0]
    for i in range(len(N_ren)):
        for j in range(len(N_ren)):
                N[i] +=  N_ren[j]*T[i][j] ##Might need to remove N_ren?? Already used on the bag??
    #5by5byNboot+1 times
    if type(R) == np.ndarray and type(B) == np.ndarray:	
        ren_B, ren_phys_B, B_phys = MultBZ(B,ZschemeC,ZA,ZP,T,N)
        ren_R, ren_phys_R = MultRZ(R,ZschemeC,T)
    elif type(R) == np.ndarray:
        ren_R, ren_phys_R = MultRZ(R,ZschemeC,T)
        B_phys=[]
    elif type(B) == np.ndarray:
        ren_B, ren_phys_B, B_phys = MultBZ(B,ZschemeC,ZA,ZP,T,N)
        ren_phys_R = []

    return ren_phys_R, B_phys
###################################################################################################    
