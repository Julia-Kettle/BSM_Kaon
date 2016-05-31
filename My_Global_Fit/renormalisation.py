#In here we need to read in all the Z factors & apply them.

import numpy as np
import time

def Read_Z(filename):
    nch=5
    nboot=100
    fi = open(filename,'r')
    lines = fi.readlines()
    print "lines", len(lines)
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
            
        
def ZA_ZP(ikin,ibeta):
    ZA=[0.71273,0.74404]
    ZP=[[0.57732,0.57733],[0.69512,0.69464],[0.65635,0.65832],[0.65635*1.0499,0.65832*1.0499],[0.69513,0.69471],[0.69513*1.0157,0.69471*1.0157],[0.69513,0.69471],[0.65635,0.65832],[0.65635*1.0499,0.65832*1.0499]]
    return ZA[ibeta], ZP[ikin][ibeta]

def MultRZ(R,ZschemeB,T):
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
                            ren_phys_R[ich,iml,ims,imq,iboot] += T[ich,jch]*ren_R[jch,iml,ims,imq,iboot]
                            print T[ich,jch], ren_R[jch,iml,ims,imq,iboot]
                            #time.sleep(1)
                        print ich, iml, ims, imq, iboot, ren_R[ich,iml,ims,imq,iboot], ren_phys_R[ich,iml,ims,imq,iboot]
    #time.sleep(3)
    return ren_R, ren_phys_R


def MultBZ(B,ZschemeB,ZA,ZP,T,N):
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
                            ren_B[ich,iml,ims,imq,iboot] += ZschemeB[ich,jch]*B[jch,iml,ims,imq,iboot]*ZA*ZA 
                        for jch in range(5):
                            ren_phys_B[ich,iml,ims,imq,iboot] += T[ich,jch]*ren_B[jch,iml,ims,imq,iboot]*ZA*ZA
                        B_phys[ich,iml,ims,imq,iboot] = ren_phys_B[ich,iml,ims,imq,iboot]/(ZP*ZP*N[ich])
    return ren_B, ren_phys_B, B_phys


def Renormalise(ibas,ikin,ibeta,Zfilename,R,B):
    #Need to know are we in renorm basis or susy basis?
    N_ren = [8/3,4/3,-2,5/3,1] #This is the factor in the bag parameter
    ZA, ZP = ZA_ZP(ikin,ibeta)
    #Zfilename = "data/Z" + scheme_name + "_boot_mu_match_" + volname + name_kin +"_block.out"
    ZschemeB, ZschemeC, ZschemeE = Read_Z(Zfilename)
    if ibas==1:#susy
        T=np.array([[1,0,0,0,0],[0,0,0,1,0],[0,0,0,-0.5,0.5],[0,0,1,0,0],[0,-0.5,0,0,0]])
    elif ibas==0:#renorm
        T=np.identity(5)
    N=[0,0,0,0,0]
    for i in range(len(N_ren)):
        for j in range(len(N_ren)):
            N[i] += N_ren[j]*T[j][i]
    print N
    #5by5byNboot+1 times
    ren_R, ren_phys_R = MultRZ(R,ZschemeC,T)
    ren_B, ren_phys_B, B_phys = MultBZ(B,ZschemeC,ZA,ZP,T,N)

    return ren_phys_R, B_phys
    
    
