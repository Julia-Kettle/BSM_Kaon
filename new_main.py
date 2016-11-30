from new_class_renorm import Renormalisation
from new_class_fitting import GlobalFit
import numpy as np


for ischeme in range(2):
    for ikin in range(1):
        for ibas in range(1):
            for iBeta in range(4): #loop through each lattice
                print "----------------------setting up as----------------------------"
                #read lattice inverse - again where do we use this?? Currently not using bootstraps for ainve. Should I be?

            ikin = 4
            IdB = []
            IdR = []
            IdM = []
            IdF = []
            storeIR=[]
            storeIB=[]
            storeIM=[]
            storeIF=[]
            storea2=[]
            storem2f2=[]
            betas=[]
            for ibeta in range(4):
                '''
                mpsq_fpsq_phys = pow(0.13957018/(4*pi*0.13),2)
                mp_phys = 0.13957018
                fp_phys = 0.13
                '''


                test= Renormalisation(ibeta,ischeme,ikin,ibas)
                B,R,msq,m,f,a2 = test.Run_All()
                if  ibeta < 2:
                    IB,dB = test.data_Interp(B)
                    IdB.append(test.slope_Interp(dB))
                    #print np.shape(IB), np.shape(IdB[0])
                    IR,dR = test.data_Interp(R)
                    IdR.append(test.slope_Interp(dR))
                    #IM,dM = test.data_Interp(m)
                    #IdM.append(test.slope_Interp(dM))
                    #IF,dF = test.data_Interp(f)
                    #IdF.append(test.slope_Interp(dF))
                else:
                    IB = test.phys_Interp(B,IdB[ibeta-2])
                    IR = test.phys_Interp(R,IdR[ibeta-2])
                    #IM = test.phys_Interp(m,IdM[ibeta-2])
                    #IF = test.phys_Interp(f,IdF[ibeta-2])
                if ibeta == 0:
                    lenIR=len(IR[1])-1
                else:
                    lenIR = len(IR[1])
                print np.shape(IR),np.shape(IB),np.shape(a2),np.shape(m)
                mat1,mat2 = test.CalculateMatrixElem(IB,IR)
                for i in range(5):
                    print i
                    if i == 0:
                        print mat1[i,:,:,-1]
                    else:
                        print mat1[i,:,:,-1], mat2[i-1,:,:,-2]
                for i in range(lenIR):
                    if ibeta >1:
                        storeIR.append(R[:,i,0,:,:])
                        storeIB.append(B[:,i,0,:,:])
                    else:
                        storeIR.append(IR[:,i,:,:])
                        storeIB.append(IB[:,i,:,:])
                    storem2f2.append(msq[:,i,-1,:,:])
                    storea2.append(a2)
                    betas.append(ibeta)
    


            test2 = GlobalFit(betas,storeIR,storem2f2,storea2,"R",ischeme,ikin,ibas,0)
            test2.SetUp()
            test3 = GlobalFit(betas,storeIB,storem2f2,storea2,"B",ischeme,ikin,ibas,0)
            test3.SetUp()

#Now include the fitting class

#Need to ensure this is all done in the correct order
#initialise
#read in all the data into B and R arrays
#read in the Zs
#renormalise 
