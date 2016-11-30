from renormalisation.renorm_class import Renormalisation
from fitting.global_fit import GlobalFit
import numpy as np
from datetime import datetime
import sys

#sys.path.insert(0, 'Users/s1035546/My_Global_Fit/dependencies')

class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open(datetime.now().strftime('logfiles/logfile_%H_%M_%d_%m_%Y.log'), "w")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass
old_stdout = sys.stdout
log_file = open("message.log","w")
sys.stdout = Logger()


if __name__ == "__main__":

    for ischeme in range(2):
        for ikin in range(5):
            for ibas in range(2):

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


                    renorm= Renormalisation(ibeta,ischeme,ikin,ibas)
                    B,R,msq,m,f,a2 = renorm.Run_All()
                    if  ibeta < 2:
                        IB,dB = renorm.data_Interp(B)
                        IdB.append(renorm.slope_Interp(dB))
                        IR,dR = renorm.data_Interp(R)
                        IdR.append(renorm.slope_Interp(dR))
                    else:
                        IB = renorm.phys_Interp(B,IdB[ibeta-2])
                        IR = renorm.phys_Interp(R,IdR[ibeta-2])
                    if ibeta == 0:
                        n_ml=len(IR[1])-1 #we discard heaviest term.
                    else:
                        n_ml = len(IR[1])
                    #mat1,mat2 = renorm.CalculateMatrixElem(IB,IR)
                    """
                    for i in range(5):
                        print i
                        if i == 0:
                            print mat1[i,:,:,-1]
                        else:
                            print mat1[i,:,:,-1], mat2[i-1,:,:,-2]
                    """
                    for i in range(n_ml):#lenIR gives the number of mlights
                        if ibeta >1:
                            storeIR.append(IR[:,i,:])
                            storeIB.append(IB[:,i,:])
                        else:
                            storeIR.append(IR[:,i,:])
                            storeIB.append(IB[:,i,:])
                        storem2f2.append(msq[:,i,-1,:])
                        storea2.append(a2)
                        betas.append(ibeta)
                        print "Interpolated B is - ",IB[:,i,-1]
                storeIB=np.array(storeIB)
                storeIR=np.array(storeIR)
                storem2f2=np.array(storem2f2)
                storea2=np.array(storea2)
                fitRlin = GlobalFit(betas,storeIR,storem2f2,storea2,"R",ischeme,ikin,ibas,0)
                fitRlin.SetUp()
                fitRlin.Run()
                fitBlin = GlobalFit(betas,storeIB,storem2f2,storea2,"B",ischeme,ikin,ibas,0)
                fitBlin.SetUp()
                fitBlin.Run()


                fitRchir = GlobalFit(betas,storeIR,storem2f2,storea2,"R",ischeme,ikin,ibas,1)
                fitRchir.SetUp()
                fitRchir.Run()
                fitBchir = GlobalFit(betas,storeIB,storem2f2,storea2,"B",ischeme,ikin,ibas,1)
                fitBchir.SetUp()
                fitBchir.Run()

sys.stdout = old_stdout
log_file.close()
#Now include the fitting class

#Need to ensure this is all done in the correct order
#initialise
#read in all the data into B and R arrays
#read in the Zs
#renormalise 
