from global_fit import GlobalFit
import numpy as np
from lattice import Lattice

ischeme=1
ikin=2
ibas=0

#datapoints = ['24 0.005 local old','24 0.01 local old','32 0.004 local old','32 0.006 local old','32 0.008 local old']#'48 0.00078 smeared new','64 0.000678 local new']
datapoints = ['24 0.005 local new','24 0.01 local new','32 0.004 local new','32 0.006 local new','48 0.00078 smeared new','64 0.000678 local new','48 0.002144 local new fine']
betas = [0,0,1,1,2,3,4]
#betas = [2,3]
B=[]
R=[]
m2=[]
a2=[]

for entry in datapoints:
    tokens = entry.split()
    smearbool = True if tokens[2]=='smeared' else False
    newold = "_old" if tokens[3]=='old' else ""
    finebool = True if len(tokens) > 4 and tokens[4]=='fine' else False
    lattice = Lattice(int(tokens[0]),smearbool,fine=finebool)
    print "lattice name - ",lattice.name
    filenameend = lattice.name+"ml_"+tokens[1]+"_ms_"+str(lattice.m_s_phys[0])+newold+".bin"

    Bin = "../renormalised/B_"+filenameend
    Rin = "../renormalised/R_"+filenameend
    m2in = "../renormalised/m_4pif_sq_"+filenameend
    #a2in = "../unrenormalised/m_4pif_sq"+lattice.name+".bin"
    fo = open(Rin,'r')
    fo.close()


    B.append(np.load(Bin))
    R.append(np.load(Rin))
    m2.append(np.load(m2in))
    print "asquared is - ",lattice.a2[-1]
    a2.append(lattice.a2)

storeR=np.squeeze(np.array(R))
storeB=np.squeeze(np.array(B))
storem2=np.squeeze(np.array(m2))
storea2=np.squeeze(np.array(a2))
            
print storem2[:,-1]
print storeR[:,0,-1]
fitR = GlobalFit(betas,storeR,storem2,storea2,"R",ischeme,ikin,ibas,2)
fitR.SetUp()
fitR.Run()
fitB = GlobalFit(betas,storeB,storem2,storea2,"B",ischeme,ikin,ibas,2)
fitB.SetUp()
fitB.Run()
