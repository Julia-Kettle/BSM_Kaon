import numpy as np                
                
                
def calcBag(matel_3p_2p):
    dim = np.shape(matel_3p_2p)
    bag = np.zeros(dim)
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                for l in range(dim[3]):
                    for m in range(dim[4]):
                        bag[i,j,k,l,m] = 3*8*matel_3p_2p[i,j,k,l,m]
    return bag

    '''
    Bsm = <K|O1|K>/((8/3)*2(Mk*Fk)^2) = 3pt/(8/3)*2pt*2pt
    Bbsmi = <K|Oi|K> / Ni*(sqrt(2)FkMk^2/(ms+md))^2
    
    R = (Fk/Mk)^2_exper * (Mk/Fk)^2*<P|O1|P>/<P|Oi|P>

    Gi are the golden combinations - no logs so chiral extrapolation easier

    

    '''

def calc_R(Fk,Mk,Rat):
    """
    Calculates the ratio of each channel 3pt over sm 3pt.
    timesed by a factor (Fk/Mk)^2phys*(Mpi/Fpi)^2
    """
    dim = []
    dim1 = np.shape(Rat)
    dim2 = np.shape(Fk)
    for i in range(len(dim1)):
        dim.append(min(dim1[i],dim2[i]))
    dim[0] = 5
    R = np.zeros(dim)
    Mkphys = 0.5*(493.677+497.614)*pow(10,-3)
    Fkphys = 156.2*pow(10,-3)
    #Mkphys = 139.57018
    #Fkphys=130
    #print "R"
    #print "R"
    #print "R"
    #print "R"
    for i in range(dim[1]):
        for j in range(dim[2]):
            for k in range(dim[3]):
                for m in range(dim[4]):
                    for l in range(1,dim[0]):
                        R[0,i,j,k,m] = pow((Fkphys/Mkphys),2)*pow((Mk[0,i,j,k,m]/Fk[0,i,j,k,m]),2)
                        R[l,i,j,k,m] = pow((Fkphys/Mkphys),2)*pow((Mk[0,i,j,k,m]/Fk[0,i,j,k,m]),2)*Rat[l-1,i,j,k,m]
                             
    return R
                        
def renR(R,Z):
    dimZ = np.shape(Z)
    dimR = np.shape(R)


def calc_G(B):
        #the B are the bag paramaters 
        #we pass the 5 channels so B[5][nboot+1] - loop throught the masses and ibeta outside
    G = np.zeros(np.shape(B))
    for i in range(len(B[0])):
        G[0][i] = B[2][i]/B[3][i] #23
        G[1][i] = B[4][i]/B[5][i] #45
        G[2][i] = B[2][i]*B[4][i] #24
        G[3][i] = B[12][i]/B[1][i] #21
        G[4][i] = B[3][i]/B[1][i] #31
        
    return G














