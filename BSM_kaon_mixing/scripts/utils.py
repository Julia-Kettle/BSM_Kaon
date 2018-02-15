import numpy as np
import os
'''
def str5sf(number):
    num_str = str(number)
    while len(num_str) <= 5:
        num_str = num_str+'0'
    return num_str

def my_reshape(M):
    # change order to nc,mq,ms,mu,nboot
    dim = np.shape(M)
    M2 = np.zeros([dim[0],dim[3],dim[2],dim[1],dim[4]])
    for i in range(dim[0]):
        for j in range(dim[1]):
            for k in range(dim[2]):
                for l in range(dim[3]):
                    for m in range(dim[4]):
                        M2[i,l,k,j,m] = M[i,j,k,l,m]
    return M2

def matmul(A,B):
    n = A.shape[0]
    m = B.shape[1]
    l1 = A.shape[1]
    l2 = B.shape[0]

    if l1==l2:
        C = np.zeros([n,m])
        for i in range(n):
            for j in range(m):
                sum = 0
                for k in range(l1):
                    sum = sum + A[i][k]*B[k][j]
                C[i,j] = sum

        return C
    else:
        print "Error, matrices dimensions are wrong"
        print "Error, matrices dimensions are wrong"
        return -1

def transpose(A):
    l1 = A.shape[0]
    l2 = A.shape[1]
    
    B=np.zeros([l2,l1])
    for i in range(l1):
        for j in range(l2):
            B[j,i] = A[i,j]
    
    return B

def diagonalize(A):
    dim = A.shape[0]
    if(dim != A.shape[1]):
        print "Matrix not square"
        print "Matrix not square"
        return -1
    else:
        B = np.zeros([dim,1])
        for i in range(dim):
            B[i] = A[i,i]
        return B
'''
def chiral_log_coeffs(param,ansatz,basis):
    #ugly code to return chiral log coeff
    if ansatz == 'linear':
        coeff = [0,0,0,0,0]
    elif ansatz == 'chiral':
        if param=='R':
            if basis=='SUSY':
                coeff=[1.5,1.5,1.5,2.5,2.5]
            elif basis=='Lattice':
                coeff=[1.5,2.5,2.5,1.5,1.5]
        elif param=='B':
            if basis=='SUSY':
                coeff=[-0.5,-0.5,-0.5,0.5,0.5]
            elif basis=='Lattice':
                coeff=[-0.5,0.5,0.5,-0.5,-0.5]
    return coeff

def arg_check(param,basis,scheme,projscheme,ansatz):
    boolExit=False #check variables and exit if at least one wrong
    if param not in ['R','B']:
        print("Error: param is " + str(param) + ": must be 'R' or 'B'")
        boolExit=True
    if basis not in ['SUSY','Lattice']:
        print("Error: basis is " + str(basis) + ": must be 'SUSY' or 'Lattice'")
        boolExit=True
    if scheme not in ['MOM','ms']:
        print("Error: schem is " + str(scheme) + ": must be 'MOM' or 'ms'")
        boolExit=True
    if projscheme not in ['qq','gg']:
        print("Error: projscheme is " + str(projscheme) + ": must be 'qq' or 'gg'")
        boolExit=True
    if boolExit==True:
        print("Exiting program")
        sys.exit(-1)

def createDir(savelocation):
    #check the save location exists and if not creat it
    if not os.path.exists(os.path.dirname(savelocation)):
        try:
            os.makedirs(os.path.dirname(savelocation))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise



