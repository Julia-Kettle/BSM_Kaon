import numpy as np

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
        return -1
    else:
        B = np.zeros([dim,1])
        for i in range(dim):
            B[i] = A[i,i]
        return B
            
'''
def lineargenfit(xdata,ydata,yerr,pdeg,fun):
    Want to modify to be able to include other types of functions of x, still 
    linear with regards to the parameters though

    """
    funntion to calculate the best fit of a polynomial of degree pdeg.
    xdata, ydata and the error in y, yerr is inputted.
    chisq minimisation is done.
    chisq = sum_i (ydatai - sum_j (pj*xi^j)/yerri)^2
    differentiate to minimise gives

    0 = sum_i( 1/yerri^2 * (ydatai - sum_j(pj*xi^j))*xi^j-1

    In matrix form can be expressed
    AT*y = AT*A*p
    where A_ij = xi^j/yerr_i, pj is parameter, yi is ydatai/yerri
    p = inv(AT*A) *AT*y = C*AT*y
    where C is the covariance matrix
    the diagonals of C give the square of standard errors on params
    """
    
    A = np.zeros([len(ydata),pdeg+1])
    yT=np.zeros([len(ydata)])
    for j in range(pdeg+1):
        sum = 0
        for i in range(len(xdata)):
            #calculate A & y(transpose)
            A[i][j] = pow(xdata[i],(j))/yerr[i]
            yT[i] = ydata[i]/yerr[i]

    
    AT = transpose(A)
    ATA = matmul(AT,A)
    ATAinv = np.linalg.inv(ATA) 
    yT = np.mat(yT) 
    y = transpose(yT)
    B = matmul(ATAinv,AT) # matrix B = inv(A*AT) * AT
    BT = transpose(B)
    C = ATAinv #covariance
    e = diagonalize(C) #get vector of errors^2
    for i in range(len(e)):
        e[i] = pow(e[i],0.5) 
    p = matmul(B,y)

    chisq = 0 #initialise
    for i in range(len(ydata)):
        temp = y[i]
        for j in range(p.shape[0]):
            temp = temp - p[j]*A[i,j] #sum over params*xdata
        chisq = chisq + pow(temp,2)  

    return p, C, e, chisq


def linearpolyfit(xdata,ydata,yerr,pdeg):

    """
    funntion to calculate the best fit of a polynomial of degree pdeg.
    xdata, ydata and the error in y, yerr is inputted.
    chisq minimisation is done.
    chisq = sum_i (ydatai - sum_j (pj*xi^j)/yerri)^2
    differentiate to minimise gives

    0 = sum_i( 1/yerri^2 * (ydatai - sum_j(pj*xi^j))*xi^j-1

    In matrix form can be expressed
    AT*y = AT*A*p
    where A_ij = xi^j/yerr_i, pj is parameter, yi is ydatai/yerri
    p = inv(AT*A) *AT*y = C*AT*y
    where C is the covariance matrix
    the diagonals of C give the square of standard errors on params
    """
    
    A = np.zeros([len(ydata),pdeg+1])
    yT=np.zeros([len(ydata)])
    for j in range(pdeg+1):
        sum = 0
        for i in range(len(xdata)):
            #calculate A & y(transpose)
            A[i][j] = pow(xdata[i],(j))/yerr[i]
            yT[i] = ydata[i]/yerr[i]

    
    AT = transpose(A)
    ATA = matmul(AT,A)
    ATAinv = np.linalg.inv(ATA) 
    yT = np.mat(yT) 
    y = transpose(yT)
    B = matmul(ATAinv,AT) # matrix B = inv(A*AT) * AT
    BT = transpose(B)
    C = ATAinv #covariance
    e = diagonalize(C) #get vector of errors^2
    for i in range(len(e)):
        e[i] = pow(e[i],0.5) 
    p = matmul(B,y)

    chisq = 0 #initialise
    for i in range(len(ydata)):
        temp = y[i]
        for j in range(p.shape[0]):
            temp = temp - p[j]*A[i,j] #sum over params*xdata
        chisq = chisq + pow(temp,2)  

    return p, C, e, chisq
'''  
'''    
import random

x = range(1000)
y = []
for i in x:
    y.append(2*pow(i,2) - 4*i +11)
    y[i] = y[i] + random.randrange(100)/500.0
print y
err = random.sample(xrange(1,1001), 1000)
for i in range(len(err)):
    err[i] = err[i]/1000.0


mat = np.mat(x)



parameters, covariance, errors, chisq = linearpolyfit(x,y,err,2)

print parameters
print covariance
print errors
print chisq/1000.0
'''
        
            
        
        
