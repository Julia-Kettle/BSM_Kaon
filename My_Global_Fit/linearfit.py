import numpy as np
from utils import *
import time


            

def lineargenfit(xdata,ydata,yerr,pdeg,fun):

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

def Interp(xdata, ydata, yerr, x0, deg):
    params, cov, err, chisq = linearpolyfit(xdata,ydata, yerr, deg)
    y0 = 0
    for i in range(deg+1):
        y0 = y0 + params[i]*pow(x0,i)

    return y0
        
        
            
def R_Interp(R, m_s, m_s_phys):
    print "Set up R for Interp"
    dimR = np.shape(R)
    print dimR
    #shape channels, sea, mval2, boot
    IR = np.zeros([dimR[0],dimR[1],dimR[3],dimR[4]])
    for ich in range(dimR[0]): #channel
        for isea in range(dimR[1]): #sea mass
            for iv2 in range(dimR[3]): #valence 2
                yerr=np.std(R[ich,isea],axis=2)
                yerr = yerr[:,0]
                for ib in range(dimR[4]): #boot
                    ydata = []
                    for ivs in range(dimR[2]):
                        print ich, isea, iv2, ivs, ib, R[ich,isea,ivs,iv2,ib]
                             #strange val - #last entry is light
                        ydata.append(R[ich,isea,ivs,iv2,ib])
                    IR[ich,isea,iv2,ib]=Interp(m_s,ydata,yerr,m_s_phys,1)
                    
    return IR
                    
