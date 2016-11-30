import numpy as np
from utils import *
import time


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
    AT = np.transpose(A)
    ATA = np.dot(AT,A)
    ATAinv = np.linalg.inv(ATA) 
    y = np.transpose(yT)
    B = np.dot(ATAinv,AT) # matrix B = inv(A*AT) * AT
    C = ATAinv #covariance
    e = np.diag(C).copy() #get vector of errors^2
    for i in range(len(e)):
        e[i] = pow(e[i],0.5) 

    p = np.dot(B,y)

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
        y0 += params[i]*pow(x0,i)
    return y0, params[-1]

def Interp_run(x0,params):
    y = 0
    for i in range(len(params)):
        y = y + params[i]*pow(x0,i)
    return y

            

