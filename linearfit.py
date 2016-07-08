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
    y0=Interp_run(x0,params)
    return y0, params[-1]
        
def Interp_run(x0,params):
	y = 0
	for i in range(len(params)):
		y = y + params[i]*pow(x0,i)
	return y
        
def R_Interp(R, m_s, m_s_phys, m_l):
    print "-----Interpolating non-physical mass data-----"
    dimR = np.shape(R)
    #shape channels, sea, mval2, boot
    IR = np.zeros([dimR[0],dimR[1],dimR[3],dimR[4]])
    store_dRdms=np.zeros([dimR[0],dimR[1],dimR[3],dimR[4]])
    for ich in range(dimR[0]): #channel
        print '\n'
        for isea in range(dimR[1]): #sea mass
            for iv2 in range(dimR[3]): #valence 2
                yerr=np.std(R[ich,isea],axis=2)
                yerr = yerr[:,0]
                store_params=[]
                for ib in range(dimR[4]): #boot
                    ydata = []
                    for ivs in range(dimR[2]):
						#strange val - #last entry is light
						ydata.append(R[ich,isea,ivs,iv2,ib])
                    IR[ich,isea,iv2,ib], params=Interp(m_s,ydata,yerr,m_s_phys,1)
                    store_dRdms[ich,isea,iv2,ib]=params
                print "channel -",ich, "\t ml -",  m_l[isea]
                for ivs in range(dimR[2]):
                    print "x({}):  {}\ty:   {} +/- {}".format(ivs,m_s[ivs],R[ich,isea,ivs,iv2,-1],np.std(R[ich,isea,ivs,iv2,:500]))
                print "x(p):  {}\typ:  {} +/- {}".format(m_s_phys,IR[ich,isea,iv2,-1],np.std(IR[ich,isea,iv2,:500]))
                print "dR/dms:  {} +/- {}\n".format(store_dRdms[ich,isea,iv2,-1],np.std(store_dRdms[ich,isea,iv2,:500]))
    return IR, store_dRdms
            
def dR_Interp(dR, m_l, m_l_phys):
	print "Interpolate dR/dms to the near physical point data"
	dimdR = np.shape(dR)
	#shape channels, sea, mval2, boot
	IdR = np.zeros([dimdR[0],dimdR[2],dimdR[3]])
	for ich in range(dimdR[0]): #channel
		for iv2 in range(dimdR[2]): #valence 2
			yerr1=np.std(dR[ich,:,:],axis=2)
			for ib in range(dimdR[3]): #boot
				ydata = []
				yerr = []
				for il in range(dimdR[1]):	
					#strange val - #last entry is light
					ydata.append(dR[ich,il,iv2,ib])
					yerr.append(yerr1[il,0])
				IdR[ich,iv2,ib], params=Interp(m_l,ydata,yerr,m_l_phys,1)
				
	return IdR, params

def phys_Interp(y_sim,dRdms,x_p,x_sim):
    print "Interpolating physical pt data to true strange physical point"
    dimy = np.shape(y_sim) #(4/5,1,1,1,501)
    y_p = np.zeros([dimy[0],dimy[1],dimy[2],dimy[4]])
    for ich in range(dimy[0]):
        for isea in range(dimy[1]):
            for iv2 in range(dimy[2]):
                for istr in range(dimy[3]):
                    for iboot in range(dimy[4]):
                        y_p[ich,isea,iv2,iboot]=y_sim[ich,isea,iv2,istr,iboot] + dRdms[ich,iv2,iboot]*(x_p-x_sim)
                print "channel -" , ich
                print "x_sim:  {}\ty_sim:  {}\tdy/dx:  {}\txp:  {}\typ:  {} +/- {}".format(x_sim,y_sim[ich,isea,iv2,istr,-1],dRdms[ich,iv2,-1],x_p,y_p[ich,isea,iv2,-1],np.std(y_p[ich,isea,iv2,:500]))
    return y_p	
