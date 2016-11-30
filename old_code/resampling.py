'''
Contains Functions to Bootstrap data
We have N data points, want to resample (with replacement) to get set number Nb 
of bootstraps, each with N data points
'''
import random
import numpy as np

def bootstrap(ydata, xdata, yerr, nb):
    n = len(ydata)

    bydata = np.zeros([nb,n])
    bxdata = np.zeros([nb,n])
    byerr = np.zeros([nb,n])

    for i in range(nb):
        for j in range(n):
            choice = random.randint(0,n-1)
            bydata[i][j] = ydata[choice]
            bxdata[i][j] = xdata[choice]
            byerr[i][j] = yerr[choice]


    return bydata, bxdata, byerr

def jackknife(ydata, xdata, yerr):

    nj = len(ydata)
    jydata = np.zeros([nj,nj-1])
    jxdata = np.zeros([nj,nj-1])
    jyerr = np.zeros([nj,nj-1])

    for i in range(nj):
        for j in range(nj):
            if(j<i):
                jydata[i][j] = ydata[j]
                jxdata[i][j] = xdata[j]
                jyerr[i][j] = yerr[j]
            if(j>i):
                jydata[i][j-1] = ydata[j]
                jxdata[i][j-1] = xdata[j]
                jyerr[i][j-1] = yerr[j]

    return jydata, jxdata, jyerr
                
