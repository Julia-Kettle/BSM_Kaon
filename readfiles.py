'''
File containing all the functions to read files needed
'''

from scipy import stats
from struct import *
import numpy as np

def read_datacol_file(filename):
    fo = open(filename,'r')
    lines = fo.readlines()
    data = []
    for line in lines:
        if type(line)==str:
            line = float(line)
            data.append(line)
    return data

def store_col_3D(x,y,data):
    col = len(data)/(x*y)
    array = np.zeros(x,y,col)
    c=0
    for i in range(x):
        for j in range(y):
            for k in range(col):
                array[i,j,k] = data[c]
                c=c+1
    return data

            
def read_Z(filename,nchan,nboot):
    indx=0
    Z = np.zeros([nchan,nchan,nboot+1])
    data = read_datacol_file(filename)
    for i in range(nchan):
        for j in range(nchan):
            for k in range(nboot+1):
                Z[i,j,k] = data[indx] 
    return Z
                
            

def read_bootstraps(filename):
    '''
    bootstrap files are files with nboot+1 LINES.
    contain bootstraps + error for one variable
    '''
    data = read_datacol_file(filename)
    error = stats.sem(data)
    return data, error

def read_results(filename,nc,n_mseal,n_mval1,n_mval2,nboot):
    print "Reading  " + filename
    print "Reading  " + filename
    data=read_datacol_file(filename)
    cnt = 0
    array = np.zeros([nc,n_mseal,n_mval1,n_mval2,nboot])
    for j in range(n_mseal):
        for k in range(n_mval1):
            for l in range(n_mval2):
                for i in range(nc):
                    for m in range(nboot):
                        array[i,j,k,l,m] = data[cnt]
                        cnt=cnt+1
    return array
    
            

def read_results_Jamie(filename):
    '''
    Need to read in binary files. consist of (4byte int) n, nboot, 
    (8byte double) nboot+1 times dummy data, (4byte int) n, nboot, 
    (8byte double) nboot+1 times real data.
    Use struct module to read in 4 or 8 bytes(depending on file posn)
    then convert back to int or float
    When converting use >, to indicate order, bigendian
    Values gained are bootstraps + central of either fk, mk or mk0fk(crossterms)
    '''
    fo = open(filename,"rb")
    print "Opening " + filename
    print "Opening " + filename
    sect = fo.read(4) #read 4bytes of data (1st n)
    count=0;
    dummy = []
    data = []
    while sect !="":
        count=count+1
        if count<=1: #first n
            n = unpack('>i',sect)[0]#convert from bin to int (bugendian)
            sect = fo.read(4) 
        elif count <=2: #first nboot
            nboot = unpack('>i',sect)[0]
            sect = fo.read(8)
        elif count <= nboot+2: #the dummy variables except last
            dummy.append((unpack('>d',sect))[0])
            sect = fo.read(8)
        elif count<= nboot+3: #last dummy variable - switch back to 4 byte
            dummy.append(unpack('>d',sect)[0])
            sect = fo.read(4)
        elif count <= nboot+4:#2nd n
            n = unpack('>i',sect)[0]
            sect = fo.read(4)
        elif count <= nboot+5:#2nd nboot
            nboot = unpack('>i',sect)[0]
            sect = fo.read(8)
        else: #data
            data.append(unpack('>d',sect)[0])
            sect = fo.read(8)

    return n, nboot, data
            
            
    



    

    
    
