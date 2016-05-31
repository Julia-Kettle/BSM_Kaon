'''
File containing all the functions to read files needed
'''
from struct import *

def readbootstraps(filename):
    '''
    bootstrap files are files with nboot+1 LINES.
    contain bootstraps + error for one variable
    '''

    fo = open(filename,'r')
    lines = fo.readlines()
    bootstrap = lines[0:-1]
    error = lines[-1]

    return bootstrap, error

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
            print nboot+1
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
            
            
    



    

    
    
