
import numpy as np

def read_file_list(filename):
    "reads in data from a file in list format"
    data=[]
    fo = open(filename,'r')
    lines = fo.readlines()
    for line in lines:
        if type(line)==str:
            line = float(line)
            data.append(line)
    return data

def store_list_array(data,array):
    "rearranges list of data into np array"
    counter=0
    for index, val in np.ndenumerate(array):
        array[index] = data[counter]
        counter+=1
    return 1

def read_bootstraps(filename):
    '''
    bootstrap files are files with nboot+1 LINES.
    contain bootstraps + error for one variable
    '''
    #do I need this? does barely anything
    data = read_file_list(filename)
    error = stats.sem(data)
    return data, error

def read_file_array(filename,array):#nc,n_mseal,n_mval1,n_mval2):
    #read the data in a filelist into an array
    data=read_file_list(filename)
    counter = 0
    for index,val in np.ndenumerate(array):
        array[index] = data[counter]
        counter+=1
    return array


def read_results_Jamie(filename):
    from struct import unpack
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

def writeout(data,filename):
    fo=open(filename,'w')
    fo.write(str(np.shape(data)))
    fo.write(str("\n"))
    for index,val in np.ndenumerate(data):
        fo.write(str(val))
        fo.write("\n")
    return 1
