'''
File containing all the functions to read files needed
'''



from scipy import stats
from struct import *
import numpy as np
from copy import deepcopy
from decimal import *
from utils import *
from math import pi
from global_settings import *


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
            
            
    
def return_3pts(iBeta,iParams,nboots):
    ###########################################################################################
    #Start reading 3pt functions
    print "Reading 3pts..."
    msea = deepcopy(iParams.m_sea_l)
    ms = deepcopy(iParams.m_val_s)
    latticedim = iParams.latticename
    filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p_over_2p.txt'
    matel_3p_2p = read_results(filename,5,len(msea),len(ms),1,nboots+1) #3pt/(2pt*2pt) ~ Bag
    filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p.txt'
    matel_3p=read_results(filename,5,len(msea),len(ms),1,nboots+1) #3pt
    filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p_rat.txt'
    matel_3p_rat = read_results(filename,4,len(msea),len(ms),1,nboots+1) #3pt1/3pti ~ R
    print "Done" + "\n"
    return matel_3p_2p, matel_3p, matel_3p_rat
    ###########################################################################################



    
###################################################################################################################################################################
def return_myBag(iBeta,iParams,nboots):
    print " Reading Physical point Bag data..."
    ms  = deepcopy(iParams.m_val_s)
    msea = deepcopy(iParams.m_sea_l)
    kaon = 's' + str(ms[0]) + '-l'+ str(msea[0])
    lat = iParams.latticename
    Ninv = [3.0/8,3.0/4,-0.5,3.0/5,1.0] #removes scaling to be consistent with Nicolas, re-added in later.
    B = np.zeros([5,1,1,1,nboots+1])
    R = np.zeros([4,1,1,1,nboots+1])    
    filename=[]
    if iBeta == 2:
        dt = 40
    elif iBeta == 3:
        dt=52
    #loop throught the channels
    for ic in [0,1,2, 3, 4]:
        # 3 & 4 need to switched - convention choice
        if ic == 2:
            filenameR='../Fits/' + lat  + '/3ptrat/channel' +str(2*ic+2) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
            filenameB='../Fits/' + lat  + '/bag/channel' +str(2*ic+2) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
        elif ic ==3:
            filenameR='../Fits/' + lat  + '/3ptrat/channel' +str(2*ic-2) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
            filenameB='../Fits/' + lat  + '/bag/channel' +str(2*ic-2) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
        else:
            filenameR='../Fits/' + lat  + '/3ptrat/channel' +str(2*ic) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
            filenameB='../Fits/' + lat  + '/bag/channel' +str(2*ic) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) +'-dt'+str(dt)+'_boots.dat'
        dataB,error=read_bootstraps(filenameB)
        dataR, error=read_bootstraps(filenameR)
        for iboot in range(nboots+1):
            if ic == 0:
                B[ic,0,0,0,iboot] = dataB[iboot]/Ninv[ic]
            else:
                B[ic,0,0,0,iboot] = dataB[iboot]/Ninv[ic]
                R[ic-1,0,0,0,iboot] = dataR[iboot]        
    print "Done"+"\n"
    return B, R
##################################################################################################################################################################

    
    
#########################################################################################################################################
def return_myMeson(iBetam,iParams,nboots,meson):
    print " Reading Physical point Mass data..."
    #setup masses & lat
    ms = deepcopy(iParams.m_val_s)
    msea = deepcopy(iParams.m_sea_l)
    lat = iParams.latticename
    mval=np.hstack([msea,ms]) # list of light,  strange val quarks
    if meson == 'pion':
        meson = 'l' + str(mval[0]) + '-l' + str(msea[0])
    elif meson == 'kaon':
        meson = 's' + str(mval[1]) + '-l'+ str(msea[0])
    else:
        print 'error - must be a pion or kaon'
    #set up matrices
    b_M = np.zeros([1,1,1,1,nboots+1])
    b_F =  np.zeros([1,1,1,1,nboots+1])
    b_msq_fsq = np.zeros([1,1,1,1,nboots+1])
    #set up filenames & read data
    filenameF='../Fits/' + lat + '/mass/mass-' + meson + '/decay-' + meson +'_boots.dat'
    filenameM='../Fits/' + lat + '/mass/mass-' + meson + '/mass-' + meson +'_boots.dat'
    dataM , error = read_bootstraps(filenameM)
    dataF , error = read_bootstraps(filenameF)
    for iboot in range(nboots): #Temporarily fill with the true phys data, before we have run and fit the lattice data
        b_M[0,0,0,0,iboot] = dataM[iboot]
        b_F[0,0,0,0,iboot] = dataF[iboot]
        b_msq_fsq[0,0,0,0,iboot] = pow(b_M[0,0,0,0,iboot]/(4*pi*b_F[0,0,0,0,iboot]),2)
    b_M[0,0,0,0,-1] = dataM[-1]
    b_F[0,0,0,0,-1] = dataF[-1]
    b_msq_fsq[0,0,0,0,-1] = pow(b_M[0,0,0,0,-1]/(4*pi*b_F[0,0,0,0,-1]),2)
    print "Done" + "\n"
    return b_M, b_F, b_msq_fsq

def return_Jamie_data(iBeta,iParams,nboot):
    #Start reading Jamie's data (2pts)
    print "Reading 2pts..."
    """
    Need to open Jamie's data here in data_Jamie/(24/32)/(amp_par)mu(mu)_mq(mq)_ms(ms)(filend_basis)
    24:  mu = 0.0050,0.0100,0.0200; mq = mu,0.0300,0.0350,0.0400; ms = mq, 0.0300,0.0350,0.0400
    All unitary, includes ss, ll, and ls - only want ll & ls
    3 lights, 3 strange
    32:  mu = 0.0040,0.0060,0.0080; mq = 0.0040, 0.0060, 0.0080; ms = 0.0250,0.0300,mq
    Includes quenched but want only unitary, only includes ls,ll
    3 lights 2 strange
    read_results_Jamie(filename)
    the data returned is a*amp_par
    """

    #get quark masses + physical
    msea = deepcopy(iParams.m_sea_l) #mv1=msea always
    ms = deepcopy(iParams.m_val_s)
    ms_phys = deepcopy(iParams.m_sea_s_phys)
    latticedim = iParams.latticename
    filend_basis = ['_9.22', '_12.52']

    #set up matrices to store results
    b_FK = np.zeros([1,len(msea),len(ms)+1,1,nboot+1])
    b_MK = np.zeros([1,len(msea),len(ms)+1,1,nboot+1])
    b_mpsq_fpsq = np.zeros([1,len(msea),len(ms)+1,1,nboot+1])



    for i in range(len(msea)):
        mv1 = [deepcopy(msea[i])] #1st valence l + sea

        for k in range(len(mv1)):
                mv2 = np.hstack((ms,msea[i])) # 2nd valence, either s or l
                for j in range(len(mv2)):

                    filenameFK = 'data_Jamie/'+latticedim + '/FK_4CHAN_mu' +str5sf(msea[i])+'_mq'+str5sf(mv1[k])+'_ms'+str5sf(mv2[j]) + filend_basis[iBeta]
                    filenameMK = 'data_Jamie/'+latticedim + '/MK_4CHAN_mu' +str5sf(msea[i])+'_mq'+str5sf(mv1[k])+'_ms'+str5sf(mv2[j]) + filend_basis[iBeta]
                    n, nboot, dataFK = read_results_Jamie(filenameFK)
                    n, nboot, dataMK = read_results_Jamie(filenameMK)
                    for l in range(len(dataMK)):
                            b_FK[0,i,j,k,l] = dataFK[l]
                            b_MK[0,i,j,k,l] = dataMK[l]
                            b_mpsq_fpsq[0,i,j,k,l]= pow((b_MK[0,i,j,k,l]/(4*pi*b_FK[0,i,j,k,l])),2)
    print "Done"+"\n"
    return b_FK, b_MK, b_mpsq_fpsq
    ######################################################################################################################################################
