'''
Perform the Global Fit. 
Read Jamie's 2pt data for 24, 32
Read Nicolas' 3pt data for 24, 32
Calculate ratios
Renormalise
Interpolate to physical strange
Perform global fit against pion mass/decay sqr ratio and lattice spacing


Need to include the new physical point data too. 
'''


import nonlinearfit as nlf
import linearfit as lf
import numpy as np
import random
import readfiles as rf
import resampling as rsamp
import setup
from utils import *
from copy import deepcopy
from decimal import *
import scipy as sp
from calculations import *
from math import *
from renormalisation import *



def return_Jamie_data(iBeta,iParams,nboot):
	#Start reading Jamie's data (2pts)
    	print "Reading 2pts..."
    	"""
    	Need to open Jamie's data here in data_Jamie/(24/32)/(amp_par)mu(mu)_mq(mq)_ms(ms)(filend_basis)
    	24:  mu = 0.0050,0.0100,0.0200; mq = mu,0.0300,0.0350,0.0400; ms = 0.0300,0.0350,0.0400
    	32:  mu = 0.0040,0.0060,0.0080; mq = 0.0040, 0.0060, 0.0080; ms = 0.0250,0.0300,mq
    	read_results_Jamie(filename)
    	the data returned is a*amp_par
    	"""

    	#get quark masses + physical
    	msea = deepcopy(iParams.m_sea_l)
    	ms = deepcopy(iParams.m_val_s)
    	ms_phys = deepcopy(iParams.m_sea_s_phys)
    
    	#set up matrices to store results
    	b_FK = np.zeros([1,len(msea),len(ms)+1,1,nboot+1])
    	b_MK = np.zeros([1,len(msea),len(ms)+1,1,nboot+1])
    	b_mpsq_fpsq = np.zeros([1,len(msea),len(ms)+1,1,nboot+1])

	
    
    	for i in range(len(msea)): 
        	mv1 = [deepcopy(msea[i])] #1st valence l + sea

        	for k in range(len(mv1)):
            		mv2 = np.hstack((ms,msea[i])) # 2nd valence, either s or l
            		for j in range(len(mv2)):

                		filenameFK = 'data_Jamie/'+latticedim[iBeta]+amp[1]+'mu'+str5sf(msea[i])+'_mq'+str5sf(mv1[k])+'_ms'+str5sf(mv2[j]) + filend_basis[iBeta]
                		filenameMK = 'data_Jamie/'+latticedim[iBeta]+amp[0]+'mu'+str5sf(msea[i])+'_mq'+str5sf(mv1[k])+'_ms'+str5sf(mv2[j]) + filend_basis[iBeta]
                		n, nboot, dataFK = rf.read_results_Jamie(filenameFK)
                		n, nboot, dataMK = rf.read_results_Jamie(filenameMK)
                		for l in range(len(dataMK)):
                    			b_FK[0,i,j,k,l] = dataFK[l]
                    			b_MK[0,i,j,k,l] = dataMK[l]
                    			b_mpsq_fpsq[0,i,j,k,l]= pow((b_MK[0,i,j,k,l]/(4*pi*b_FK[0,i,j,k,l])),2)
	return b_FK, b_MK, b_mpsq_fpsq
    ######################################################################################################################################################


def return_3pts(iBeta,iParams,nboots):
    	###########################################################################################
    	#Start reading 3pt functions
    	print "Reading 3pts..."	
    	msea = deepcopy(iParams.m_sea_l)
    	ms = deepcopy(iParams.m_val_s)
    	filename = 'data/boot_matel_'+latticedim[iBeta]+'cubed_IW_3p_over_2p.txt'
    	matel_3p_2p = rf.read_results(filename,5,len(msea),len(ms),1,nboot+1) #3pt/(2pt*2pt) ~ Bag
    	filename = 'data/boot_matel_'+latticedim[iBeta]+'cubed_IW_3p.txt'
    	matel_3p=rf.read_results(filename,5,len(msea),len(ms),1,nboot+1) #3pt
    	filename = 'data/boot_matel_'+latticedim[iBeta]+'cubed_IW_3p_rat.txt'
   	matel_3p_rat = rf.read_results(filename,4,len(msea),len(ms),1,nboot+1) #3pt1/3pti ~ R
    	return matel_3p_2p, matel_3p, matel_3p_rat
	###########################################################################################

def return_myBag(iBeta,iParams,nboots):
	print " Reading Physical point Bag data..."
	ms  = deepcopy(iParams.m_val_s)
	msea = deepcopy(iParams.m_sea_l)
	for ic in [1,2, 3, 4, 5]:
		filename = 'Julia_data/bag/s' + str(ms) + '-l'+ str(msea) + '/dt40/s'+ str(ms) + '-l' + str(msea) + '-channel'+ic+'.dat'
		print filename
		rf.read_bootstraps(filename)



latticedim = ['24','32']
volname = ['24cubed','32cubed']
beta = [0,1]
amp = ['/MK_4CHAN_','/FK_4CHAN_']
conv = 2 #l-s & s-s
name_kin = ['E','gg','qq','qg''gq']
name_scheme = ['MOM','ms']
name_basis = ['Lattice', 'SUSY']
filend_basis = ['_9.22', '_12.52']


params = [] #holds class with all data for each lattice.
b_ainvs=[]
b_ZV=[]
nboot = 500

#Initilasie lists which hold all lattices data
IRlst=[]
IBlst=[]
Rlst = []
Blst = []
mlst=[]
flst = []
m2_f2lst=[]
ainv=[]

###################
# set up kinematics 
#eventually will need to loop through kinematics - put whole file within function, passing these values??
kin=0
ikin=0
scheme=0
ibas=0
###################

for iBeta in range(2): #loop through each lattice
    
    #read lattice inverse - again where do we use this?? Currently not using bootstraps for ainve. Should I be?
    filename = 'common_data/boot_ainv_'+latticedim[iBeta]+'cubed_IW_'+str(nboot)
    bootainv , err = rf.read_bootstraps(filename)
    b_ainvs.append(bootainv)

    #read ZV of 2pt - why is this here? Do we not do this again within renormalisation??
    filename = 'common_data/boot_ZV_'+latticedim[iBeta]+'cubed_IW_'+str(nboot)
    bootzv, err = rf.read_bootstraps(filename)
    b_ZV.append(bootzv)

mpsq_fpsq_phys = pow(0.13957018/(4*pi*0.13),2)


for iBeta in range(2): #loop through lattice
	
    	#class cotaining paramaters set up    
    	iParams = setup.Params(iBeta)
    	params.append(iParams)   
	ainv.append(deepcopy(iParams.ainv))	

	if iBeta < 2:
		b_MK, b_FK, b_mpsq_fpsq = return_Jamie_data(iBeta,iParams,nboot)
		matel_3p_2p, matel_3p, matel_3p_rat = return_3pts(iBeta,iParams,nboot)
	else:

		return_myBag(iBeta,iParams,nboot)
		#read all the new physical files here
		#need to store the bootstraps first - txt files should do
    	####################################################################################################################################################



    	##########################################
    	#Calculate the Bag paramater and Ratio R
    	#R = 
    	print "Calculate bag and R ..."
    	B = calcBag(matel_3p_2p)
    	R = calc_R(b_FK,b_MK,matel_3p_rat)
    	##########################################
    	#print R[:,:,0:2,:,-1]


    	##########################################################################################################################
    	#Renormlise B and R
    	print name_scheme[scheme],latticedim[iBeta],name_kin[kin]
    	Zfilename = "data/Z" + name_scheme[scheme] + "_boot_mu_match_" + volname[iBeta] + "_" + name_kin[kin] + "_block.out" 
    	Rren_phys, Bren_phys = Renormalise(ibas,kin,iBeta,Zfilename,R,B)
    	##########################################################################################################################

    
    	###########################################################
    	#Interpolate results to physical masses
    	
    	msea = deepcopy(iParams.m_sea_l)
    	ms = deepcopy(iParams.m_val_s)
    	ms_phys = deepcopy(iParams.m_sea_s_phys)

	IR = lf.R_Interp(Rren_phys[:,:,0:2,:,:], ms[0:2], ms_phys)
    	IB = lf.R_Interp(Bren_phys[:,:,0:2,:,:], ms[0:2], ms_phys)
    	###########################################################
  

    	#############################
    	#Append matrices of results
    	#to lists
    	Blst.append(B)
    	Rlst.append(R)
    	IRlst.append(IR)
    	IBlst.append(IB)
    	mlst.append(b_MK)
    	flst.append(b_FK)
    	m2_f2lst.append(b_mpsq_fpsq)
    	#############################

    

p0 = [1,1,1]



#p, cov, bp, bcov, bchisq = nlf.bootglobalfit(y1[0,0,0],y2[0,0,0],1/ainv[0],1/ainv[1],m1[0,0,0],m2[0,0,0],500,p0,mpsq_fpsq_phys)
#actually do the global fit (looping through the sea masses)


dimIR = np.shape(IR)
dimIB = np.shape(IB)
Rphys=[]
Bphys=[]


print "ainv", ainv

for ic in range(dimIR[0]): #loop channels
	for iv2 in range(dimIR[2]): #loop valence 2 (=1 here, unitary=sea)
        	IB1_store = np.zeros([dimIB[1]-1,dimIB[3]]) 
     	   	IB2_store = np.zeros([dimIB[1],dimIB[3]])
     		IR1_store = np.zeros([dimIR[1]-1,dimIR[3]])
        	IR2_store = np.zeros([dimIR[1],dimIR[3]])
        	m1_store = np.zeros([dimIR[1]-1,dimIR[3]])  
        	m2_store = np.zeros([dimIR[1],dimIR[3]])
        	for iboot in range(dimIR[3]): #loop boots
      			for isea in range(dimIR[1]): #loops sea quark
                		if isea <= 1:
                    			IB1_store[isea,iboot] = (IBlst[0])[ic,isea,iv2,iboot]
                    			IR1_store[isea,iboot] = (IRlst[0])[ic,isea,iv2,iboot]
                    			m1_store[isea,iboot]=(m2_f2lst[0])[0,isea,-1,iv2,iboot]
                		IB2_store[isea,iboot] = (IBlst[1])[ic,isea,iv2,iboot]
                		IR2_store[isea,iboot] = (IRlst[1])[ic,isea,iv2,iboot]
                		m2_store[isea,iboot]=(m2_f2lst[1])[0,isea,-1,iv2,iboot]
        
        	nameR="R" + str(ic) + "_"+ name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]
        	nameB="B"+str(ic)+"_"+name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]
       		Rpcent, covcent, Rbp, Rbcov, Rbchisq, Ryp = nlf.bootglobalfit(IR1_store,IR2_store,1/(ainv[0]*ainv[0]),1/(ainv[1]*ainv[1]),m1_store,m2_store,500,p0,mpsq_fpsq_phys,r"$m^2_{\pi}/(4\pi f_{\pi})^2$",r"$R_{"+str(ic+1)+"}$",nameR)
       		Bpcent, Bcovcent, Bbp, Bbcov, Bbchisq, Byp = nlf.bootglobalfit(IB1_store,IB2_store,1/(ainv[0]*ainv[0]),1/(ainv[1]*ainv[1]),m1_store,m2_store,500,p0,mpsq_fpsq_phys,r"$m^2_{\pi}/(4\pi f_{\pi})^2$",r"$B_{"+str(ic+1)+"}$",nameB)
       		Rphys.append(Ryp)
       		Bphys.append(Byp)


fo = open("R_B_"+name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]+".dat",'w')
for i in range(5):
    	j=i+1
    	fo.write(str(j) + "\t" + str(Rphys[i]) + "\t" + str(Bphys[i]) + "\n")

fo.close()
