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

from random import *
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
	print "Done"+"\n"
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
    	print "Done" + "\n"
	return matel_3p_2p, matel_3p, matel_3p_rat
	###########################################################################################

def return_myBag(iBeta,iParams,nboots):
	print " Reading Physical point Bag data..."
	ms  = deepcopy(iParams.m_val_s)
	msea = deepcopy(iParams.m_sea_l)
	meson = 's' + str(ms[0]) + '-l'+ str(msea[0])
	lat = iParams.latticename
	#Ninv = [3.0/8,3.0/4,-0.5,3.0/5,1.0] #Need to remove the scaling for consistency with Nicolas' data - is added in later. 
	Ninv = [3.0/8,1.0,1.0,1.0,1.0]
	B = np.zeros([5,1,1,1,nboots+1])
	b_MK = np.zeros([1,1,1,1,nboots+1])
	b_FK =  np.zeros([1,1,1,1,nboots+1])
	b_mpsq_fpsq = np.zeros([1,1,1,1,nboots+1])
	for j in range(nboots): #Temporarily fill with the true phys data, before we have run and fit the lattice data
		b_MK[0,0,0,0,j] = 0.13957018 + (randrange(0,99))/1000.0
		b_FK[0,0,0,0,j] = 0.13 + (randrange(0,100))/1000.0
		b_mpsq_fpsq[0,0,0,0,j] = pow(b_MK[0,0,0,0,j]/(4*pi*b_FK[0,0,0,0,j]),2)
	b_MK[0,0,0,0,-1] = 0.13957018
	b_FK[0,0,0,0,-1] = 0.13 
	b_mpsq_fpsq[0,0,0,0,-1] = pow(b_MK[0,0,0,0,-1]/(4*pi*b_FK[0,0,0,0,-1]),2)
	for ic in [1,2, 3, 4, 5]:
		filename = 'Julia_data/bag/' + lat  + '/bag-' + meson + '/dt40/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) + '-channel'+str(ic)+'.dat'
		data,error=rf.read_bootstraps(filename)
		for iboot in range(nboots+1):
			B[ic-1,0,0,0,iboot] = data[iboot]/Ninv[ic-1]
	print "Done"+"\n"
	return B, b_MK, b_FK, b_mpsq_fpsq



latticedim = ['24','32','48','48','64']
volname = ['24cubed','32cubed','48cubed','48cubedfine','64cubed']
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


for iBeta in range(2): #loop through 24, and 32 (Non-Physical)
	
    	#class cotaining paramaters set up    
    	iParams = setup.Params(iBeta)
    	params.append(iParams)   
	ainv.append(deepcopy(iParams.ainv))	

	b_MK, b_FK, b_mpsq_fpsq = return_Jamie_data(iBeta,iParams,nboot) #read and return meson mass + decay
	matel_3p_2p, matel_3p, matel_3p_rat = return_3pts(iBeta,iParams,nboot) #read and return 3pt functions

	##########################################
    	#Calculate the Bag paramater and Ratio R 
    	print "Calculate bag and R ..."
    	B = matel_3p_2p
    	R = calc_R(b_FK,b_MK,matel_3p_rat)
    	##########################################


    	##########################################################################################################################
    	#Renormlise B and R
    	Zfilename = "data/Z" + name_scheme[scheme] + "_boot_mu_match_" + volname[iBeta] + "_" + name_kin[kin] + "_block.out" 	
	Rren_phys,Bren_phys = Renormalise(ibas,kin,iBeta,Zfilename,R=R,B=B)
    	##########################################################################################################################

    
    	###########################################################
    	#Interpolate results to physical masses
    	
    	msea = deepcopy(iParams.m_sea_l)
    	ms = deepcopy(iParams.m_val_s)
    	ms_phys = deepcopy(iParams.m_sea_s_phys)

	#print Bren_phys[:,:,0:2,:,:]

	IR = lf.R_Interp(Rren_phys[:,:,:,:,:], ms[:], ms_phys)
    	IB = lf.R_Interp(B[:,:,:,:,:], ms[:], ms_phys)
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

print "Finished reading + interpolating Nicolas' data"

for iBeta in [2,3]: #loop through the first 2 phsy data lattices (reuse Z frm 24 & 32)  	
    	iParams = setup.Params(iBeta)
    	params.append(iParams)   
	ainv.append(deepcopy(iParams.ainv))	
	B_new, b_MK, b_FK, b_mpsq_fpsq = return_myBag(iBeta,iParams,nboot)
	#read all the new physical files here		
	#still need the R data + pion mass + decays
    	
	##########################################################################################################################
    	#Renormlise B and R
    	Zfilename = "data/Z" + name_scheme[scheme] + "_boot_mu_match_" + volname[iBeta-2] + "_" + name_kin[kin] + "_block.out" 
	Rdummy, Bren_phys = Renormalise(ibas,kin,iBeta,Zfilename,B=B_new)
    	##########################################################################################################################

	#no need to interpolate here, already at the physical strange

	Blst.append(Bren_phys)
	#Rlst.append(R)
	#IRlst.append(R)
	IBlst.append(B_new)
	#IBlst.append(Bren_phys)
	mlst.append(b_MK)
	flst.append(b_FK)
	m2_f2lst.append(b_mpsq_fpsq)
	
p0 = [1,1,1]



#p, cov, bp, bcov, bchisq = nlf.bootglobalfit(y1[0,0,0],y2[0,0,0],1/ainv[0],1/ainv[1],m1[0,0,0],m2[0,0,0],500,p0,mpsq_fpsq_phys)
#actually do the global fit (looping through the sea masses)


#dimIR = np.shape(IR)
dimIB = np.shape(IB)
#Rphys=[]
Bphys=[]

a_store = np.zeros([2])
a_store[0] = 1/(ainv[2]*ainv[2])
a_store[1] = 1/(ainv[3]*ainv[3])


for ic in range(dimIB[0]): #loop channels
	IB3_store = np.zeros([2,nboot+1])
	m3_store = np.zeros([2,nboot+1])

	for iboot in range(nboot+1):
		IB3_store[0,iboot] = (IBlst[2])[ic,0,0,0,iboot]
		IB3_store[1,iboot] = (IBlst[3])[ic,0,0,0,iboot]
		m3_store[0,iboot] = (m2_f2lst[2])[0,0,0,0,iboot]
		m3_store[1,iboot] = (m2_f2lst[3])[0,0,0,0,iboot]

	for iv2 in range(dimIB[2]): #loop valence 2 (=1 here, unitary=sea)
        	IB1_store = np.zeros([dimIB[1]-1,dimIB[3]]) 
     	   	IB2_store = np.zeros([dimIB[1],dimIB[3]])
     		#IR1_store = np.zeros([dimIR[1]-1,dimIR[3]])
        	#IR2_store = np.zeros([dimIR[1],dimIR[3]])
        	m1_store = np.zeros([dimIB[1]-1,dimIB[3]])  
        	m2_store = np.zeros([dimIB[1],dimIB[3]])
        	for iboot in range(dimIB[3]): #loop boots
      			for isea in range(dimIB[1]): #loops sea quark
                		if isea <= 1:
                    			IB1_store[isea,iboot] = (IBlst[0])[ic,isea,iv2,iboot]
                    			#IR1_store[isea,iboot] = (IRlst[0])[ic,isea,iv2,iboot]
                    			m1_store[isea,iboot]=(m2_f2lst[0])[0,isea,-1,iv2,iboot]
                		IB2_store[isea,iboot] = (IBlst[1])[ic,isea,iv2,iboot]
                		#IR2_store[isea,iboot] = (IRlst[1])[ic,isea,iv2,iboot]
                		m2_store[isea,iboot]=(m2_f2lst[1])[0,isea,-1,iv2,iboot]
        
        	nameR="R" + str(ic) + "_"+ name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]
        	nameB="B"+str(ic)+"_"+name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]
       		#Rpcent, covcent, Rbp, Rbcov, Rbchisq, Ryp = nlf.bootglobalfit(IR1_store,IR2_store,1/(ainv[0]*ainv[0]),1/(ainv[1]*ainv[1]),m1_store,m2_store,500,p0,mpsq_fpsq_phys,r"$m^2_{\pi}/(4\pi f_{\pi})^2$",r"$R_{"+str(ic+1)+"}$",nameR)
       		Bpcent, Bcovcent, Bbp, Bbcov, Bbchisq, Byp = nlf.bootglobalfit(IB1_store,IB2_store,1/(ainv[0]*ainv[0]),1/(ainv[1]*ainv[1]),m1_store,m2_store,500,p0,mpsq_fpsq_phys,r"$m^2_{\pi}/(4\pi f_{\pi})^2$",r"$B_{"+str(ic+1)+"}$",nameB,IB3_store,a_store,m3_store)
       		#Rphys.append(Ryp)
       		Bphys.append(Byp)


fo = open("R_B_"+name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]+".dat",'w')
for i in range(5):
    	j=i+1
    	#fo.write(str(j) + "\t" + str(Rphys[i]) + "\t" + str(Bphys[i]) + "\n")

fo.close()
