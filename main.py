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
	latticedim = iParams.latticename
    	filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p_over_2p.txt'
    	matel_3p_2p = rf.read_results(filename,5,len(msea),len(ms),1,nboots+1) #3pt/(2pt*2pt) ~ Bag
    	filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p.txt'
    	matel_3p=rf.read_results(filename,5,len(msea),len(ms),1,nboots+1) #3pt
    	filename = 'data/boot_matel_'+latticedim +'cubed_IW_3p_rat.txt'
   	matel_3p_rat = rf.read_results(filename,4,len(msea),len(ms),1,nboots+1) #3pt1/3pti ~ R
    	print "Done" + "\n"
	return matel_3p_2p, matel_3p, matel_3p_rat
	###########################################################################################


#########################################################################################################################################
def return_myMeson(iBetam,iParams,nboots):
	print " Reading Physical point Mass data..."
	#setup masses & lat
	ms = deepcopy(iParams.m_val_s)
	msea = deepcopy(iParams.m_sea_l)
	lat = iParams.latticename
	mval=np.hstack([msea,ms]) # list of light,  strange val quarks
	
	pion = 'l' + str(mval[0]) + '-l' + str(msea[0])
	kaon = 's' + str(mval[1]) + '-l'+ str(msea[0])

	#set up matrices
	b_MK = np.zeros([1,1,2,1,nboots+1])
	b_FK =  np.zeros([1,1,2,1,nboots+1])
	b_mpsq_fpsq = np.zeros([1,1,2,1,nboots+1])
	#set up filenames & read data
	filename=[]
	filename.append('Julia_data/mass/' + lat + '/mass-' + pion + '/mass-' + pion +'.dat')
	filename.append('Julia_data/mass/' + lat + '/mass-' + kaon + '/mass-' + kaon +'.dat')
	for ival in range(len(mval)):
		data , error = rf.read_bootstraps(filename[ival])
		#assign data to arrays
		for iboot in range(nboots): #Temporarily fill with the true phys data, before we have run and fit the lattice data
			b_MK[0,0,ival,0,iboot] = data[iboot]
			b_FK[0,0,ival,0,iboot] = 0.15551 + (randrange(0,100))/10000.0
			b_mpsq_fpsq[0,0,ival,0,iboot] = pow(b_MK[0,0,ival,0,iboot]/(4*pi*b_FK[0,0,ival,0,iboot]),2)
		b_MK[0,0,ival,0,-1] = data[-1]
		b_FK[0,0,ival,0,-1] = 0.13
		print b_MK[0,0,0,0,-1], b_FK[0,0,0,0,-1] 
		b_mpsq_fpsq[0,0,ival,0,-1] = pow(b_MK[0,0,ival,0,-1]/(4*pi*b_FK[0,0,ival,0,-1]),2)
	print "Done" + "\n"
	return b_MK, b_FK, b_mpsq_fpsq
############################################################################################################################################


###################################################################################################################################################################
def return_myBag(iBeta,iParams,nboots):
	print " Reading Physical point Bag data..."
	ms  = deepcopy(iParams.m_val_s)
	msea = deepcopy(iParams.m_sea_l)
	kaon = 's' + str(ms[0]) + '-l'+ str(msea[0])
	lat = iParams.latticename
	Ninv = [3.0/8,3.0/4,-0.5,3.0/5,1.0] #removes scaling to be consistent with Nicolas, re-added in later.
	B = np.zeros([5,1,1,1,nboots+1])
	R = np.zeros([5,1,1,1,nboots+1])
	#loop throught the channels
	for ic in [1,2, 3, 4, 5]:
		# 3 & 4 need to switched - convention choice
		if ic == 3:			
			filename = 'Julia_data/bag/' + lat  + '/bag-' + kaon + '/dt40/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) + '-channel'+str(ic+1)+'.dat'
		elif ic ==4:
			filename = 'Julia_data/bag/' + lat  + '/bag-' + kaon + '/dt40/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) + '-channel'+str(ic-1)+'.dat'
		else:
			filename = 'Julia_data/bag/' + lat  + '/bag-' + kaon + '/dt40/bag-s'+ str(ms[0]) + '-l' + str(msea[0]) + '-channel'+str(ic)+'.dat'
		data,error=rf.read_bootstraps(filename)
		for iboot in range(nboots+1):
			if ic == 1:
				B[ic-1,0,0,0,iboot] = data[iboot]/Ninv[ic-1]
				R[ic-1,0,0,0,iboot] = 1.0
			else:
				B[ic-1,0,0,0,iboot] = data[iboot]/Ninv[ic-1]
				R[ic-1,0,0,0,iboot] = data[iboot]/Ninv[ic-1]

	print "Done"+"\n"
	return B, R
##################################################################################################################################################################

def return_trueM(aM,a):
	dim_aM = np.shape(aM)

def main(scheme,kin,bas):

	latticedim = ['24','32','48','48','64']
	volname = ['24cubed','32cubed','48cubed','48cubedfine','64cubed']
	beta = [0,1]
	amp = ['/MK_4CHAN_','/FK_4CHAN_']
	#conv = 2 #l-s & s-s
	name_kin = ['E','gg','qq','qg','gq']
	name_scheme = ['MOM','ms']
	name_basis = ['Lattice', 'SUSY']


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
	#kin=0
	#ikin=0
	#scheme=0
	#ibas=0
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
	mp_phys = 0.13957018

	for iBeta in range(2): #loop through 24, and 32 (Non-Physical)
		
		#class cotaining paramaters set up    
		iParams = setup.Params(iBeta)
		params.append(iParams)   
		ainv.append(deepcopy(iParams.ainv))	

		b_FK, b_MK, b_mpsq_fpsq = return_Jamie_data(iBeta,iParams,nboot) #read and return meson mass + decay
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
		IB = lf.R_Interp(Bren_phys[:,:,:,:,:], ms[:], ms_phys)
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
		B_new, R_new = return_myBag(iBeta,iParams,nboot)
		b_MK, b_FK, b_mpsq_fpsq = return_myMeson(iBeta,iParams,nboot)
		#read all the new physical files here		
		#still need the R data + pion mass + decays
		
		##########################################################################################################################
		#Renormlise B and R
		Zfilename = "data/Z" + name_scheme[scheme] + "_boot_mu_match_" + volname[iBeta-2] + "_" + name_kin[kin] + "_block.out" 
		Rren_phys, Bren_phys = Renormalise(ibas,kin,iBeta,Zfilename,R=R_new,B=B_new)
		##########################################################################################################################

		#no need to interpolate here, already at the physical strange

		Blst.append(Bren_phys)
		Rlst.append(R)
		IRlst.append(Rren_phys)
		#IBlst.append(B_new)
		IBlst.append(Bren_phys)
		mlst.append(b_MK)
		flst.append(b_FK)
		m2_f2lst.append(b_mpsq_fpsq)
		
	p0 = [1,1,1]



	#p, cov, bp, bcov, bchisq = nlf.bootglobalfit(y1[0,0,0],y2[0,0,0],1/ainv[0],1/ainv[1],m1[0,0,0],m2[0,0,0],500,p0,mpsq_fpsq_phys)
	#actually do the global fit (looping through the sea masses)


	dimIR = np.shape(IR)
	dimIB = np.shape(IB)
	Rphys=[]
	Bphys=[]

	a_store = np.zeros([2])
	a_store[0] = 1/(ainv[2]*ainv[2])
	a_store[1] = 1/(ainv[3]*ainv[3])


	for ic in range(dimIB[0]): #loop channels
		IB3_store = np.zeros([2,nboot+1])
		m3_store = np.zeros([2,nboot+1])
		IR3_store = np.zeros([2,nboot+1])
		for iboot in range(nboot+1):
			IB3_store[0,iboot] = (IBlst[2])[ic,0,0,0,iboot]
			IB3_store[1,iboot] = (IBlst[3])[ic,0,0,0,iboot]
			m3_store[0,iboot] = pow((mlst[2])[0,0,0,0,iboot]*ainv[2],2)
			m3_store[1,iboot] = pow((mlst[3])[0,0,0,0,iboot]*ainv[3],2)
			IR3_store[0,iboot] = (IRlst[2])[ic,0,0,0,iboot]
			IR3_store[1,iboot] = (IRlst[3])[ic,0,0,0,iboot]
		for iv2 in range(dimIB[2]): #loop valence 2 (=1 here, unitary=sea)
			IB1_store = np.zeros([dimIB[1]-1,dimIB[3]]) 
			IB2_store = np.zeros([dimIB[1],dimIB[3]])
			IR1_store = np.zeros([dimIB[1]-1,dimIB[3]])
			IR2_store = np.zeros([dimIB[1],dimIB[3]])
			m1_store = np.zeros([dimIB[1]-1,dimIB[3]])  
			m2_store = np.zeros([dimIB[1],dimIB[3]])
			for iboot in range(dimIB[3]): #loop boots
				for isea in range(dimIB[1]): #loops sea quark
					if isea <= 1:
						IB1_store[isea,iboot] = (IBlst[0])[ic,isea,iv2,iboot]
						IR1_store[isea,iboot] = (IRlst[0])[ic,isea,iv2,iboot]
						m1_store[isea,iboot]=pow((mlst[0])[0,isea,-1,iv2,iboot]*ainv[0],2)
					IB2_store[isea,iboot] = (IBlst[1])[ic,isea,iv2,iboot]
					IR2_store[isea,iboot] = (IRlst[1])[ic,isea,iv2,iboot]
					m2_store[isea,iboot]=pow((mlst[1])[0,isea,-1,iv2,iboot]*ainv[1],2)
			nameR="./plots/R" + str(ic) + "_"+ name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]
			nameB="./plots/B"+str(ic)+"_"+name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]
			print m1_store[:,-1], m2_store[:,-1], m3_store[:,-1], pow(mp_phys,2)
			print IR1_store[:,-1], IR2_store[:,-1], IR3_store[:,-1]
			Rpcent, covcent, Rbp, Rbcov, Rbchisq, Ryp = nlf.bootglobalfit(IR1_store,IR2_store,1/(ainv[0]*ainv[0]),1/(ainv[1]*ainv[1]),m1_store,m2_store,500,p0,pow(mp_phys,2),r"$m^2_{\pi}$",r"$R_{"+str(ic+1)+"}$",nameR,IR3_store,a_store,m3_store)
			Bpcent, Bcovcent, Bbp, Bbcov, Bbchisq, Byp = nlf.bootglobalfit(IB1_store,IB2_store,1/(ainv[0]*ainv[0]),1/(ainv[1]*ainv[1]),m1_store,m2_store,500,p0,pow(mp_phys,2),r"$m^2_{\pi}$",r"$B_{"+str(ic+1)+"}$",nameB,IB3_store,a_store,m3_store)
			Rphys.append(Ryp)
			Bphys.append(Byp)

	fo = open("./results/R_B_"+name_basis[ibas]+"_"+name_scheme[scheme]+"_"+name_kin[kin]+".dat",'w')
	for i in range(5):
		j=i+1
		print str(j) + "\t" + str((Rphys[i])[-1]) + "\t"  + str(np.std(Rphys[i][0:-1])) + "\t" +  str(Bphys[i][-1]) + "\t" + str(np.std(Bphys[i][0:-1])) + "\n"

		fo.write(str(j) + "\t" + str(Rphys[i][-1]) + "\t"  + str(np.std(Rphys[i][0:-1])) + "\t" +  str(Bphys[i][-1]) + "\t" + str(np.std(Bphys[i][0:-1])) + "\n")

	fo.close()

if __name__ == "__main__":
	
	for ischeme in range(2):
		for ikin in range(5):
			for ibas in range(2):
				main(ischeme,ikin,ibas)
