import numpy as np

class Params:
    	'''
    	defines parameters depending on lattice type. 
    	taken from Nicholas matlab code
    	Need to figure out how to return these values easily...
    	Maybe put into main code??    
	'''
	def __init__(self, beta):

        	if beta == 0:
			self.latticename = '24'
            		self.ainv = 1.7845
			self.m_sea_s = 0.04
            		self.m_sea_s_phys = 0.03224
            		self.m_sea_s_pq = 0.03
            		self.m_sea_l = [0.0050, 0.0100, 0.0200]
            		self.m_val_s = [0.0300, 0.0350, 0.0400]
            		self.n_l_m = 1 
            		self.n_uni_pi = 2
            		self.kaon_f = 2
            		self.kaon_l = 4
            		self.n_m_val_l = 1
            		self.n_m_val_s = len(self.m_val_s)
            		self.n_m_sea_l = len(self.m_sea_l)
            		self.norm = pow(24,3)
            		self.twall = 1
            		self.cps_conv = 1 ##???
            		self.n_chan = 5;
        
        	elif beta == 1:
			self.latticename = '32'
			self.ainv = 2.3883
            		self.m_sea_s = 0.04
            		self.m_sea_s_phys = 0.02477
            		self.m_sea_s_pq = 0.025
            		self.m_sea_l = [0.0040, 0.0060, 0.0080]
            		self.m_val_s = [0.0250, 0.0300]
            		self.n_l_m = 10
            		self.n_uni_pi = 3
            		self.kaon_f = 11
           	 	self.kaon_l = 19
            		self.n_m_val_l = 1
            		self.m_m_val_s = len(self.m_val_s)
            		self.n_m_sea_l = len(self.m_sea_l)
            		self.norm = pow(32,3)
            		self.twall = 1
            		self.cps_conv = 3.2 ##???
            		self.n_chan = 5;


        	elif beta == 1:
			self.latticename = '48'
			self.ainv = 2.3883
            		self.m_sea_s = 0.0362
            		self.m_sea_s_phys = 0.02477
            		self.m_sea_s_pq = 0.025
            		self.m_sea_l = [0.00078]
            		self.m_val_s = [0.0362]
            		self.n_l_m = 10
            		self.n_uni_pi = 1
            		self.kaon_f = 11
           	 	self.kaon_l = 19
            		self.n_m_val_l = 1
            		self.m_m_val_s = len(self.m_val_s)
            		self.n_m_sea_l = len(self.m_sea_l)
            		self.norm = pow(32,3)
            		self.twall = 1
            		self.cps_conv = 3.2 ##???
            		self.n_chan = 5;
        	
		elif beta == 3:
			self.latticename = '48fine'
			self.ainv = 2.3883
            		self.m_sea_s = 0.02661
            		self.m_sea_s_phys = 0.02661
            		self.m_sea_s_pq = 0.02661
            		self.m_sea_l = [0.002661]
            		self.m_val_s = [0.02661]
            		self.n_l_m = 10
            		self.n_uni_pi = 3
            		self.kaon_f = 11
           	 	self.kaon_l = 19
            		self.n_m_val_l = 1
            		self.m_m_val_s = len(self.m_val_s)
            		self.n_m_sea_l = len(self.m_sea_l)
            		self.norm = pow(32,3)
            		self.twall = 1
            		self.cps_conv = 3.2 ##???
            		self.n_chan = 5;
        	
		elif beta == 4:
			self.latticename = '64'
			self.ainv = 2.3883
            		self.m_sea_s = 0.02114
            		self.m_sea_s_phys = 0.02114
            		self.m_sea_s_pq = 0.02114
            		self.m_sea_l = [0.000678]
            		self.m_val_s = [0.02114]
            		self.n_l_m = 10
            		self.n_uni_pi = 3
            		self.kaon_f = 11
           	 	self.kaon_l = 19
            		self.n_m_val_l = 1
            		self.m_m_val_s = len(self.m_val_s)
            		self.n_m_sea_l = len(self.m_sea_l)
            		self.norm = pow(32,3)
            		self.twall = 1
            		self.cps_conv = 3.2 ##???
            		self.n_chan = 5;




		self.name_chan= ['bPS_bPS','bPS_bA0','bPS_pPS','bPS_pA0','bA0_pA0']
		self.sgn_chan = [1,-1,1,-1,1]
       		self.n_src = 2;
		self.name_src = ['sr','sk']

	def show(self,values="all"):
        	variables = self.__dict__
        	if values == "all":
            		for key in variables:
                		print key + " - ", variables[key]
        	elif type(values) is list:
            		for i in values:
                		print i + " - ", variables[i]
        	elif type(values) is str:
            		print values + " - ", variables[values]
        	else:
            		print 'Error, wrong data types passes, must pass show "all", a string or a list of strings' 

 
 
def masses(m_val_s,m_sea_l, convention):

     '''
     set up indices and val q_m for given l seas q_m
     2 conventions - light-strange only
     & light-light + light-strange
     '''
     
     if convention==1:
         m_val1 = m_val_s
         m_val2 = m_sea_l #taking a light valence quark to have same mass as
         #sea quark
     if convention==2:
         m_val1 = np.hstack((m_sea_l,m_val_s))
         m_val2 = m_sea_l

     return m_val1, m_val2


def set_name_basis(i):
    if i==1:
        name_bas = 'Lattice basis'
    elif i==2:
        name_bas = 'SUSY basis'

    return name_bas
