import nonlinearfit
import numpy as np
import random
import readfiles as rf
import resampling as rsamp
import setup

"""
Jamie's data 
in dataJamie/(latticepts)/(Amptype)_mu(msealight)_mq(mql=msealight)_ms(mvalstrange)_(ibetatype)
latticepts = 24 or 32

"""

latticedim = ['24','32']
amp_par = ['/MK_4CHAN_','/FK_4CHAN_']
conv = 2 #l-s & s-s
name_kin = ['E','gg','qq','qg''gq']
name_scheme = ['MOM','ms']
name_basis = ['Lattice basis', 'SUSY basis']
filend_basis = ['_9.22', '_12.52']
#set up 2d matix of tex names later

for kin in range(len(name_kin)):
    for scheme in range(len(name_scheme)):
        for basis in range(len(name_basis)):
            params = []
            for lattice in range(len(latticedim)): # should I move this for loop to outside? params stay same during. Or do a separate loop to set up the params classes?

                iParams = setup.Params(latticedim[lattice])
                params.append(iParams)
                """
                Need to open Jamie's data here in data_Jamie/(24/32)/(amp_par)mu(mu)_mq(mq)_ms(ms)(filend_basis)
                24:  mu = 0.0050,0.0100,0.0200; mq = mu,0.0300,0.0350,0.0400; ms = 0.0300,0.0350,0.0400
                32:  mu = 0.0040,0.0060,0.0080; mq = 0.0040, 0.0060, 0.0080; ms = 0.0250,0.0300,mq
                read_results_Jamie(filename)
                the data returned is a*amp_par
                """
                mu = iParams.m_sea_l
                for i_mu in mu:
                    mq = iParams.m_val_s
                    mq.append(mu):
                    ms = iParams.m_val_s
                    
                #filename = 'data_Jamie/'+latticedim+amp_par+'mu'+str(ml)+'_mq'+str(mq)+'_ms'+str(ms) 
                
                
    #compute f_pi
