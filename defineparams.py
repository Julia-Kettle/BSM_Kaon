import numpy as np

def params(lattice):
    '''
    defines parameters depending on lattice type. 
    taken from Nicholas matlab code
    '''
    if lattice == '24cubedIW':
        m_sea_s = 0.04
        m_seas_s_phys = 0.03224
        m_sea_s_pq = 0.03
        m_sea_l = [0.005, 0.010, 0.02]
        m_val_s = [0.030, 0.035, 0.040]
        n_l_m = 1
        n_uni_pi = 2
        kaon_f = 2
        kaon_l = 4
        n_m_val_l = 1
        m_m_val_s = len(m_val_s)
        n_m_sea_l = len(m_sea_l)
        norm = pow(24,3)
        twall = 1
        cps_conv = 1 ##???
        n_chan = 5;
        
        #b-box, p-point
        name_chan= ['bPS_bPS','bPS_bA0','bPS_pPS','bPS_pA0','bA0_pA0']
        sgn_chan = [1,-1,1,-1,1]
        n_src = 2;
        name_src = ['sr','sk']
        
    
