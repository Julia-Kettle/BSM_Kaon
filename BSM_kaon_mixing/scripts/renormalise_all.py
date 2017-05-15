from renorm_class import Renormalisation
from renorm_class import StrangeAdjustment
from global_fit import GlobalFit
import numpy as np
import sys
import utils
from lattice import Lattice


def main(lat,smeared,ikin,ibas,ischeme,dR=None,dB=None,old=False,fine=False):

    lattice=Lattice(lat,smeared,old=old,fine=fine)
    lat = lattice.name #if smeared == True else lattice.name+"cubed"
    print "-------------------------------------"
    print lattice.name
    print "-------------------------------------"
    if old == False:
        m = np.load("../unrenormalised/m_"+str(lat)+".bin")
        f = np.load("../unrenormalised/f_"+str(lat)+".bin")
        m_4pif_sq = np.load("../unrenormalised/m_4pif_sq"+str(lat)+".bin")
        B = np.load("../unrenormalised/B_"+str(lat)+".bin")
        R = np.load("../unrenormalised/R_"+str(lat)+".bin")
    else:
        m = np.load("../unrenormalised/m_"+str(lat)+"old.bin")
        f = np.load("../unrenormalised/f_"+str(lat)+"old.bin")
        m_4pif_sq = np.load("../unrenormalised/m_4pif_sq"+str(lat)+"old.bin")
        B = np.load("../unrenormalised/B_"+str(lat)+"old.bin")
        R = np.load("../unrenormalised/R_"+str(lat)+"old.bin")

    print "B - ",B[:,:,:,-1]
    print "R - ",R[:,:,:,-1]
    renorm=Renormalisation(lattice,ischeme,ikin,ibas)
    R=renorm.Renorm_R(R)
    B=renorm.Renorm_B(B)

    #bfo = open("renormalised/B_"+str(lat)+".bin",'w')
    #np.save(bfo,B)
    #bfo.close()
    #rfo = open("renormalised/R_"+str(lat)+".bin",'w')
    #np.save(rfo,R)
    #rfo.close()



    if len(lattice.m_val_s)==1:
        if lattice.m_val_s != lattice.m_s_phys:
            if dB == None or dR ==None:
                print "Error - need slope"
                sys.exit(-1)
            else:
                print "adjusting the phys pt data"
                sa=StrangeAdjustment(lattice)
                R = sa.extrap_1pt(R,dR)
                B = sa.extrap_1pt(B,dB)

        else:
            B=B[:,:,0,:]
            R=R[:,:,0,:]
            '''
            bfo = open("renormalised/B_int"+str(lat)+".bin",'w')
            np.save(bfo,B)
            bfo.close()
            rfo = open("renormalised/R_int"+str(lat)+".bin",'w')
            np.save(rfo,R)
            rfo.close()
            '''
    else:
        #still need to interpolate
        sa=StrangeAdjustment(lattice)
        B,dB = sa.interpolation(B)
        R,dR = sa.interpolation(R)
        print dR[:,:,-1]
        #dB = dB[:,0,:]
        #dR = dR[:,0,:]
        dB = sa.slope_extrap(dB,lattice.m_l_phys)
        dR = sa.slope_extrap(dR,lattice.m_l_phys)

        '''
        bfo = open("renormalised/B_int"+str(lat)+".bin",'w')
        np.save(bfo,B)
        bfo.close()
        rfo = open("renormalised/R_int"+str(lat)+".bin",'w')
        np.save(rfo,R)
        rfo.close()
        '''
    
    for iml in range(len(lattice.m_sea_l)):
        if old == False:
            rfo = open("../renormalised/R_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+".bin",'w')
            bfo = open("../renormalised/B_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+".bin",'w')
            mfo = open("../renormalised/m_4pif_sq_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+".bin",'w')
        else:
            rfo = open("../renormalised/R_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_old.bin",'w')
            bfo = open("../renormalised/B_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_old.bin",'w')
            mfo = open("../renormalised/m_4pif_sq_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_old.bin",'w')

        np.save(bfo,B[:,iml,:])
        np.save(rfo,R[:,iml,:])
        np.save(mfo,m_4pif_sq[0,iml,-1,:])


    return B,R,dB,dR

    '''
    for ml in lattice.m_sea_l:
    
    
    bfo = open("renormalised/B_"+str(lat)+"cubed.bin",'w')
    rfo = open("renormalised/R_"+str(lat)+"cubed.bin",'w')
    np.save(bfo,B)
    np.save(rfo,R)

    bfo = open("renormalised/B_int"+str(lat)+"cubed.bin",'w')
    rfo = open("renormalised/R_int"+str(lat)+"cubed.bin",'w')
    np.save(bfo,B[:,:,0,:])
    np.save(rfo,R[:,:,0,:])
    '''

    '''
        Old code for strange adjustment - already done for 64 - not needed elsewhere now
    sa = StrangeAdjustment(lattice)
    print lat
    if lat == 32:
        print "32 hit"
        R, dR = sa.interpolation(R)
        B, dB = sa.interpolation(B)
        dR = sa.slope_extrap(dR,0.000678)
        dB = sa.slope_extrap(dB,0.000678)
    elif lat == 64:
        R = sa.extrap_1pt(R,dR)
        B = sa.extrap_1pt(B,dB)
    else:
        R = R[:,:,0,:]
        B = B[:,:,0,:]
    
    bfo = open("renormalised/B_int"+str(lat)+"cubed.bin",'w')
    rfo = open("renormalised/R_int"+str(lat)+"cubed.bin",'w')
    np.save(bfo,B)
    np.save(rfo,R)
    
    if lat == 64:
        bfo = open("renormalised/B_int"+str(lat)+"cubed.bin",'r')
        rfo = open("renormalised/R_int"+str(lat)+"cubed.bin",'r')
        B=np.load(bfo)
        R=np.load(rfo)
    else:
        print str(lat) + str(ibas) + str(ikin) +str(ischeme)
        
        B = np.load("../unrenormalised/B_"+str(lat)+"cubed.bin")
        R = np.load("../unrenormalised/R_"+str(lat)+"cubed.bin")

        
        bfo = open("renormalised/B_"+str(lat)+"cubed.bin",'r')
        rfo = open("renormalised/R_"+str(lat)+"cubed.bin",'r')
        B=np.load(bfo)
        R=np.load(rfo)
        
    for i in range(len(lattice.m_sea_l)):
        betas.append(count)
        storeB.append(np.squeeze(B[:,i,:]))
        storeR.append(np.squeeze(R[:,i,:]))
        storem2.append(np.squeeze(m_4pif_sq[:,i,-1,:]))
        storea2.append(np.squeeze(lattice.a2))
    '''


for ibas in range(1):
    for ikin in range(1):
        for ischeme in range(1):
            ikin=1
            ibas=0
            ischeme=0


            main(24,False,ikin,ibas,ischeme)
            main(24,True,ikin,ibas,ischeme)
            main(32,False,ikin,ibas,ischeme)
            B,R,dB,dR=main(24,False,ikin,ibas,ischeme,old=True)
            main(48,False,ikin,ibas,ischeme,dB=dB,dR=dR)
            main(48,False,ikin,ibas,ischeme,fine=True)
            main(48,True,ikin,ibas,ischeme)
            B,R,dB,dR=main(32,False,ikin,ibas,ischeme,old=True)
            main(64,False,ikin,ibas,ischeme,dB=dB,dR=dR)
