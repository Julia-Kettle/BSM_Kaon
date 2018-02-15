from renorm_class import Renormalisation
from renorm_class import StrangeAdjustment
from global_fit import GlobalFit
import numpy as np
import sys
import utils
from lattice import Lattice


def main(lat,smeared,ikin,ibas,ischeme,dR=None,dB=None,old=False,fine=False,am48=False,ainv=1):

    lattice=Lattice(lat,smeared,old=old,fine=fine)
    lat = lattice.name #if smeared == True else lattice.name+"cubed"
    print "-------------------------------------"
    print lattice.name
    print "-------------------------------------"
    kin=['E','qq','gg'][ikin] 
    bas=['SUSY','Lattice'][ibas]
    scheme=['MOM','ms'][ischeme]
    try:
        if old == False:
            m = np.load("../unrenormalised/m_"+str(lat)+".bin")
            f = np.load("../unrenormalised/f_"+str(lat)+".bin")
            m_4pif_sq = np.load("../unrenormalised/m_4pif_sq_"+str(lat)+".bin")
            B = np.load("../unrenormalised/B_"+str(lat)+".bin")
            R = np.load("../unrenormalised/R_"+str(lat)+".bin")
        else:
            m = np.load("../unrenormalised/m_"+str(lat)+"old.bin")
            f = np.load("../unrenormalised/f_"+str(lat)+"old.bin")
            m_4pif_sq = np.load("../unrenormalised/m_4pif_sq_"+str(lat)+"old.bin")
            B = np.load("../unrenormalised/B_"+str(lat)+"old.bin")
            R = np.load("../unrenormalised/R_"+str(lat)+"old.bin")
    except:
        print "failed"
        return

    try:
        renorm=Renormalisation(lattice,ischeme,ikin,ibas)
    except:
        return
    R=renorm.renorm_R(R)
    B=renorm.renorm_B(B)

    if len(lattice.m_val_s)==1:
        if lattice.m_val_s != lattice.m_s_phys:
            if dB == None or dR ==None:
                print "Error - need slope"
                sys.exit(-1)
            else:
                print "adjusting the phys pt data"
                sa=StrangeAdjustment(lattice)
                print "R before - ", R[:,:,:,-1]
                R = sa.extrap_1pt(R,dR)
                B = sa.extrap_1pt(B,dB)
                print "R after -", R[:,:,-1]

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
        #dB = dB[:,0,:]
        #dR = dR[:,0,:]
        if am48 == True :
            lattice.m_l_phys[0] = 0.002144*2.774/ainv
        dB = sa.slope_extrap(dB,lattice.m_l_phys)
        dR = sa.slope_extrap(dR,lattice.m_l_phys)

    for iml in range(len(lattice.m_sea_l)):
        if old == False:
            rfo = open("../renormalised/R_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+".bin",'w')
            bfo = open("../renormalised/B_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+".bin",'w')
            mfo = open("../renormalised/m_4pif_sq_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+".bin",'w')
            ffo = open("../renormalised/f_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+".bin",'w')
        else:
            rfo = open("../renormalised/R_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+"_old.bin",'w')
            bfo = open("../renormalised/B_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+"_old.bin",'w')
            mfo = open("../renormalised/m_4pif_sq_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+".bin",'w')
            ffo = open("../renormalised/f_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+"_old.bin",'w')

        np.save(bfo,B[:,iml,:])
        np.save(rfo,R[:,iml,:])
        np.save(mfo,m_4pif_sq[iml,-1,:])
        np.save(ffo,f[iml,-1,:])

    return B,R,dB,dR

for ibas in range(2):
    for ikin in range(2):
        for ischeme in range(2):
            ikin=2
            ibas=0
            print ['E','qq','gg'][ikin], ['SUSY','Lattice'][ibas], ['MOM','ms'][ischeme]
            '''
            main(24,False,ikin,ibas,ischeme)
            main(32,False,ikin,ibas,ischeme)
            main(48,True,ikin,ibas,ischeme)
            main(64,True,ikin,ibas,ischeme)
            main(32,True,ikin,ibas,ischeme)
            main(48,False,ikin,ibas,ischeme,fine=True)
            '''
            main(48,True,ikin,ibas,ischeme,fine=True)
