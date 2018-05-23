from renorm_class import Renormalisation
from renorm_class import StrangeAdjustment
from global_fit import GlobalFit
import numpy as np
import sys
import utils
from lattice import Lattice


def getdata(lat,datadir,old=False):
    fileend=".bin" if not(old) else "_old.bin"
    m = np.load(datadir+"/m_"+str(lat)+fileend)
    f = np.load(datadir+"/f_"+str(lat)+fileend)
    m_4pif_sq = np.load(datadir+"/m_4pif_sq_"+str(lat)+fileend)
    B = np.load(datadir+"/B_"+str(lat)+fileend)
    R = np.load(datadir+"/R_"+str(lat)+fileend)

    return m,f,m_4pif_sq,B,R

def savetext(filename,array):
    fotext  =   open(filename,'w')
    for entry in np.nditer(array):
            fotext.write(str(entry)+"\n")
    fotext.close()

def savealltotext(savedir,fileend,R,B,m_4pif,f):
        savetext(savedir+"/textformat/R_"+fileend,R)
        savetext(savedir+"/textformat/B_"+fileend,B)
        savetext(savedir+"/textformat/m_4pif_"+fileend,m_4pif)
        savetext(savedir+"/textformat/f_"+fileend,f)

def savebinary(filename,array):
    fo  =   open(filename,'w')
    np.save(fo,array)
    fo.close()

def savealltobinary(savedir,fileend,R,B,m_4pif,f):        
        savebinary(savedir+"/R_"+fileend,R)      
        savebinary(savedir+"/B_"+fileend,B)      
        savebinary(savedir+"/m_4pif_"+fileend,m_4pif)      
        savebinary(savedir+"/f_"+fileend,f)

def strangeinterpolate(lattice,B,R,dB,dR):
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
    else:
        #still need to interpolate
        sa=StrangeAdjustment(lattice)
        B,dB = sa.interpolation(B)
        R,dR = sa.interpolation(R)
        if am48 == True :
            lattice.m_l_phys[0] = 0.002144*2.774/ainv
        dB = sa.slope_extrap(dB,lattice.m_l_phys)
        dR = sa.slope_extrap(dR,lattice.m_l_phys)
    return B,R,dB,dR

def main(lat,smeared,ikin,ibas,ischeme,dR=None,dB=None,old=False,fine=False,am48=False,ainv=1):

    lattice=Lattice(lat,smeared,old=old,fine=fine)
    lat = lattice.name #if smeared == True else lattice.name+"cubed"
    print "-------------------------------------"
    print lattice.name
    print "-------------------------------------"
    kin=['E','qq','gg'][ikin] 
    bas=['SUSY','Lattice'][ibas]
    boolSUSY=True if bas=='SUSY' else False
    scheme=['MOM','ms'][ischeme]

    datadir =   "/Users/s1035546/Global_Fitting_Environment/BSM_kaon_mixing/unrenormalised"
    baredir =   datadir
    dirRenorm="/Users/s1035546/Global_Fitting_Environment/BSM_kaon_mixing/data"
    dirZBK=dirRenorm
    outdir="/Users/s1035546/Global_Fitting_Environment/BSM_kaon_mixing/renormalised"

    #try:
    m,f,m_4pif_sq,B,R   =   getdata(lat,datadir,old)
    print "renormalising - "+scheme+" "+kin+" "+bas    
    #except:
    #print "failed"
    #return
    

    for iml in range(len(lattice.m_sea_l)):
        fileend=str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+".txt"
        savealltotext(baredir,fileend,R[:,iml,:],B[:,iml,:],m_4pif_sq[iml,-1,:],f[iml,-1,:])



    #read the renorm files
    zbk=Renormalisation.readZBK(lattice.renormname,scheme,kin,dirZBK)
    zs=Renormalisation.readZS(lattice.renormname,scheme,kin,dirRenorm)
    zv=Renormalisation.readZV(lattice.renormname,dirRenorm)

    B,R,dB,dR=strangeinterpolate(lattice,B,R,dB,dR)

    renorm = Renormalisation(zbk,zs,zv,boolSUSY)
    R=renorm.renorm_R(R)
    B=renorm.renorm_B(B)



    for iml in range(len(lattice.m_sea_l)):
        print "saving to ","../renormalised/B_"+str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_"+scheme+"_"+kin+".bin",'w'
        if not(old):
            fileendbin=str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+".bin"
            fileendtxt=str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+".txt"
        else:
            fileendbin=str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_old.bin"
            fileendtxt=str(lattice.name)+"ml_"+str(lattice.m_sea_l[iml])+"_ms_"+str(lattice.m_s_phys[0])+"_"+bas+"_old.txt"

        savealltotext(outdir,fileendtxt,R[:,iml,:],B[:,iml,:],m_4pif_sq[iml,-1,:],f[iml,-1,:]) 
        savealltobinary(outdir,fileendbin,R[:,iml,:],B[:,iml,:],m_4pif_sq[iml,-1,:],f[iml,-1,:]) 

    return dB,dR

for ibas in range(2):
    for ikin in range(1,3):
        for ischeme in range(1):
            #ikin=2
            #ibas=0
            print ['E','qq','gg'][ikin], ['SUSY','Lattice'][ibas], ['MOM','ms'][ischeme]
            main(24,False,ikin,ibas,ischeme)
            main(32,False,ikin,ibas,ischeme)
            main(48,True,ikin,ibas,ischeme)
            main(64,True,ikin,ibas,ischeme)
            main(32,True,ikin,ibas,ischeme)
            main(48,False,ikin,ibas,ischeme,fine=True)
            main(48,True,ikin,ibas,ischeme,fine=True)
