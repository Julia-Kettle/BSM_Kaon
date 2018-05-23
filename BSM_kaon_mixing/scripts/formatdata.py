from copy import deepcopy
import utils
import file_io
import numpy as np
from lattice import Lattice
import sys

def calcR(matel_3p_rat,bdecay,bmass):
    """
    Calculates the ratios defined as R=(fK/mK)^2_exp * (m_K/f_K)^2_lat * <Oi>/<O1>
    """

    Mkphys = 0.5*(493.677+497.614)*pow(10,-3)
    Fkphys = 156.2*pow(10,-3)

    R=np.empty_like(matel_3p_rat)

    #only the kaon data needed - discard pion
    bmK = bmass[:,:-1,:]
    bfK = bdecay[:,:-1,:]

    for ich in range(np.size(R,0)):
        R[ich] = pow(Fkphys/Mkphys,2)*pow(bmK/bfK,2)*matel_3p_rat[ich]
    return R

def numpysave(filename,array):
    #save array to binary file using numpy
    fo=open(filename,'w')
    np.save(fo,array)
    fo.close()

def save_olddata(lattice):
    """
    Reads in the 3 point functions into suitably sized arrays
    """
    msea = lattice.m_sea_l
    ms = lattice.m_val_s
    mval = np.hstack([msea,ms])
    nboots=500

    filename = '../data/boot_matel_'+lattice.name +'cubed_IW_3p_over_2p.txt'
    matel_3p_2p = np.zeros([len(msea),len(ms),5,nboots+1])
    file_io.read_file_array(filename,matel_3p_2p) #3pt/(2pt*2pt) ~ Bag

    filename = '../data/boot_matel_'+lattice.name +'cubed_IW_3p.txt'
    matel_3p=np.empty_like(matel_3p_2p)
    file_io.read_file_array(filename,matel_3p) #3pt

    filename = '../data/boot_matel_'+lattice.name +'cubed_IW_3p_rat.txt'
    matel_3p_rat = np.ones_like(matel_3p_2p)
    #fill raw ratios array with data, leaving channel 0 as one. 
    file_io.read_file_array(filename,matel_3p_rat[:,:,1:,:]) #3pt1/3pti ~ R

    b_F = np.zeros([len(msea),len(ms)+1,nboots+1])
    b_M = np.zeros_like(b_F)
    b_msq_fsq = np.zeros_like(b_F)

    filend_basis = '_9.22' if lattice.name=='24' else '_12.52'
    for i in range(len(msea)):
            mv2 = np.hstack((ms,msea[i])) # 2nd valence, either s or l
            for j in range(len(mv2)):
                filenameF = '../data_Jamie/'+lattice.name + '/FK_4CHAN_mu' +utils.str5sf(msea[i])+'_mq'+utils.str5sf(msea[i])+'_ms'+utils.str5sf(mv2[j]) + filend_basis
                filenameM = '../data_Jamie/'+lattice.name + '/MK_4CHAN_mu' +utils.str5sf(msea[i])+'_mq'+utils.str5sf(msea[i])+'_ms'+utils.str5sf(mv2[j]) + filend_basis
                n, nboot, dataF = file_io.read_results_Jamie(filenameF)
                n, nboot, dataM = file_io.read_results_Jamie(filenameM)
                for l in range(len(dataM)):
                    b_F[i,j,l] = dataF[l]
                    b_M[i,j,l] = dataM[l]
                    b_msq_fsq[i,j,l]= pow((b_M[i,j,l]/(4*np.pi*b_F[i,j,l])),2)
    R=calcR(matel_3p_rat,b_F,b_M)
    try:
        numpysave("../unrenormalised/B_"+lattice.name+".bin",matel_3p_2p)
    except:
        pass
    numpysave("../unrenormalised/R_"+lattice.name+".bin",R)
    numpysave("../unrenormalised/m_"+lattice.name+".bin",b_M)
    numpysave("../unrenormalised/f_"+lattice.name+".bin",b_F)
    numpysave("../unrenormalised/m_4pif_sq_"+lattice.name+".bin",b_msq_fsq)


def save_mydata(lattice,latdir=""):
    """
    Reads in the new physical point bags and bare ratios
    """
    print "-------------------------------------------"
    print lattice.name
    print "-------------------------------------------"

    ms  = deepcopy(lattice.m_val_s)
    msea = deepcopy(lattice.m_sea_l)
    mval = np.hstack([msea,ms])

    lat = lattice.name+'cubed' if lattice.smeared ==False and lattice.fine  == False else lattice.name
    Ninv = [3.0/8,3.0/4,-0.5,3.0/5,1.0] #removes scaling to be consistent with Nicolas, re-added in later.
    #Ninv = [1.0,1.0,1.0,1.0,1.0]
    nboots=500

    matel_3p_2p = np.zeros([5,len(msea),1,nboots+1])
    matel_3p_rat = np.ones_like(matel_3p_2p)

    b_M=np.zeros([len(msea),2,nboots+1])
    b_F=np.zeros_like(b_M)
    b_msq_fsq=np.zeros_like(b_M)

    for isea in range(len(msea)):
        #pion and kaon name
        mesons = ['s'+str(ms[0])+'-l'+str(msea[isea]),'l'+str(msea[isea])+'-l'+str(msea[isea])]

        #Read the 3pt data looping through the channels 
        for ic in range(5):
            # 3 & 4 need to switched - convention choice
            if ic == 2:     chan=str(2*ic+2)
            elif ic ==3:    chan=str(2*ic-2)
            else:           chan=str(2*ic)
            filenameR='/Users/s1035546/FIT/'+latdir+'/fits/final/binned/uncorrelated/3ptrat/channel'+chan+'/3ptrat-'+mesons[0]+'/3ptrat-'+ mesons[0] +'.txt'
            filenameB='/Users/s1035546/FIT/'+latdir+'/fits/final/binned/uncorrelated/bag/channel'+chan+'/bag-'+mesons[0]+'/bag-'+ mesons[0] +'.txt'
            #filenameR='/Users/s1035546/Fits/' + lat  + '/3ptrat/channel' +chan + '/3ptrat-' + mesons[0] + '/3ptrat-'+ mesons[0] +'_boots.dat'
            #filenameB='/Users/s1035546/Fits/' + lat  + '/bag/channel' + chan + '/bag-' + mesons[0] + '/bag-' + mesons[0] +'_boots.dat'
            #Read bag parameters and remove renormalisation - matches Nicolas' convention.
            file_io.read_file_array(filenameB,matel_3p_2p[ic,isea])
            matel_3p_2p[ic,isea] = matel_3p_2p[ic,isea]/Ninv[ic]
            print matel_3p_2p[ic,isea,:,-1]
            #Read ratio parameters for the BSM leaving SM channel 1 as 1. 
            if ic > 0:
                file_io.read_file_array(filenameR,matel_3p_rat[ic,isea])
        #read pion and kaon masses and decay constants
        for imeson in range(len(mesons)):
            #set up filenames & read data for mass and decay constant
            filenamef='/Users/s1035546/FIT/'+latdir+'/fits/final/binned/uncorrelated/sim_mass/sim_mass_'+mesons[imeson]+'/decay_'+ mesons[imeson] +'.txt'
            filenamem='/Users/s1035546/FIT/'+latdir+'/fits/final/binned/uncorrelated/sim_mass/sim_mass_'+mesons[imeson]+'/sim_mass_'+ mesons[imeson] +'.txt'
            file_io.read_file_array(filenamem,b_M[isea,imeson])
            file_io.read_file_array(filenamef,b_F[isea,imeson])
            print b_M[:,:,-1]
            print b_F[:,:,-1]
            #calculate ratio of m/f for x axis data
            for iboot in range(nboots+1):
                b_msq_fsq[isea,imeson,iboot] = pow((b_M[isea,imeson,iboot]/(4*np.pi*b_F[isea,imeson,iboot])),2)

    R=calcR(matel_3p_rat,b_F,b_M)
    numpysave("../unrenormalised/B_"+lattice.name+".bin",matel_3p_2p)
    numpysave("../unrenormalised/R_"+lattice.name+".bin",R)
    numpysave("../unrenormalised/m_"+lattice.name+".bin",b_M)
    numpysave("../unrenormalised/f_"+lattice.name+".bin",b_F)
    numpysave("../unrenormalised/m_4pif_sq_"+lattice.name+".bin",b_msq_fsq)

    for i in range(501):
        print R[:,:,:,i]


    print R[:,:,-1]
    print matel_3p_2p[:,:,-1]
    print b_M[:,:,-1]
    print b_F[:,:,-1]

def main():
    lat48finesmeared=Lattice(48,True,fine=True)
    save_mydata(lat48finesmeared,'48cubedfinesmeared')
    lat24=Lattice(24,False)
    save_mydata(lat24,'24cubed')
    print "new 24 done"
    #lat32=Lattice(32,False)
    #save_mydata(lat32)
    print "new 32 done"
    #lat32=Lattice(32,False,True)
    #save_olddata(lat32)
    print "old 32 done"
    #lat24=Lattice(24,False,True)
    #save_olddata(lat24)
    print "old 24 done"
    print "48 done"
    lat48=Lattice(48,True)
    save_mydata(lat48,'48cubedsmeared')
    print "48 smeared done"
    lat64=Lattice(64,True)
    save_mydata(lat64,'64cubedsmeared')
    #print "64 done"
    lat32smeared=Lattice(32,True)
    save_mydata(lat32smeared,'32cubedsmeared')

if __name__ == "__main__":
    main()
