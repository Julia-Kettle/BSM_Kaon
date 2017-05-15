from copy import deepcopy
import utils
import file_io
import numpy as np
from lattice import Lattice

def calcR(matel_3p_rat,b_F,b_M):
    '''
    Calculates the  matrix (5,ml,ms,boots) of the ratios
    '''
    dim = []
    dim1 = np.shape(matel_3p_rat)
    dim2 = np.shape(b_F)
    for i in range(len(dim1)):
        dim.append(min(dim1[i],dim2[i]))
    dim[0] = 5
    R = np.zeros(dim)
    Mkphys = 0.5*(493.677+497.614)*pow(10,-3)
    Fkphys = 156.2*pow(10,-3)
    for index, val in np.ndenumerate(R):
        if index[0] == 0:
            R[index] = pow((Fkphys/Mkphys),2)*pow((b_M[index])/b_F[index],2)
        else:
            R[index] = pow((Fkphys/Mkphys),2)*pow((b_M[(0,)+index[1:]])/b_F[(0,)+index[1:]],2)*matel_3p_rat[(index[0]-1,) + index[1:]]
    return R

def save_olddata(lattice):
    """
    Reads in the 3 point functions into suitably sized arrays
    """
    msea = lattice.m_sea_l
    ms = lattice.m_val_s
    mval = np.hstack([msea,ms])
    latticedim = lattice.name
    nboots=500

    filename = '../data/boot_matel_'+latticedim +'cubed_IW_3p_over_2p.txt'
    matel_3p_2p = np.zeros([len(msea),len(ms),5,nboots+1])
    file_io.read_file_array(filename,matel_3p_2p) #3pt/(2pt*2pt) ~ Bag
    matel_3p_2p = np.swapaxes(matel_3p_2p,0,2)
    matel_3p_2p = np.swapaxes(matel_3p_2p,1,2)

    filename = '../data/boot_matel_'+latticedim +'cubed_IW_3p.txt'
    matel_3p=np.zeros([len(msea),len(ms),5,nboots+1])
    file_io.read_file_array(filename,matel_3p) #3pt
    matel_3p = np.swapaxes(matel_3p,0,2)
    matel_3p = np.swapaxes(matel_3p,1,2)

    filename = '../data/boot_matel_'+latticedim +'cubed_IW_3p_rat.txt'
    matel_3p_rat = np.zeros([len(msea),len(ms),4,nboots+1])
    file_io.read_file_array(filename,matel_3p_rat) #3pt1/3pti ~ R
    matel_3p_rat = np.swapaxes(matel_3p_rat,0,2)
    matel_3p_rat = np.swapaxes(matel_3p_rat,1,2)

    b_F = np.zeros([1,len(msea),len(ms)+1,nboots+1])
    b_M = np.zeros([1,len(msea),len(ms)+1,nboots+1])
    b_msq_fsq = np.zeros([1,len(msea),len(ms)+1,nboots+1])

    filend_basis = '_9.22' if lattice.name=='24' else '_12.52'
    for i in range(len(msea)):
            mv2 = np.hstack((ms,msea[i])) # 2nd valence, either s or l
            for j in range(len(mv2)):
                filenameF = '../data_Jamie/'+lattice.name + '/FK_4CHAN_mu' +utils.str5sf(msea[i])+'_mq'+utils.str5sf(msea[i])+'_ms'+utils.str5sf(mv2[j]) + filend_basis
                filenameM = '../data_Jamie/'+lattice.name + '/MK_4CHAN_mu' +utils.str5sf(msea[i])+'_mq'+utils.str5sf(msea[i])+'_ms'+utils.str5sf(mv2[j]) + filend_basis
                n, nboot, dataF = file_io.read_results_Jamie(filenameF)
                n, nboot, dataM = file_io.read_results_Jamie(filenameM)
                for l in range(len(dataM)):
                    b_F[0,i,j,l] = dataF[l]
                    b_M[0,i,j,l] = dataM[l]
                    b_msq_fsq[0,i,j,l]= pow((b_M[0,i,j,l]/(4*np.pi*b_F[0,i,j,l])),2)
    print b_F[0,:,:,-1]
    print b_M[0,:,:,-1]
    R=calcR(matel_3p_rat,b_F,b_M)

    bag_out = "../unrenormalised/B_"+latticedim+"old.bin"
    ratio_out="../unrenormalised/R_"+latticedim+"old.bin"
    mass_out="../unrenormalised/m_"+latticedim+"old.bin"
    decay_out="../unrenormalised/f_"+latticedim+"old.bin"
    mass_decay_out="../unrenormalised/m_4pif_sq"+latticedim+"old.bin"
    fo=open(bag_out,'w')
    np.save(fo,matel_3p_2p)
    fo.close()

    fo=open(ratio_out,'w')
    np.save(fo,R)
    fo.close()

    fo=open(mass_out,'w')
    np.save(fo,b_M)
    fo.close()

    fo=open(decay_out,'w')
    np.save(fo,b_F)
    fo.close()

    fo=open(mass_decay_out,'w')
    np.save(fo,b_msq_fsq)
    fo.close()

def save_mydata(lattice,dt):
    """
    Reads in the new physical point bags and bare ratios
    """
    print lattice.name
    ms  = deepcopy(lattice.m_val_s)
    msea = deepcopy(lattice.m_sea_l)
    print len(msea)
    mval = np.hstack([msea,ms])
    #lat = lattice.lattice_name+["cubed","cubed","smeared","cubed"][self.ibeta]
    lat = lattice.name+'cubed' if lattice.smeared ==False and lattice.fine  == False else lattice.name
    Ninv = [3.0/8,3.0/4,-0.5,3.0/5,1.0] #removes scaling to be consistent with Nicolas, re-added in later.
    nboots=500
    matel_3p_2p = np.zeros([5,len(msea),1,nboots+1])
    matel_3p_rat = np.zeros([4,len(msea),1,nboots+1])
    for isea in range(len(msea)):
        kaon = 's' + str(ms[0]) + '-l'+ str(msea[isea])
        for ic in range(5):
            try:
                # 3 & 4 need to switched - convention choice
                if ic == 2:
                    filenameR='/Users/s1035546/Fits/' + lat  + '/3ptrat/channel' +str(2*ic+2) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                    filenameB='/Users/s1035546/Fits/' + lat  + '/bag/channel' +str(2*ic+2) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                elif ic ==3:
                    filenameR='/Users/s1035546/Fits/' + lat  + '/3ptrat/channel' +str(2*ic-2) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                    filenameB='/Users/s1035546/Fits/' + lat  + '/bag/channel' +str(2*ic-2) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                else:
                    filenameR='/Users/s1035546/Fits/' + lat  + '/3ptrat/channel' +str(2*ic) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                    filenameB='/Users/s1035546/Fits/' + lat  + '/bag/channel' +str(2*ic) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'

                file_io.read_file_array(filenameB,matel_3p_2p[ic,isea])
                if ic > 0:
                    file_io.read_file_array(filenameR,matel_3p_rat[ic-1,isea])
            except:
                dt=26
                if ic == 2:
                    filenameR='/Users/s1035546/Fits/' + lat  + '/3ptrat/channel' +str(2*ic+2) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                    filenameB='/Users/s1035546/Fits/' + lat  + '/bag/channel' +str(2*ic+2) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                elif ic ==3:
                    filenameR='/Users/s1035546/Fits/' + lat  + '/3ptrat/channel' +str(2*ic-2) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                    filenameB='/Users/s1035546/Fits/' + lat  + '/bag/channel' +str(2*ic-2) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                else:
                    filenameR='/Users/s1035546/Fits/' + lat  + '/3ptrat/channel' +str(2*ic) + '/3ptrat-' + kaon + '-dt'+str(dt)+'/3ptrat-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'
                    filenameB='/Users/s1035546/Fits/' + lat  + '/bag/channel' +str(2*ic) + '/bag-' + kaon + '-dt'+str(dt)+'/bag-s'+ str(ms[0]) + '-l' + str(msea[isea]) +'-dt'+str(dt)+'_boots.dat'

                file_io.read_file_array(filenameB,matel_3p_2p[ic,isea])
                if ic > 0:
                    file_io.read_file_array(filenameR,matel_3p_rat[ic-1,isea])
            for iboot in range(nboots+1):
                matel_3p_2p[ic,isea,:,iboot] = matel_3p_2p[ic,isea,:,iboot]/Ninv[ic]
    
    b_M=np.zeros([1,len(msea),2,nboots+1])
    b_F=np.zeros([1,len(msea),2,nboots+1])
    b_msq_fsq=np.zeros([1,len(msea),2,nboots+1])
    for isea in range(len(msea)):
        kaon = 's' + str(ms[0]) + '-l'+ str(msea[isea])
        print msea[isea]
        if lat == '24smeared':
            kaon = 'l'+ str(msea[isea]) + '-s' + str(ms[0])
        for imv2 in range(2):
            if imv2 == 1:
                meson = 'l' + str(msea[isea]) + '-l' + str(msea[isea])
            elif imv2 == 0:
                meson = kaon
            print meson
            #set up matrices
            #set up filenames & read data
            filenameF='/Users/s1035546/Fits/' + lat + '/mass/mass-' + meson + '/decay-' + meson +'_boots.dat'
            filenameM='/Users/s1035546/Fits/' + lat + '/mass/mass-' + meson + '/mass-' + meson +'_boots.dat'
            #filenameF='/Users/s1035546/Fits/' + lat + '/sim_mass/sim_mass-' + meson + '/decay-' + meson +'_boots.dat'
            #filenameM='/Users/s1035546/Fits/' + lat + '/sim_mass/sim_mass-' + meson + '/sim_mass-' + meson +'_boots.dat'
            print filenameM
            file_io.read_file_array(filenameM,b_M[0,isea,imv2,:])
            file_io.read_file_array(filenameF,b_F[0,isea,imv2,:])
            for iboot in range(nboots+1):
                b_msq_fsq[0,isea,imv2,iboot] = pow((b_M[0,isea,imv2,iboot]/(4*np.pi*b_F[0,isea,imv2,iboot])),2)

    R=calcR(matel_3p_rat,b_F,b_M)
    print "R - ", R[:,:,-1]
    bag_out = "../unrenormalised/B_"+lattice.name+".bin"
    ratio_out="../unrenormalised/R_"+lattice.name+".bin"
    mass_out="../unrenormalised/m_"+lattice.name+".bin"
    decay_out="../unrenormalised/f_"+lattice.name+".bin"
    mass_decay_out="../unrenormalised/m_4pif_sq"+lattice.name+".bin"

    fo=open(bag_out,'w')
    np.save(fo,matel_3p_2p)
    fo.close()

    fo=open(ratio_out,'w')
    np.save(fo,R)
    fo.close()

    fo=open(mass_out,'w')
    np.save(fo,b_M)
    fo.close()

    fo=open(decay_out,'w')
    np.save(fo,b_F)
    fo.close()

    fo=open(mass_decay_out,'w')
    np.save(fo,b_msq_fsq)
    fo.close()
'''
lat24=Lattice(24,False)
save_mydata(lat24,24)
print "new 24 done"
lat24=Lattice(24,True)
save_mydata(lat24,32)
print "new 24 smeared done"
lat32=Lattice(32,False)
save_mydata(lat32,32)
print "new 32 done"
lat32=Lattice(32,False,True)
save_olddata(lat32)
print "old 32 done"
lat24=Lattice(24,False,True)
save_olddata(lat24)
print "old 24 done"
'''
#lat48=Lattice(48,False)
#save_mydata(lat48,40)
print "48 done"
lat48=Lattice(48,True)
save_mydata(lat48,40)
#lat48fine=Lattice(48,False,fine=True)
#save_mydata(lat48fine,40)
print "48 smeared done"
#lat64=Lattice(64,False)
#save_mydata(lat64,52)
print "64 done"
