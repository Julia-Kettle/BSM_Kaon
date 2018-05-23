import utils
import file_io

class Lattice:
    def __init__(self,l,smeared,old=False,fine=False):
        self.dim = str(l)+"^3"
        self.smeared = smeared #boolean value
        self.fine = fine
        self.name = str(l)
        self.name = self.name+"fine" if fine==True else self.name
        self.name = self.name+"smeared" if self.smeared==True else self.name

        if fine == True:
            self.renormname="48fine"
        elif l==24 or l==48:
            self.renormname="24cubed"
        elif l==32 or l==64:
            self.renormname="32cubed"

        if old == True:
            fi = open("../lattice_setup/lattice_setup_" + self.name+"_old.txt",'r')
        else:
            fi = open("../lattice_setup/lattice_setup_" + self.name+".txt",'r')
        line =fi.readline()
        self.m_sea_l=[float(i) for i in line.split()]
        line =fi.readline()
        self.m_val_s=[float(i) for i in line.split()]
        line=fi.readline()
        self.m_l_phys =[float(i) for i in line.split()]
        line =fi.readline()
        self.m_s_phys=[float(i) for i in line.split()]
        #read in the lattice spacing bootstraps and set central values.
        filename = '../common_data/boot_ainv_'+str(l)+'cubed_IW_'+str(500) if fine==False else '../common_data/boot_ainv_'+str(l)+'cubedfine_IW_'+str(500)
        self.bootainv = file_io.read_file_list(filename)
        ainv_cent = {'24':1.7848,'32':2.3833,'48':1.7295,'64':2.3586,'48fine':2.774}
        self.bootainv[-1] = ainv_cent[str(l)] if fine==False else ainv_cent[str(l)+"fine"]
        self.a2 = []
        for b in range(len(self.bootainv)):
            self.a2.append(pow(self.bootainv[b],-2))




class RenormParams():
    def __init__(self,lattice,ischeme,ikin,ibasis):
        "Initialise the global variables, and those specific to this lattice"
        self.ikin = ikin
        self.ischeme = ischeme
        self.ibasis = ibasis
        self.scheme_name = ['MOM','ms'][ischeme]
        self.basis_name = ['SUSY','Lattice'][ibasis]
        self.kin_name = ['E','qq','qg','gq','gg'][ikin]
        self.T = [ [np.array([[1,0,0,0,0],[0,0,0,1,0],[0,0,0,-0.5,0.5],[0,0,1,0,0],[0,-0.5,0,0,0]])
], [np.identity(5) ] ][ibasis]

