class Renormalisation:



    def Renormalise(self,doR,doB):
        #Pass the filename with the Zfactors, basis, kinematics and vol to renormalise the ratio and the bag
        self.ZA_ZP()
        Zfilename = "data/Z" + self.scheme_name + "_boot_mu_match_" + self.lattice_name_Zs + "cubed_" + self.kin_name + "_block.out"
        self.Read_Z(Zfilename)
        self.N=np.squeeze(np.dot(self.T,self.N_ren))
        if doB == 1:
            self.MultBZ()
        if doR ==1:
            self.MultRZ()


    def MultRZ(self,ratio,z_bag,basis):
        """
        Renormalise the Ratios by multiplying with the Z_bag
        bootstrap by bootstrap
        """
        #set up the renormalised ratio arrays
        dimr = np.shape(ratio)
        ren_ratio = np.zeros(dimratio)
        ren_phys_ratio = np.zeros(dimratio)
        #swap last bootstrap axis with channel axis to allow matrix multiplication, bootstrap by bootstrap
        ren_ratio = np.swapaxes(self.ren_ratio,0,-1)
        ren_phys_ratio = np.swapaxes(self.ren_phys_ratio,0,-1)
        ratio = np.swapaxes(self.ratio,0,-1)

        #loop through ratio, bootstraps and mass indices
        #multiply ratios by Z factor matrix and the basis matrix 
        for index in np.ndindex(np.shape(ratio)[:-1]):
            ren_ratio[index] = np.dot(ratio[index],Z_bag)/self.bZ_scheme[0][0]
            ren_phys_ratio[index] = np.dot(Basis,ren_ratio[index])
        ren_ratio = np.swapaxes(ren_ratio,0,-1)
        ren_phys_ratio = np.swapaxes(ren_phys_ratio,0,-1)
        ratio = np.swapaxes(ratio,0,-1)
        #print "Renormalised R is- \n",self.ren_phys_R[:,0,0,-1]
        return ren_phys_ratio

    def MultBZ(self,bag,z_bag,zp,za,basis,normalisation):
        """
        Renormalise the Bag parameters 
        """
        #set up arrays for the renormalised bags
        #swap axes so that channel axis last
        dimbag = np.shape(bag)
        bag = np.swapaxes(bag,0,-1)
        ren_bag = np.swapaxes(np.zeros(dimbag),0,-1)
        ren_phys_bag = np.swapaxes(np.zeros(dimbag),0,-1)
        bag_phys = np.swapaxes(np.zeros(dimbag),0,-1)
        #loop through masses and boots then carry out matrix multiplication with renorm factors
        for index in np.ndindex(np.shape(self.bag)[:-1]):
            self.ren_bag[index] = (za*za*(np.dot(z_bag[:,:],bag[index]))
            self.ren_phys_bag[index] = np.dot(basis,ren_bag[index])
            for ich in range(1,5):
                bag_phys[index][ich] = ren_phys_bag[index][ich]/(zp*zp*normalisation[ich])
            bag_phys[index][0] = ren_phys_bag[index][0]/(za*za*normalisation[0])
        bag = np.swapaxes(bag,0,-1)
        ren_bag = np.swapaxes(ren_bag,0,-1)
        ren_phys_bag = np.swapaxes(ren_phys_bag,0,-1)
        bag_phys = np.swapaxes(bag_phys,0,-1)
        #print "Renormalised B is- \n",self.B_phys[:,0,0,-1]
        return bag_phys
