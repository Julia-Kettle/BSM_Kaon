def Global_Params(ibeta)
    latticedim = ['24','32','48','64']
    volname = ['24cubed','32cubed','48cubed','64cubed']
    name_kin = ['E','gg','qq','qg','gq']
    name_scheme = ['MOM','ms']
    name_basis = ['Lattice', 'SUSY']
    return latticedim[ibeta],volname[ibeta],name_kin[ibeta],name_scheme[ibeta],name_basis[ibeta]
