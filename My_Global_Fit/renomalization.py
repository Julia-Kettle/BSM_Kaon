import scipy.numpy as np
import readfiles as rf

def Read_Zs():
	rf.read_datacol_file(zfilename)
	for iop in range(nop):
		for jop in range(nop):
			for iboot in range(nboot):
				i=i+1
				Z[iop,jop,iboot]=data[i]
	
