import sys
import random

def bootstrapGen(central,error,nboot,filename):
    fout = open(filename,'w')

    for i in range(nboot):
        b = random.gauss(central,error)
        fout.write(str(b))
        fout.write("\n")
    fout.write(str(central))

if __name__=="__main__":
    ainv = [1.730,2.359]
    aerr = [0.004,0.007]
    filename = ["boot_ainv_48cubed_IW_500","boot_ainv_64cubed_IW_500"]
    for i in range(2):
        bootstrapGen(ainv[i],aerr[i],500,filename[i])
