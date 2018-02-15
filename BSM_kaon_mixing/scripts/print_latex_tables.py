import numpy as np
import formatting as fm


def print_latex_table(val,name,yphyslin,eyphyslin,chidoflin,echidoflin,yphyschir,eyphyschir,chidofchir,echidofchir):
    filename="../latex_tables/"+name+".tex"
    fo=open(filename,'w')

    fo.write("\\begin{table} \n")
    fo.write("\\begin{tabular}{ c |  c c | c c  | c } \n") 
    fo.write("\hline \hline \n")
    fo.write("\& linear fit & $\chi^2$ /dof & chiral PT fit & $\chi^2$ /dof & RBC/UKQCD16 \\\ \n")
    fo.write("\hline \n")
    for ich in range(np.size(yphyslin)):
        fo.write("$"+str(val)+"_"+str(ich+1)+"$ & ")
        fo.write(fm.format_val_err(yphyslin[ich],eyphyslin[ich]) + " & " + str(round(chidoflin[ich],2)) +" & ")
        fo.write(fm.format_val_err(yphyschir[ich],eyphyschir[ich]) + " & " + str(round(chidofchir[ich],2)) + " & ")
        fo.write("- \\\ \n") #this needs to be replaced with values from Nicolas' paper
    fo.write("\hline \n")
    fo.write("\end{tabular} \n")
    fo.write("\end{table} \n")
######################################################################################################################################################

