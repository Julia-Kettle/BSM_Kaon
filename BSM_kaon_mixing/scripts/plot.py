import numpy as np
import matplotlib.pyplot as plt

def plot_interp(y_data,y_err,x_data,y_int,y_int_err,x_int,slope,name,xlabel,ylabel):
    plt.figure()
    plt.errorbar(x_data,y_data,yerr=y_err,color='r',ecolor='r',fmt='o')
    x = np.arange(0,1.05*max(x_data),(max(x_data)-x_int)/1000.0)
    y = y_int + slope*(x-x_int)
    plt.plot(x,y,'k')
    plt.errorbar(x_int,y_int,yerr=y_int_err,color='b',ecolor='b',fmt='v')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig("./plots/interpolation/"+name+".pdf",format='pdf')
    plt.close()
