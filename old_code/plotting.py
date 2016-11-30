import matplotlib.pyplot as plt
import numpy as np
from math import ceil

def plot_Global(p,function,y1,y2,yerr1, yerr2,yperr,m1,m2,mp,a1,a2,xlabel,ylabel,name,y3=[],yerr3=[],m3=[],a3=[]):
    '''
    call this from boot strap global
    plot R or B against (m/f)^2
    plot function fitted
    plot physical mass + continuum
    '''
    plt.figure()
    plt.errorbar(m1[:,-1],y1[:,-1],yerr=yerr1,ecolor='r',fmt='v',color='r',markersize=5,label='24cubed48')
    plt.errorbar(m2[:,-1],y2[:,-1],yerr=yerr2,ecolor='b',fmt='v',color='b',markersize=5,label='32cubed64')
    if type(y3)==np.ndarray and type(m3)==np.ndarray and type(yerr3)==np.ndarray:
        phys_lab=['48cubed64 (phys pt)','64cubed128 (phys pt)']
        phys_col=['g','y']
        for i in range(len(y3)):
            plt.errorbar(m3[i,-1],y3[i,-1],yerr=yerr3[i],ecolor=phys_col[i],fmt='v',color=phys_col[i],markersize=5,label=phys_lab[i])
   #Plotted all the data points with errorbars

    xmax=plt.xlim()[1]
    xmin=plt.xlim()[0]
    #save limits 

    x = range(0,int(ceil(xmax*1.1*100000)))
    #set range over which function plotted.
    f1=[]
    f2=[]
    f3=[]
    f4=[]
    #for i in range(len(x)):
    #    x[i] = x[i]/100000.0
    #    f1.append(function(p,x[i],1.0/(a1[-1]*a1[-1])))
    #    f2.append(function(p,x[i],1.0/(a2[-1]*a2[-1])))
    #    f3.append(function(p,x[i],1.0/(a3[0][-1]*a3[0][-1])))
    #    f4.append(function(p,x[i],1.0/(a3[1][-1]*a3[1][-1])))
#find the continuum results
    #yp = function(p,mp,0)
    #p1 = function(p,mp,1.0/(a1[-1]*a1[-1]))
    #p2 = function(p,mp,1.0/(a2[-1]*a2[-1]))

    #plot the fit lines
    #plt.plot(x,f1,'r-')
    #plt.plot(x,f2,'b-') 
    #plt.plot(x,f3,'g-')
    #plt.plot(x,f4,'y-')
    #plt.errorbar(mp,yp,yerr=yperr,fmt='mo',ecolor='m',markersize=6.5,label='continuum limit phys pt')
    #plt.plot(mp,p1,'rx',markersize=5)
    #plt.plot(mp,p2,'bx',markersize=5)


    ylims=plt.ylim()
    plt.vlines(mp,ylims[0],ylims[1],colors='k', linestyles='dashed',label='physical pion mass')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
        ncol=3, mode="expand", borderaxespad=0.,prop={'size':7})
    plt.xlabel(xlabel,fontsize=16)
    plt.ylabel(ylabel,fontsize=20)
    plt.ylim(ylims)
    plt.xlim([0,1.05*xmax])
    #plt.show()
    plt.savefig(name+'.eps',format='eps')
    plt.close()
    
