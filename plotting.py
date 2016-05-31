import matplotlib.pyplot as plt
import numpy as np

def plot_Global(p,function,y1,y2,yerr1, yerr2,m1,m2,mp,a1,a2,xlabel,ylabel,name,y3=[],yerr3=[],m3=[]):
    '''
    call this from boot strap global
    plot R or B against (m/f)^2
    plot function fitted
    plot physical mass + continuum
    '''

    x = range(0,10000)
    f1=[]
    f2=[]
    for i in range(len(x)):
        x[i] = x[i]/100000.0
        f1.append(function(p,x[i],a1))
        f2.append(function(p,x[i],a2))
        #print x[i], f1[i], f2[i]
    yp = function(p,mp,0)
    p1 = function(p,mp,a1)
    p2 = function(p,mp,a2)
    #print x
    #print f1
    plt.figure()
    #print len(m1[0:-1]), len(y1), len(yerr1)
    plt.errorbar(m1[:,-1],y1[:,-1],yerr=yerr1)
    plt.errorbar(m2[:,-1],y2[:,-1],yerr=yerr2)
    plt.plot(x,f1,'r--')
    plt.plot(x,f2,'b--')
    plt.plot(mp,yp,'kx',)
    #plt.plot(mp,p1,'kx',)
    #plt.plot(mp,p2,'kx',)

    if y3 and m3 and yerr3:
        for i in range(len(y3)):
            plt.errorbar(m3[i],y3[i],'go',yerr=yerr3[i])

    plt.xlabel(xlabel,fontsize=16)
    plt.ylabel(ylabel,fontsize=20)
    #plt.show()
    plt.savefig(name,format='pdf')

    return yp
