import numpy as np
from collections import OrderedDict
from globaldefs import *

def readlabelsizes(fi):
    line=fi.readline()
    tick,x,y = line.split()
    return int(tick), int(x), int(y)

def readsubplotadj(fi):
    line=fi.readline()
    l,b,r,t=line.split()
    return float(l),float(b),float(r),float(t)

def readSkipComment(fi):
    while(True):
        line=fi.readline()
        if line.startswith('#'):
            continue
        else:
            break
    return line

def readlegend(fi):
    #read legend loc
    line=readSkipComment(fi)
    loc=int(line)
    #legend bbox
    line=readSkipComment(fi)
    bboxloc=line.split()
    for i in range(len(bboxloc)):
        bboxloc[i]=float(bboxloc[i])
    bboxloc=tuple(bboxloc)
    #expand
    line=readSkipComment(fi)
    boolExpand=bool(int(line))
    mode='expand' if boolExpand else ''
    #borderpad
    line=readSkipComment(fi)
    borderPad=int(line)
    #columns
    line=readSkipComment(fi)
    ncol=int(line)

    return loc,bboxloc,mode,borderPad,ncol


def readPlotParams(infilename):
    fi=open(infilename,'r')

    eof = False
    while not(eof):
        line = fi.readline()
        if line.startswith('#End'):
            eof==True
            break
        elif line.startswith('#legend'):
            loc,bboxloc,mode,borderPad,ncol=readlegend(fi)
        elif line.startswith('#labelsizes'):
            ticksize,xsize,ysize=readlabelsizes(fi)
        elif line.startswith('#subplotadjust'):
            ladj,badj,radj,tadj=readsubplotadj(fi)
    return loc,bboxloc,mode,borderPad,ncol,ticksize,xsize,ysize,ladj,badj,radj,tadj


def plotVerticalLine(ax,x,label=''):
    ymin,ymax=ax.get_ylim()
    ax.plot(mfsqPhys*np.ones([1000]),np.linspace(ymin,ymax,1000),color='k',linestyle=':',label=r'$(\mathrm{m_\pi/4 \pi f_\pi)_{phys}^2}$')
    ax.set_ylim(ymin,ymax)
    return ax

def drawLegend(ax,loc,bboxloc,borderPad,ncol,mode):
    
    handles, labels = ax.get_legend_handles_labels()
    # remove the errorbars
    #for h in handles:
        #print h
    #handles = [h[0] for h in handles]
    by_label = OrderedDict(zip(labels, handles))
    #adjust axis location to fit in legend
    chartBox = ax.get_position()
    ax.set_position([chartBox.x0, chartBox.y0*1.15, chartBox.width, chartBox.height*0.85])
    ax.legend(by_label.values(), by_label.keys(),loc=loc,bbox_to_anchor=bboxloc, ncol=ncol, mode=mode, borderaxespad=borderPad,prop={'size':14},fontsize=12,numpoints=1)

    return ax

