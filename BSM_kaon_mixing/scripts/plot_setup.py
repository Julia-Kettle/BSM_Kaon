import numpy as np
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
    loc=line
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


def plotVerticalLine(ax,x):
    ymin,ymax=ax.get_ylim()
    ax.plot(mfsqPhys*np.ones([1000]),np.linspace(ymin,ymax,1000),color='k')
    ax.set_ylim(ymin,ymax)
    return ax