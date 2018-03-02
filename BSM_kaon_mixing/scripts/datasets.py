import numpy as np 
from matplotlib import pyplot as plt 
import globaldefs

class DataPoint(object):
    #object containing plotting params and values for data type y=f(x1,x2)
    def __init__(self,group,label,x1,x2,y,ey=0,):
        #initialise data point - contains values y,x1,x2 and plotting params
        self.y=y
        self.x1=x1
        self.x2=x2
        self.ey=ey
        self.group=group
        self.label=label
        #set default colour and marker
        self.marker='o'
        self.colour='m'
        self.msize='12'
        self.fillstyle='none'
        self.linestyle='-'
        #self.no_points=len(y)

class DataSet(object):
    #object containing lists of plotting params and values for a number of data points of type y=f(x1,x2)
    #contains functionality to create, append, create a set of bootstraps of Datasets and return the minimums and maximums
    
    def __init__(self,name,y,yerr,x1,x2,group,label,marker,colour,msize,fillstyle,linestyle):
        #initialise data set - contains lists of datapoint values y,x1,x2 and plotting params
        #label of dataset in latex formatting
        self.texname=name
        #number of points
        self.no_points=len(y)
        # values of points
        self.y=y
        self.x1=x1
        self.x2=x2
        self.ey=yerr

        #allows to set "group" of each point (i.e. if 2 points from same lattice spacing cane label)
        self.group=group
        #label each datapoint - for legend on plot
        self.label=label

        #plotting params
        self.marker=marker
        self.colour=colour
        self.msize=msize
        self.fillstyle=fillstyle
        self.linestyle=linestyle

    def append(self,point):
        #append a to list of data points from data point type
        self.y=np.append(self.y,point.y)
        self.ey=np.append(self.ey,point.ey)
        self.x1=np.append(self.x1,point.x1)
        self.x2=np.append(self.x2,point.x2)
        self.group=np.append(self.group,point.group)
        self.label=np.append(self.label,point.label)
        self.colour=np.append(self.colour,point.colour)
        self.marker=np.append(self.marker,point.marker)
        self.msize=np.append(self.msize,point.msize)
        self.fillstyle=np.append(self.fillstyle,point.fillstyle)
        self.linestyle=np.append(self.linestyle,point.linestyle)
        self.no_points+=1

    def __str__(self):
        #define print statement for a dataset
        #prints table of labels and values
        string = "label       x1           x2           y          err\n"
        for i in range(self.no_points):
            string+=str(self.label[i]) + "      " + str(round(self.x1[i],5))+ "      " +str(round(self.x2[i],5)) + "      " + str(round(self.y[i],5))+ "      " + str(round(self.ey[i],5))+"\n"
        return string
    
    def plot(self,ax,axis=0):
            #plot errorbar of dataset to axis
            x=self.x2 if (axis==1) else self.x1
            for i in range(self.no_points):
                ax.errorbar(x[i],self.y[i],yerr=self.ey[i],ecolor=self.colour[i],color=self.colour[i],
                fmt=self.marker[i],fillstyle=self.fillstyle[i],markersize=self.msize[i],label=self.label[i])
            return ax

    def min(self,axis):
            #return minimum along chosen x1,x2 or y
        if axis == 0:
            return min(self.x1)
        elif axis == 1:
            return min(self.x2)
        elif axis == 2:
            return min(self.y)
        else:
            print("Error: axis must be 0,1 or 2 to represent x1,x2 or y")
            return 0

    def max(self,axis):
            #return maximum along chosen x1,x2 or y
        if axis == 0:
            return max(self.x1)
        elif axis == 1:
            return max(self.x2)
        elif axis == 2:
            return max(self.y)
        else:
            print("Error: axis must be 0,1 or 2 to represent x1,x2 or y")
            return 0
    
    @staticmethod 
    def create_from_datapoints(self,datapoints,name):
        #create dataset from list of datapoints
        y=[]
        ey=[]
        x1=[]
        x2=[]
        group=[]
        label=[]
        marker=[]
        msize=[]
        fillstyle=[]
        for point in range(len(datapoints)):
            y.append(point.y)
            ey.append(point.ey)
            x1.append(point.x1)
            x2.append(point.x2)
            group.append(point.group)
            label.append(point.label)
            marker.append(point.marker)
            msize.append(point.msize)
            fillstyle.append(point.fillstyle)
        return DataSet(name,y,x1,x2,group,label,marker,colour,msize,fillstyle)
   
    @staticmethod
    def create_bootstraps(name,y,x1,x2,group,label,marker,colour,msize,fillstyle,linestyle):
        #create n datasets from no_points by n arrays of y,x1,x2 and no_points lists of plot params
        nboots=(np.size(y[0]))
        bdatasets=[]
        ey=np.std(y[:,:-1],axis=1)
        for i in range(nboots):
            bdatasets.append(DataSet(name,y[:,i],ey,x1[:,i],x2[:,i],group,label,marker,colour,msize,fillstyle,linestyle))
        return bdatasets