#import python modules
import sys
import math
import os
#import my code
from bootstrapfit2D import BootstrapFit2D
from datasets import *
from file_io import *
from fit_functions import *
from plot_setup import *
import globaldefs
from formatting import *
from utils import *
from plot_utils import *

'''
Julia Kettle Feb 2018
main code to loop through that different renorm/basis/parameter options and for each channel carry out bootstraps of nonlinear fit
fit has form y(x1,x2) = p[0](1+ p[1]( x1 + coeff*log(x1) ) + p[3]*x2 ) with x1 = (mll/4pi*fll)^2 & x2 = a^2
fit functions defined in fit_funnctions.py

datapoints included and plotting style of point determined from an input file with form
#name  ml  ms   gr  lab col mar siz fill line
where each line is 1 datapoint

data is then read in from renormalised data stored using fetch Data
bootstraps of datasets created then for each boot non-linear fit performed. 

a plot is done and saved with plotting params defined in a user inputted file. see plot_params.txt
bootstraps saved to txt files

a template tex file is 

'''

#read the input file and return lists of datapoint properties
def fetchInput(inputfile):

    #read data from an input file of format
    # name ml ms group colours markers markersize fillstyle linestyle 
    name=[]
    ml=[]
    ms=[]      
    group=[]
    label=[]  
    colours=[]
    markers=[]
    msizes=[]
    fillstyles=[]
    linestyles=[]
    legendlab=[]

    #check at start of file
    inputfile.seek(0)

    for line in inputfile.readlines():
        if line.startswith('#'): #skip commented lines
            continue
        else:              
            data=line.split()
            name.append(data[0])
            ml.append(data[1])
            ms.append(data[2])
            group.append(int(data[3]))
            label.append(data[4])
            colours.append(data[5])
            markers.append(data[6])
            msizes.append(data[7])
            fillstyles.append(data[8])
            linestyles.append(data[9])
            legendlab.append(data[10])

    return name,ml,ms,group,label,colours,markers,msizes,fillstyles,linestyles,legendlab

#read R/B and m,a renormalised data and return in arrays of [nchan,ndp,nboots]
def fetchData(param,name,ml,ms,basis,scheme,projscheme,do2GeV=False):
   
    data=[]
    mfsq=[]
    ainv=[]
    a2=[]

    

    for i in range(len(name)):
        #set up filenames for reading the data
        fileend=name[i]+"ml_"+ml[i]+"_ms_"+ms[i]+"_"+basis+"_"+scheme+"_"+projscheme+".bin"
        afilename="../common_data/boot_ainv_"+name[i][:2]+"cubedfine_IW_500" if ("fine" in name[i]) else "../common_data/boot_ainv_"+name[i][:2]+"cubed_IW_500"
        

        
        #read in data. format is lists of numpy arrays
        if not(do2GeV):
            print("fetching data from - ../renormalised/"+param+"_"+fileend)
            print("fetch lattice spacing from - "+afilename+"\n")
            data.append(np.load(open("../renormalised/"+param+"_"+fileend,'r')))
            mfsq.append(np.load(open("../renormalised/m_4pif_sq_"+fileend,'r')))
            ainv.append(np.array(read_file_list(afilename)))
            a2.append(1.0/(ainv[i]*ainv[i]))
        else:
            print("fetching data from - ../renormalised/2GeV/"+param+"_"+fileend)
            print("fetch lattice spacing from - "+afilename+"\n")
            data.append(np.load(open("../renormalised/2GeV/"+param+"_"+fileend,'r')))
            mfsq.append(np.load(open("../renormalised/2GeV/m_4pif_sq_"+fileend,'r')))
            ainv.append(np.array(read_file_list(afilename)))
            a2.append(1.0/(ainv[i]*ainv[i]))
        
    #convert to arrays of no_points x nchan x nboots or no_points x nboots
    data=np.squeeze(np.array(data))
    #set to dimenstions [nchan,no_points,nboots]
    data=data.swapaxes(0,1)
    mfsq=np.squeeze(np.array(mfsq))
    a2=np.squeeze(np.array(a2))    
    print"\n"

    return mfsq,a2,data

#plots the fits (using params in paramfile) and saves to savefile
def plotFit(fitObj,otherAnsatzFitObj,xlabel,ylabel,paramfile,savefile):
    print "In plotFit"
    fig = plt.figure()
    ax=fig.add_subplot(111)

    print "subplot added"
    #plot central values of dataset as errorbars and fitlines
    ax=fitObj.dslist[-1].plot(ax)
    ax=fitObj.plot_fitlines(ax,0)
    print "main fit plotted"
    #plot other ansatz physical point results + fit line
    ax = otherAnsatzFitObj.dslist[-1].plotOnePoint(ax,'phys')                        
    ax = otherAnsatzFitObj.plot_fitline(ax,'phys',0)
    print "other plotted"
    
    #read figure parameters from file
    loc,bboxloc,mode,borderPad,ncol,ticklabelsize,xlabelsize,ylabelsize,subplt_adj_l,subplt_adj_b,subplt_adj_r,subplt_adj_t=readPlotParams(paramfile)

    #draw the legend
    ax=plotVerticalLine(ax,mfsqPhys)
    drawLegend(ax,loc,bboxloc,borderPad,ncol,mode)
    print "legend drawn"
    #ax=setLegend(ax,paramfile)
    #set labels and fontsize
    ax.tick_params(axis='both',labelsize = ticklabelsize)
    ax.locator_params(nticks=4)
    ax.set_ylabel(ylabel,fontsize=ylabelsize)
    ax.set_xlabel(xlabel,fontsize=xlabelsize)

    #adjust subplot dimenstions
    #fig.subplots_adjust(left=subplt_adj_l,bottom=subplt_adj_b,right=subplt_adj_r,top=subplt_adj_t)

    #save figure
    fig.savefig(savefile, format='pdf')#,bbox_inches="tight")
    print "figure saved"
    plt.close(fig)
    print "figure closed"



#main 
#loops through scheme, basis, proj etc options
# carries out fit
#plots results to figures, saves bootstraps and writes results to tex files
def main():

    if len(sys.argv) < 4:
        #check for input files and output file
        print("Usage: python "+sys.argv[0]+" <fit input file> <figure param file> <save location of plots> ")
        sys.exit(-1)
    else:
        #set file/dirs from command line
        inputfile=sys.argv[1]
        figParamfile=sys.argv[2]
        savelocation=sys.argv[3]

    try:
        do2GeV=int(float(sys.argv[4]))
    except:
        do2GeV=False

    #read input file to determine which datapoints to include in fit
    name,ml,ms,group,label,colours,markers,msizes,fillstyles,linestyles,legendlab = fetchInput(open(inputfile,'r'))
    
    #read in tex template
    with open('../results/texTEMPLATE.tex','r') as contentFile:
        texTemplate=contentFile.read()
    
    with open('../results/tablesTEMPLATE.tex','r') as contentFile:
        tabTemplate=contentFile.read()

    #check location exists and create if not
    createDir(savelocation+"/bootstraps/")

    #loop through the basis, final scheme, projection scheme & parameter & fit ansatz
    # performing a bootstrapped fit for each case
    for basis in ['Lattice','SUSY']:

        #set up texTemplate buffer - one fit texfile per basis
        texTemplate_buffer=texTemplate
        tabTemplate_buffer=tabTemplate

        for scheme in ['MOM','ms']:
            for projscheme in ['gg','qq']:
                for param in ['R','B']:
                    try:
                        mfsq,a2,data=fetchData(param,name,ml,ms,basis,scheme,projscheme,do2GeV)
                    except:
                        print "fetching data for " +scheme+" "+projscheme+" "+ basis+" failed"
                        for ansatz in ['linear','chiral']:
                            for ich in range(nchan):
                                tabEntry=param+scheme+projscheme+ansatz+str(ich+1)
                                tabEntryChi=param+scheme+projscheme+"chisq"+ansatz+str(ich+1)
                                tabEntryPlot=param+scheme+projscheme+ansatz+"PLOT"+str(ich+1)
                                texTemplate_buffer=texTemplate_buffer.replace(tabEntry,'-')
                                texTemplate_buffer=texTemplate_buffer.replace(tabEntryChi,'-')
                                texTemplate_buffer=texTemplate_buffer.replace(tabEntryPlot,'../dummy.png')
                                #paper tables
                                tabTemplate_buffer=tabTemplate_buffer.replace(tabEntry,'-')
                                tabTemplate_buffer=tabTemplate_buffer.replace(tabEntryChi,'-')
                        continue

                    #store number of channels and bootstraps(inc central)
                    nchan = len(data)
                    nboots = len(data[0,0])

                    for ich in range(nchan):
                        fitList=[]
                        datasetList=[]
                        savefileList=[]
                        for ansatz in ['linear','chiral']:
                            print "--------------------------------------------------------------------------"
                            print "--------------",scheme,projscheme,param,ich,ansatz,"-----------------"
                            print "--------------------------------------------------------------------------"

                            coeff=chiral_log_coeffs(param,ansatz,basis)

                            if ansatz == 'linear':
                                anslab='lin'
                            elif ansatz == 'chiral':
                                anslab='\chi^{PT}'

                            #create data file to save fit bootstraps
                            outDataFile=open(savelocation+"/bootstraps/"+param+str(ich+1)+"_"+basis+"_"+scheme+"_"+projscheme+"_"+ansatz+'.dat','w')
                            savefileend=param+str(ich+1)+"_"+basis+"_"+scheme+"_"+projscheme+'_'+ansatz+'.pdf'
                            savefile=savelocation+"/"+savefileend
                            
                            #create bootstrap of datasets
                            dataset=DataSet.create_bootstraps('$\mathrm{'+param+'_'+str(ich+1)+'}$',data[ich],mfsq,a2,group,label,markers,colours,msizes,fillstyles,linestyles,legendlab)
                                
                            ####################################################################################
                            #initialise and do fit, append phys result to datasets + plot
                            #fit=BootstrapFit2D(dataset,globfunc_fixedlog,[1,1,1],coeff[ich])
                            fit=BootstrapFit2D(dataset,globfunc_polynomial_m,[1,1,1,1],coeff[ich])
                            fit.dofit()
                            leglabPhys = "$\mathrm{"+param+"_"+str(ich+1)+"(0,m_\pi^{phys})^{"+anslab+"}}$"
                            if ansatz == 'linear':
                                fit.appendPhysVal(mfsqPhys,aPhys,colour='k',leglab=leglabPhys)
                            else:
                                fit.appendPhysVal(mfsqPhys,aPhys,leglab=leglabPhys)
                            fitList.append(fit)
                            datasetList.append(dataset)
                            savefileList.append(savefile)

                            #save the physical value for writing to results
                            yp=fit.bootphysvalue(mfsqPhys,aPhys)
                            yperr=np.std(yp[:-1])
                            chidof=fit.bchidof
                            fitparams=fit.bparams
                            fit.getResiduals()
                            #write bootstraps to txt file
                            outDataFile.write('# {0:6}\t{1:6}\t{2:6}\t{3:6}\t{4:6}\t{5:6}\n'.format('result','err','chidof','p[0]','p[1]','p[2]'))
                            for i in range(nboots):
                                outDataFile.write('  {0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\t{4:6f}\t{5:6f}\n'.format(yp[i],yperr,chidof[i],fitparams[i,0],fitparams[i,1],fitparams[i,2]))
                            print('  {0:6f}\t{1:6f}\t{2:6f}\t{3:6f}\t{4:6f}\t{5:6f}\n'.format(yp[i],yperr,chidof[-1],fitparams[-1,0],fitparams[-1,1],fitparams[-1,2]))
                            

                            #create target string to be replaced in texString 
                            #e.g. RMOMgglinear1 for R1_phys(err) with linear ansatz, MOM scheme, gg projection 
                            tabEntry=param+scheme+projscheme+ansatz+str(ich+1)
                            tabEntryChi=param+scheme+projscheme+"chisq"+ansatz+str(ich+1)
                            tabEntryPlot=param+scheme+projscheme+ansatz+"PLOT"+str(ich+1)
                            stringDataSet="("
                            for labeldp in fit.dslist[-1].label[:-1]:
                                stringDataSet+=labeldp+","
                            stringDataSet=stringDataSet[:-1]+")"
                            #replace strings defined above with results/plots from fits
                            texTemplate_buffer=texTemplate_buffer.replace(tabEntry,format_val_err(yp[-1],yperr))
                            texTemplate_buffer=texTemplate_buffer.replace(tabEntryChi,str(chidof[-1])[:4])
                            texTemplate_buffer=texTemplate_buffer.replace(tabEntryPlot,savefileend)
                            texTemplate_buffer=texTemplate_buffer.replace("BASIS",basis)
                            texTemplate_buffer=texTemplate_buffer.replace("DATASET",stringDataSet)

                            #table template
                            tabTemplate_buffer=tabTemplate_buffer.replace(tabEntry,format_val_err(yp[-1],yperr))
                            tabTemplate_buffer=tabTemplate_buffer.replace(tabEntryChi,str(chidof[-1])[:4])

                        #allows us to plot both fit ansatz results on plot
                        for iAnsatz in range(len(['linear','chiral'])):
                            iOther = (iAnsatz+1)%2
                            #plot current fit ansatz results
                            print "iAnsatz",iAnsatz
                            plotFit(fitList[iAnsatz],fitList[iOther],r"$m_{\pi}^2/(4\pi f_{\pi}^2)$",r""+datasetList[iAnsatz][-1].texname,figParamfile,savefileList[iAnsatz])
                            #ax = fitList[iOther].plot_fitband(ax,'phys',0)

                            
                        ######################################### END OF CHANNEL LOOP #################################################
                    
                    ############################################ END OF ANSATZ LOOP ################################################

                ############################################ END OF PARAM LOOP ##################################################

            ############################################## END OF PROJECTION LOOP ############################################
            
        ############################################## END OF SCHEME LOOP ################################################
        
        # write the corrected tex file string 
        outTexFile=open(savelocation+"/fits_"+basis+".tex",'w')
        outTexFile.write(texTemplate_buffer)

               # write the corrected tex file string 
        outTabFile=open(savelocation+"/tables_"+basis+".tex",'w')
        outTabFile.write(tabTemplate_buffer)
                                
    ############################################# END OF BASIS LOOP #################################################

############################################################## END OF MAIN ########################################################################
                            

if __name__ == "__main__":
    main()


    
