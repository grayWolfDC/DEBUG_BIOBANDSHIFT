# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:09:46 2015

@author: Erdem Karako:ylu:
This is a script for plotting comparisons between bandshifted data 
and data originally acquired (via another sensor) at the target band
with time/space overlap.
"""
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

def GetCompData(fName,dataLabel='comparison',lolim = -150,hilim = 150):
    ncData = nc.Dataset(fName)
    for key in ncData.variables:
        if key.startswith(dataLabel):
            comp = ncData.variables[key]
    compAr = np.array(comp)
    flatCompAr = compAr.flatten()
    flatCompAr = flatCompAr[np.isfinite(flatCompAr)]
    ncData.close()
    lolimIdx = flatCompAr >= lolim
    hilimIdx = flatCompAr <= hilim
    goodIdx = lolimIdx * hilimIdx
    return flatCompAr[goodIdx]

def PlotComparisons(objList,colorList,plotTitle='Band-shifting comparison',
                    binRange=(-100,100)):
    bns=np.linspace(binRange[0],binRange[1],100);
    plt.figure(figsize=(14,6))
    plt.title(plotTitle,fontsize=18)
    percLo = []
    percHi = []
    percMed = []
    for i,obj in enumerate(objList):
        obj.n,obj.bins,_ = plt.hist(obj.data,bins=bns,normed='yes',
                                    histtype='stepfilled',alpha=0.45,
                                    label=obj.tag,color=colorList[i])

        lo = np.percentile(obj.data,2.5)
        hi = np.percentile(obj.data,97.5)
        med = np.percentile(obj.data,50)
        percLo.append(lo)
        percHi.append(hi)
        percMed.append(med)
    axes = plt.gca()
    ymax= max(axes.get_ylim())
    plt.ylabel('freq.',fontsize=16)
    plt.xlabel('% difference',fontsize=16)
    plt.legend(loc='best',fontsize=16)
    plt.vlines(0,0,ymax,linestyle='--',lw=3);
    #plt.vlines(-10,-10,ymax,linestyle=':',lw=2);
    #plt.vlines(10,10,ymax,linestyle=':',lw=2);
    plt.ylim(ymax=ymax)
    #plt.plot(med1,0.004,'ko')
    #plt.hlines(0.004,lo1,hi1,linestyle='--',lw=2);
    extraXticks=[-10,10]
    plt.xticks(list(plt.xticks()[0])+extraXticks,fontsize=14)
    plt.yticks(fontsize=14);
    
    
class DataCompObj(object):
    def __init__(self,fName,tag):
        self.fName = fName
        self.tag = tag
        self.data = GetCompData(self.fName)
        
        