# histogram_script.py:
# Constructs the histograms for the nebular parameters 
# (c(Hb), Te, Ne, ionic/elemental abundances and ICFs, saved in the output plots folder)
# (C) Stavros Akras


from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import array
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
from astropy.io import fits
from astropy.io.fits import getdata
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
from scipy.interpolate import spline
import matplotlib.patches as patches
import matplotlib.font_manager as font_manager
import random
import math
import numpy as np

from scipy.stats import norm
import seaborn as sns
import sys
from satellite import mean_value_script as mvs

### values and errors as input

def histogram_plots(data,data_err,name,line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered):
    index_size=param_mod_name.index("total_num_pixels_horiz")
    sizex=param_model_values[index_size]
    sizey=param_model_values[index_size]
    
    file10 = open('general_output_file.txt','a')

    
    newdata=[]
    newerror=[]
    k=0
    data2=np.nan_to_num(data)
    data_err2=np.nan_to_num(data_err)
    for j in range(0,sizey):
        for i in range(0,sizex):
            if data2[i,j]!=0:
                newdata.append(0)
                newerror.append(0)
                newdata[k]=data2[i,j]
                newerror[k]=data_err2[i,j]
                k=k+1
    
    newdata2=[]
    kk=0
    if len(newdata)!=0: 
        meanvalue=mvs.meanvalue(newdata,0)
        for i in range(0,k):
            if newdata[i]<meanvalue*3. and newdata[i]>meanvalue/3. :
                newdata2.append(0)
                newdata2[kk]=newdata[i]
                kk=kk+1
    else: 
        meanvalue=0.
        kk=1
        newdata2.append(0)
            
    fig,ax1=plt.subplots(figsize=(11.5,7.2))
    datamu, datastd = norm.fit(newdata2)
    if kk!=0: binSR=int(round(np.log10(kk)*3.3+1.))
    else: binSR=1
    plt.hist(newdata2, density=False, bins=binSR)
    plt.ylabel("Counts",fontsize=21,labelpad=6)
    ax1.tick_params(axis='x', labelsize=21, pad=6)
    ax1.tick_params(axis='y', labelsize=21, pad=6)
    ax1.ticklabel_format(axis='x',style='sci',scilimits=(1,4),labelsize=22)
    xmin, xmax = plt.xlim()
    
    
    index_name=param_estimated.index(name)
    if (param_requered[index_name]=='yes'):
        if name=="c_Hb": plt.xlabel(r'c(H$\beta$)',fontsize=21,labelpad=6)
        if name=="Te(NII6548_84)_Ne(SII6716_31)": plt.xlabel(r'Te([N II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(NII6548_84)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OI6300_63)_Ne(SII6716_31)": plt.xlabel(r'Te([O I]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(OI6300_63)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OII3727_29_7320_30)_Ne(SII6716_31)": plt.xlabel(r'Te([O II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(OII3727_29_7320_30)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OIII4959_5007)_Ne(SII6716_31)": plt.xlabel(r'Te([O III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(OIII4959_5007)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(SIII6312_9069)_Ne(SII6716_31)": plt.xlabel(r'Te([S III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(SIII6312_9069)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OII3727_29_7320_30)_Ne(OII3727_29)": plt.xlabel(r'Te([O II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(OII3727_29)_Te(OII3727_29_7320_30)": plt.xlabel(r'Ne([O II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(NII6548_84)_Ne(OII3727_29)": plt.xlabel(r'Te([N II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(OII3727_29)_Te(NII6548_84)": plt.xlabel(r'Ne([O II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OI6300_63)_Ne(OII3727_29)": plt.xlabel(r'Te([O I]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(OII3727_29)_Te(OI6300_63)":  plt.xlabel(r'Ne([O II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OIII4959_5007)_Ne(ClIII5517_38)": plt.xlabel(r'Te([O III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ClIII5517_38)_Te(OIII4959_5007)": plt.xlabel(r'Ne([Cl III]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(SIII6312_9069)_Ne(ClIII5517_38)": plt.xlabel(r'Te([S III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ClIII5517_38)_Te(SIII6312_9069)": plt.xlabel(r'Ne([Cl III]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OIII4959_5007)_Ne(ArVI4712_40)": plt.xlabel(r'Te([O III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ArVI4712_40)_Te(OIII4959_5007)": plt.xlabel(r'Ne([Ar IV]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(SIII6312_9069)_Ne(ArVI4712_40)": plt.xlabel(r'Te([S III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ArVI4712_40)_Te(SIII6312_9069)": plt.xlabel(r'Ne([Ar IV]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(NII6548_84)_Ne(ClIII5517_38)": plt.xlabel(r'Te([N II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ClIII5517_38)_Te(NII6548_84)": plt.xlabel(r'Ne([Cl III]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="He_adund_KB": plt.xlabel(r'He abundance (KB)',fontsize=21,labelpad=6)
        if name=="He_ICF_KB": plt.xlabel(r'He ICF (KB)',fontsize=21,labelpad=6)
        if name=="He_abund_DIMS": plt.xlabel(r'He abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="He_ICF_DIMS": plt.xlabel(r'He ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="N_adund_KB": plt.xlabel(r'N abundance (KB)',fontsize=21,labelpad=6)
        if name=="N_ICF_KB": plt.xlabel(r'N ICF (KB)',fontsize=21,labelpad=6)
        if name=="N_abund_DIMS": plt.xlabel(r'N abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="N_ICF_DIMS": plt.xlabel(r'N ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="O_adund_KB": plt.xlabel(r'O abundance (KB)',fontsize=21,labelpad=6)
        if name=="O_ICF_KB": plt.xlabel(r'O ICF (KB)',fontsize=21,labelpad=6)
        if name=="O_abund_DIMS":  plt.xlabel(r'O abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="O_ICF_DIMS": plt.xlabel(r'O ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="S_adund_KB": plt.xlabel(r'S abundance (KB)',fontsize=21,labelpad=6)
        if name=="S_ICF_KB": plt.xlabel(r'S ICF (KB)',fontsize=21,labelpad=6)
        if name=="S_abund_DIMS": plt.xlabel(r'S abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="S_ICF_DIMS": plt.xlabel(r'S ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="Ne_adund_KB": plt.xlabel(r'Ne abundance (KB)',fontsize=21,labelpad=6)
        if name=="Ne_ICF_KB": plt.xlabel(r'Ne ICF (KB)',fontsize=21,labelpad=6)
        if name=="Ne_abund_DIMS": plt.xlabel(r'Ne abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="Ne_ICF_DIMS": plt.xlabel(r'Ne ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="Ar_adund_KB": plt.xlabel(r'Ar abundance (KB)',fontsize=21,labelpad=6)
        if name=="Ar_ICF_KB": plt.xlabel(r'Ar ICF (KB)',fontsize=21,labelpad=6)
        if name=="Ar_abund_DIMS": plt.xlabel(r'Ar abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="Ar_ICF_DIMS": plt.xlabel(r'Ar ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="Cl_adund_KB": plt.xlabel(r'Cl abundance (KB)',fontsize=21,labelpad=6)
        if name=="Cl_ICF_KB": plt.xlabel(r'Cl ICF (KB)',fontsize=21,labelpad=6)
        if name=="Cl_abund_DIMS":  plt.xlabel(r'Cl abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="Cl_ICF_DIMS": plt.xlabel(r'Cl ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="NO_ratio_KB": plt.xlabel(r'NO abundance ratio (KB)',fontsize=21,labelpad=6) 
        if name=="NO_ratio_DIMS": plt.xlabel(r'NO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        if name=="SO_ratio_KB": plt.xlabel(r'SO abundance ratio (KB)',fontsize=21,labelpad=6)  
        if name=="SO_ratio_DIMS": plt.xlabel(r'SO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        if name=="ArO_ratio_KB": plt.xlabel(r'ArO abundance ratio (KB)',fontsize=21,labelpad=6)  
        if name=="ArO_ratio_DIMS":plt.xlabel(r'ArO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        if name=="NeO_ratio_KB": plt.xlabel(r'NeO abundance ratio (KB)',fontsize=21,labelpad=6)  
        if name=="NeO_ratio_DIMS": plt.xlabel(r'NeO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        if name=="ClO_ratio_KB": plt.xlabel(r'ClO abundance ratio (KB)',fontsize=21,labelpad=6) 
        if name=="ClO_ratio_DIMS": plt.xlabel(r'ClO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        
        plt.savefig("output_plots/hist"+name+".png")
        plt.savefig("output_plots/hist"+name+".pdf")
#    plt.show()
        print >>file10, "#######################################################"
        print >>file10, "#######################################################"
        print >>file10, name,':(mean value)',datamu, "+/- (std value)", datastd

    file10.close()
    
    return

### values as input 

def histogram_plots2(data,name,line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered):
    index_size=param_mod_name.index("total_num_pixels_horiz")
    sizex=param_model_values[index_size]
    sizey=param_model_values[index_size]
    
    file10 = open('general_output_file.txt','a')

    newdata=[]
    newerror=[]
    k=0
    data2=np.nan_to_num(data)
    data_err2=np.nan_to_num(data_err)
    for j in range(0,sizey):
        for i in range(0,sizex):
            if data2[i,j]!=0:
                newdata.append(0)
                newerror.append(0)
                newdata[k]=data2[i,j]
                k=k+1
     
    newdata2=[]
    kk=0
    if len(newdata)!=0: 
        meanvalue=mvs.meanvalue(newdata,0)
        for i in range(0,k):
            if newdata[i]<meanvalue*3 and newdata[i]>meanvalue/3 :
                newdata2.append(0)
                newdata2[kk]=newdata[i]
                kk=kk+1
    else: 
        meanvalue=0.
        kk=1
        newdata2.append(0)
    
    fig,ax1=plt.subplots(figsize=(11.5,7.2))
    datamu, datastd = norm.fit(newdata2) 
    if kk!=0: binSR=int(round(np.log10(kk)*3.3+1.))
    else: binSR=1
    plt.hist(newdata2, density=False, bins=binSR)
    plt.ylabel("Counts",fontsize=21,labelpad=6)
    ax1.tick_params(axis='x', labelsize=21, pad=6)
    ax1.tick_params(axis='y', labelsize=21, pad=6)
    ax1.ticklabel_format(axis='x',style='sci',scilimits=(1,4),labelsize=22)
    plt.xlim(xmin,xmax)
    
    index_name=param_estimated.index(name)
    if (param_requered[index_name]=='yes'):
        if name=="c_Hb": plt.xlabel(r'c(H$\beta$)',fontsize=21,labelpad=6)
        if name=="Te(NII6548_84)_Ne(SII6716_31)": plt.xlabel(r'Te([N II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(NII6548_84)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OI6300_63)_Ne(SII6716_31)": plt.xlabel(r'Te([O I]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(OI6300_63)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OII3727_29_7320_30)_Ne(SII6716_31)": plt.xlabel(r'Te([O II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(OII3727_29_7320_30)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OIII4959_5007)_Ne(SII6716_31)": plt.xlabel(r'Te([O III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(OIII4959_5007)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(SIII6312_9069)_Ne(SII6716_31)": plt.xlabel(r'Te([S III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(SII6716_31)_Te(SIII6312_9069)": plt.xlabel(r'Ne([S II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OII3727_29_7320_30)_Ne(OII3727_29)": plt.xlabel(r'Te([O II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(OII3727_29)_Te(OII3727_29_7320_30)": plt.xlabel(r'Ne([O II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(NII6548_84)_Ne(OII3727_29)": plt.xlabel(r'Te([N II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(OII3727_29)_Te(NII6548_84)": plt.xlabel(r'Ne([O II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OI6300_63)_Ne(OII3727_29)": plt.xlabel(r'Te([O I]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(OII3727_29)_Te(OI6300_63)":  plt.xlabel(r'Ne([O II]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OIII4959_5007)_Ne(ClIII5517_38)": plt.xlabel(r'Te([O III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ClIII5517_38)_Te(OIII4959_5007)": plt.xlabel(r'Ne([Cl III]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(SIII6312_9069)_Ne(ClIII5517_38)": plt.xlabel(r'Te([S III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ClIII5517_38)_Te(SIII6312_9069)": plt.xlabel(r'Ne([Cl III]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(OIII4959_5007)_Ne(ArVI4712_40)": plt.xlabel(r'Te([O III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ArVI4712_40)_Te(OIII4959_5007)": plt.xlabel(r'Ne([Ar IV]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(SIII6312_9069)_Ne(ArVI4712_40)": plt.xlabel(r'Te([S III]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ArVI4712_40)_Te(SIII6312_9069)": plt.xlabel(r'Ne([Ar IV]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="Te(NII6548_84)_Ne(ClIII5517_38)": plt.xlabel(r'Te([N II]) (K)',fontsize=21,labelpad=6)
        if name=="Ne(ClIII5517_38)_Te(NII6548_84)": plt.xlabel(r'Ne([Cl III]) (cm$^{-3}$)',fontsize=21,labelpad=6)
        if name=="He_adund_KB": plt.xlabel(r'He abundance (KB)',fontsize=21,labelpad=6)
        if name=="He_ICF_KB": plt.xlabel(r'He ICF (KB)',fontsize=21,labelpad=6)
        if name=="He_abund_DIMS": plt.xlabel(r'He abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="He_ICF_DIMS": plt.xlabel(r'He ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="N_adund_KB": plt.xlabel(r'N abundance (KB)',fontsize=21,labelpad=6)
        if name=="N_ICF_KB": plt.xlabel(r'N ICF (KB)',fontsize=21,labelpad=6)
        if name=="N_abund_DIMS": plt.xlabel(r'N abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="N_ICF_DIMS": plt.xlabel(r'N ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="O_adund_KB": plt.xlabel(r'O abundance (KB)',fontsize=21,labelpad=6)
        if name=="O_ICF_KB": plt.xlabel(r'O ICF (KB)',fontsize=21,labelpad=6)
        if name=="O_abund_DIMS":  plt.xlabel(r'O abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="O_ICF_DIMS": plt.xlabel(r'O ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="S_adund_KB": plt.xlabel(r'S abundance (KB)',fontsize=21,labelpad=6)
        if name=="S_ICF_KB": plt.xlabel(r'S ICF (KB)',fontsize=21,labelpad=6)
        if name=="S_abund_DIMS": plt.xlabel(r'S abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="S_ICF_DIMS": plt.xlabel(r'S ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="Ne_adund_KB": plt.xlabel(r'Ne abundance (KB)',fontsize=21,labelpad=6)
        if name=="Ne_ICF_KB": plt.xlabel(r'Ne ICF (KB)',fontsize=21,labelpad=6)
        if name=="Ne_abund_DIMS": plt.xlabel(r'Ne abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="Ne_ICF_DIMS": plt.xlabel(r'Ne ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="Ar_adund_KB": plt.xlabel(r'Ar abundance (KB)',fontsize=21,labelpad=6)
        if name=="Ar_ICF_KB": plt.xlabel(r'Ar ICF (KB)',fontsize=21,labelpad=6)
        if name=="Ar_abund_DIMS": plt.xlabel(r'Ar abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="Ar_ICF_DIMS": plt.xlabel(r'Ar ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="Cl_adund_KB": plt.xlabel(r'Cl abundance (KB)',fontsize=21,labelpad=6)
        if name=="Cl_ICF_KB": plt.xlabel(r'Cl ICF (KB)',fontsize=21,labelpad=6)
        if name=="Cl_abund_DIMS":  plt.xlabel(r'Cl abundance (DIMS)',fontsize=21,labelpad=6)
        if name=="Cl_ICF_DIMS": plt.xlabel(r'Cl ICF (DIMS)',fontsize=21,labelpad=6) 
        if name=="NO_ratio_KB": plt.xlabel(r'NO abundance ratio (KB)',fontsize=21,labelpad=6) 
        if name=="NO_ratio_DIMS": plt.xlabel(r'NO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        if name=="SO_ratio_KB": plt.xlabel(r'SO abundance ratio (KB)',fontsize=21,labelpad=6)  
        if name=="SO_ratio_DIMS": plt.xlabel(r'SO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        if name=="ArO_ratio_KB": plt.xlabel(r'ArO abundance ratio (KB)',fontsize=21,labelpad=6)  
        if name=="ArO_ratio_DIMS":plt.xlabel(r'ArO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        if name=="NeO_ratio_KB": plt.xlabel(r'NeO abundance ratio (KB)',fontsize=21,labelpad=6)  
        if name=="NeO_ratio_DIMS": plt.xlabel(r'NeO abundance ratio (DIMS)',fontsize=21,labelpad=6)
        if name=="ClO_ratio_KB": plt.xlabel(r'ClO abundance ratio (KB)',fontsize=21,labelpad=6) 
        if name=="ClO_ratio_DIMS": plt.xlabel(r'ClO abundance ratio (DIMS)',fontsize=21,labelpad=6)


        plt.savefig("output_plots/hist"+name+".png")
        plt.savefig("output_plots/hist"+name+".pdf")
#    plt.show()
        print >>file10, "#######################################################"
        print >>file10, "#######################################################"
        print >>file10, name,':(mean value)',datamu, "+/- (std value)", datastd

    file10.close()
    
    return
