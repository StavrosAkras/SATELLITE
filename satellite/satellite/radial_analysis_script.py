# radial_analysis_script.py:
# Computes the fluxes for various emission lines at specific PA (given by the user) for different radial distances.
# Plots the normalized to the maximum flux radial distributions of all the emission line. 
# The normalization is made by dividing the distribution with the maximum flux in the specific region 
# for r > r_mim (given by used)
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
import pyneb as pn


from satellite import radial_line_flux_script as rlfs
from satellite import find_maxvlaue_script as fms
from satellite import testing_unequal_arrays_script as tuas
from satellite import norm_flux_error_script as nfes



def radial_analysis(flux2D,flux2D_error,line_names,line_ext_error,lines_available,lines_radial,param_mod_name,param_model_values,par_plotname,par_plotymin,par_plotymax):
    
    class flux_radial:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]

    class flux_radial_error:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]
        
    class flux_radial_derred:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]

    class flux_radial_derred_error:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]
        
    
        
    file10 = open('general_output_file.txt','a')
    file2 = open('output_linesintensities_reddered_radial.txt','w')
    file22 = open('output_linesintensities_dereddered_radial.txt','w')
  

    index_pixel_scale=param_mod_name.index("pixel_scale")
    pixscale=float(param_model_values[index_pixel_scale])*0.01

    index_tot_pix=param_mod_name.index("total_num_pixels_horiz")
    index_slit_wid=param_mod_name.index("slit_width_in_spaxels")
    index_slit_len=param_mod_name.index("slit_length_in_spaxels")
    slitwidth=param_model_values[index_slit_wid]
    slitlength=param_model_values[index_slit_len]
    maxsize=param_model_values[index_tot_pix]
    
    
    index_angle_radial=param_mod_name.index("angle_for_radial_flux")    
    angle_for_radial=param_model_values[index_angle_radial]

    print >>file10, "#######################################################################################"
    print >>file10, "#######################################################################################"
    print >>file10, "#######################################################################################"
    print >>file10, "#######################################################################################"
    print >>file10, "#######################################################################################"
    print >>file10, "#######################################################################################"
    print >>file10, "#######################################################################################"
    print >>file10, "#######################################################################################"
    print >>file10, "#######################################################################################"
    fig,ax3=plt.subplots(figsize=(12,7))
    
    cHb,cHb_err,radial_dis_arc2=[],[],[]
    
    flux_radial_derred.Ha_6563,flux_radial_derred.Hb_4861,flux_radial_derred.Hg_4340,flux_radial_derred.Hd_4101,flux_radial_derred.HeIa_5876,flux_radial_derred.HeIb_6678,flux_radial_derred.HeIIa_4686,flux_radial_derred.HeIIb_5412,flux_radial_derred.OIa_5577,flux_radial_derred.OIb_6300,flux_radial_derred.OIc_6363,flux_radial_derred.OIIa_3727,flux_radial_derred.OIIb_3729,flux_radial_derred.OIIc_7320,flux_radial_derred.OIId_7330,flux_radial_derred.OIIIa_4363,flux_radial_derred.OIIIb_4959,flux_radial_derred.OIIIc_5007,flux_radial_derred.NI_5199,flux_radial_derred.NIIa_5755,flux_radial_derred.NIIb_6548,flux_radial_derred.NIIc_6584,flux_radial_derred.SIIa_6716,flux_radial_derred.SIIb_6731,flux_radial_derred.SIIIa_6312,flux_radial_derred.SIIIb_9069,flux_radial_derred.NeIIIa_3868,flux_radial_derred.NeIIIb_3967,flux_radial_derred.ArIII_7136,flux_radial_derred.ArIVa_4712,flux_radial_derred.ArIVb_4740,flux_radial_derred.ClIIIa_5517,flux_radial_derred.ClIIIb_5538,flux_radial_derred.CI_8727,flux_radial_derred.CII_6461=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
        
    flux_radial_derred_error.Ha_6563,flux_radial_derred_error.Hb_4861,flux_radial_derred_error.Hg_4340,flux_radial_derred_error.Hd_4101,flux_radial_derred_error.HeIa_5876,flux_radial_derred_error.HeIb_6678,flux_radial_derred_error.HeIIa_4686,flux_radial_derred_error.HeIIb_5412,flux_radial_derred_error.OIa_5577,flux_radial_derred_error.OIb_6300,flux_radial_derred_error.OIc_6363,flux_radial_derred_error.OIIa_3727,flux_radial_derred_error.OIIb_3729,flux_radial_derred_error.OIIc_7320,flux_radial_derred_error.OIId_7330,flux_radial_derred_error.OIIIa_4363,flux_radial_derred_error.OIIIb_4959,flux_radial_derred_error.OIIIc_5007,flux_radial_derred_error.NI_5199,flux_radial_derred_error.NIIa_5755,flux_radial_derred_error.NIIb_6548,flux_radial_derred_error.NIIc_6584,flux_radial_derred_error.SIIa_6716,flux_radial_derred_error.SIIb_6731,flux_radial_derred_error.SIIIa_6312,flux_radial_derred_error.SIIIb_9069,flux_radial_derred_error.NeIIIa_3868,flux_radial_derred_error.NeIIIb_3967,flux_radial_derred_error.ArIII_7136,flux_radial_derred_error.ArIVa_4712,flux_radial_derred_error.ArIVb_4740,flux_radial_derred_error.ClIIIa_5517,flux_radial_derred_error.ClIIIb_5538,flux_radial_derred_error.CI_8727,flux_radial_derred_error.CII_6461=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    
    
    index1=index2=0
    index_line=line_names.index("HI_6563s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.Ha_6563,flux_radial_error.Ha_6563,radial_dis_arc=rlfs.radialfluxline(flux2D.Ha_6563,flux2D_error.Ha_6563,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
 
        flux_radial.Ha_6563,flux_radial_error.Ha_6563,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.Ha_6563,flux_radial_error.Ha_6563,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("HI_6563e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraHa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.Ha_6563)):
                    flux_radial_error.Ha_6563[kk]=flux_radial_error.Ha_6563[kk]+flux_radial.Ha_6563[kk]*errextraHa
            else:
                errextraHa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.Ha_6563)):
                    flux_radial_error.Ha_6563[kk]=flux_radial.Ha_6563[kk]*errextraHa
                
            
        max_flux,rad_formax=fms.maxvalue(flux_radial.Ha_6563,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.Ha_6563/max_flux
        plt.plot(radial_dis_arc,radial,label="HI 6563")
        
        print >>file2,"dist ","Ha"
        for i in range(0,len(flux_radial.Ha_6563)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.Ha_6563[i])
            
        print >>file10, "Distance in arcsec for the peak of Ha      6563","%.3f" % rad_formax
        index1=1

    index_line=line_names.index("HI_4861s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.Hb_4861,flux_radial_error.Hb_4861,radial_dis_arc=rlfs.radialfluxline(flux2D.Hb_4861,flux2D_error.Hb_4861,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.Hb_4861,flux_radial_error.Hb_4861,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.Hb_4861,flux_radial_error.Hb_4861,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("HI_4861e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraHb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.Hb_4861)):
                    flux_radial_error.Hb_4861[kk]=flux_radial_error.Hb_4861[kk]+flux_radial.Hb_4861[kk]*errextraHb
            else:
                for kk in range(0,len(flux_radial.Hb_4861)):
                    flux_radial_error.Hb_4861[kk]=flux_radial_error.Hb_4861[kk]+flux_radial.Hb_4861[kk]*errextraHb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.Hb_4861,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.Hb_4861/max_flux
        plt.plot(radial_dis_arc,radial,label="HI 4861")
        
        print >>file2,"dist ","Hb"
        for i in range(0,len(flux_radial.Hb_4861)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.Hb_4861[i])
            
        print >>file10, "Distance in arcsec for the peak of Hb      4861","%.3f" % rad_formax
        index2=1


    index_line=line_names.index("HI_4340s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.Hg_4340,flux_radial_error.Hg_4340,radial_dis_arc=rlfs.radialfluxline(flux2D.Hg_4340,flux2D_error.Hg_4340,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.Hg_4340,flux_radial_error.Hg_4340,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.Hg_4340,flux_radial_error.Hg_4340,radial_dis_arc,pixscale,slitlength)
       
        index_extra_error=line_names.index("HI_4340e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraHg=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.Hg_4340)):
                    flux_radial_error.Hg_4340[kk]=flux_radial_error.Hg_4340[kk]+flux_radial.Hg_4340[kk]*errextraHg
            else:
                for kk in range(0,len(flux_radial.Hg_4340)):
                    flux_radial_error.Hg_4340[kk]=flux_radial_error.Hg_4340[kk]+flux_radial.Hg_4340[kk]*errextraHg
       
        max_flux,rad_formax=fms.maxvalue(flux_radial.Hg_4340,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.Hg_4340/max_flux
        plt.plot(radial_dis_arc,radial,label="HI 4340")
        
        print >>file2,"dist ","Hc"
        for i in range(0,len(flux_radial.Hg_4340)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.Hg_4340[i])
            
        print >>file10, "Distance in arcsec for the peak of Hg      4340","%.3f" % rad_formax

    index_line=line_names.index("HI_4101s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.Hd_4101,flux_radial_error.Hd_4101,radial_dis_arc=rlfs.radialfluxline(flux2D.Hd_4101,flux2D_error.Hd_4101,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)

        flux_radial.Hd_4101,flux_radial_error.Hd_4101,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.Hd_4101,flux_radial_error.Hd_4101,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("HI_4101e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraHd=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.Hd_4101)):
                    flux_radial_error.Hd_4101[kk]=flux_radial_error.Hd_4101[kk]+flux_radial.Hd_4101[kk]*errextraHd
            else:
                for kk in range(0,len(flux_radial.Hd_4101)):
                    flux_radial_error.Hd_4101[kk]=flux_radial_error.Hd_4101[kk]+flux_radial.Hd_4101[kk]*errextraHd
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.Hd_4101,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.Hd_4101/max_flux
        plt.plot(radial_dis_arc,radial,label="HI 4861")
        
        print >>file2,"dist ","Hd"
        for i in range(0,len(flux_radial.Hd_4101)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.Hd_4101[i])
        
        print >>file10, "Distance in arcsec for the peak of Hd      4101","%.3f" % rad_formax
        
    index_line=line_names.index("HeI_5876s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.HeIa_5876,flux_radial_error.HeIa_5876,radial_dis_arc=rlfs.radialfluxline(flux2D.HeIa_5876,flux2D_error.HeIa_5876,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.HeIa_5876,flux_radial_error.HeIa_5876,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.HeIa_5876,flux_radial_error.HeIa_5876,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("HeI_5876e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraHeIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.HeIa_5876)):
                    flux_radial_error.HeIa_5876[kk]=flux_radial_error.HeIa_5876[kk]+flux_radial.HeIa_5876[kk]*errextraHeIa
            else:
                for kk in range(0,len(flux_radial.HeIa_5876)):
                    flux_radial_error.HeIa_5876[kk]=flux_radial_error.HeIa_5876[kk]+flux_radial.HeIa_5876[kk]*errextraHeIa
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.HeIa_5876,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.HeIa_5876/max_flux
        plt.plot(radial_dis_arc,radial,label="He I 5876")

        print >>file2,"dist ","HeI 5876"
        for i in range(0,len(flux_radial.HeIa_5876)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.HeIa_5876[i])
            
        print >>file10, "Distance in arcsec for the peak of He I    5876","%.3f" % rad_formax    
        
    index_line=line_names.index("HeI_6678s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.HeIb_6678,flux_radial_error.HeIb_6678,radial_dis_arc=rlfs.radialfluxline(flux2D.HeIb_6678,flux2D_error.HeIb_6678,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.HeIb_6678,flux_radial_error.HeIb_6678,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.HeIb_6678,flux_radial_error.HeIb_6678,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("HeI_6678e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraHeIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.HeIb_6678)):
                    flux_radial_error.HeIb_6678[kk]=flux_radial_error.HeIb_6678[kk]+flux_radial.HeIb_6678[kk]*errextraHeIb
            else:
                for kk in range(0,len(flux_radial.HeIb_6678)):
                    flux_radial_error.HeIb_6678[kk]=flux_radial_error.HeIb_6678[kk]+flux_radial.HeIb_6678[kk]*errextraHeIb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.HeIb_6678,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.HeIb_6678/max_flux
        plt.plot(radial_dis_arc,radial,label="He I 6678")

        print >>file2,"dist ","HeI 6678"
        for i in range(0,len(flux_radial.HeIb_6678)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.HeIb_6678[i])
        
        print >>file10, "Distance in arcsec for the peak of He I    6678","%.3f" % rad_formax   
        
    index_line=line_names.index("HeII_4686s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.HeIIa_4686,flux_radial_error.HeIIa_4686,radial_dis_arc=rlfs.radialfluxline(flux2D.HeIIa_4686,flux2D_error.HeIIa_4686,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.HeIIa_4686,flux_radial_error.HeIIa_4686,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.HeIIa_4686,flux_radial_error.HeIIa_4686,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("HeII_4686e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraHeIIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.HeIIa_4686)):
                    flux_radial_error.HeIIa_4686[kk]=flux_radial_error.HeIIa_4686[kk]+flux_radial.HeIIa_4686[kk]*errextraHeIIa
            else:
                for kk in range(0,len(flux_radial.HeIIa_4686)):
                    flux_radial_error.HeIIa_4686[kk]=flux_radial_error.HeIIa_4686[kk]+flux_radial.HeIIa_4686[kk]*errextraHeIIa
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.HeIIa_4686,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.HeIIa_4686/max_flux
        plt.plot(radial_dis_arc,radial,label="He II 4686")
        
        print >>file2,"dist ","HeII 4686"
        for i in range(0,len(flux_radial.HeIIa_4686)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.HeIIa_4686[i])
        
        print >>file10, "Distance in arcsec for the peak of He II   4686","%.3f" % rad_formax     
        
    index_line=line_names.index("HeII_5412s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.HeIIb_5412,flux_radial_error.HeIIb_5412,radial_dis_arc=rlfs.radialfluxline(flux2D.HeIIb_5412,flux2D_error.HeIIb_5412,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.HeIIb_5412,flux_radial_error.HeIIb_5412,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.HeIIb_5412,flux_radial_error.HeIIb_5412,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("HeII_5412e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraHeIIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.HeIIb_5412)):
                    flux_radial_error.HeIIb_5412[kk]=flux_radial_error.HeIIb_5412[kk]+flux_radial.HeIIb_5412[kk]*errextraHeIIb
            else:
                for kk in range(0,len(flux_radial.HeIIb_5412)):
                    flux_radial_error.HeIIb_5412[kk]=flux_radial_error.HeIIb_5412[kk]+flux_radial.HeIIb_5412[kk]*errextraHeIIb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.HeIIb_5412,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.HeIIb_5412/max_flux
        plt.plot(radial_dis_arc,radial,label="He II 5412")
        
        print >>file2,"dist ","HeII 5412"
        for i in range(0,len(flux_radial.HeIIb_5412)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.HeIIb_5412[i])
        
        print >>file10, "Distance in arcsec for the peak of He II   5412","%.3f" % rad_formax
    
    index_line=line_names.index("N2_5755s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.NIIa_5755,flux_radial_error.NIIa_5755,radial_dis_arc=rlfs.radialfluxline(flux2D.NIIa_5755,flux2D_error.NIIa_5755,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale) 
        
        flux_radial.NIIa_5755,flux_radial_error.NIIa_5755,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.NIIa_5755,flux_radial_error.NIIa_5755,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("N2_5755e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraNIIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.NIIa_5755)):
                    flux_radial_error.NIIa_5755[kk]=flux_radial_error.NIIa_5755[kk]+flux_radial.NIIa_5755[kk]*errextraNIIa
            else:
                for kk in range(0,len(flux_radial.NIIa_5755)):
                    flux_radial_error.NIIa_5755[kk]=flux_radial_error.NIIa_5755[kk]+flux_radial.NIIa_5755[kk]*errextraNIIa
        
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.NIIa_5755,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.NIIa_5755/max_flux
        plt.plot(radial_dis_arc,radial,label="[N II] 5755")
        
        print >>file2,"dist ","[N II] 5755"
        for i in range(0,len(flux_radial.NIIa_5755)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.NIIa_5755[i])
        
        print >>file10, "Distance in arcsec for the peak of [N II]  5755","%.3f" % rad_formax
        
    index_line=line_names.index("N2_6548s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.NIIb_6548,flux_radial_error.NIIb_6548,radial_dis_arc=rlfs.radialfluxline(flux2D.NIIb_6548,flux2D_error.NIIb_6548,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale) 
        
        flux_radial.NIIb_6548,flux_radial_error.NIIb_6548,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.NIIb_6548,flux_radial_error.NIIb_6548,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("N2_6548e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraNIIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.NIIb_6548)):
                    flux_radial_error.NIIb_6548[kk]=flux_radial_error.NIIb_6548[kk]+flux_radial.NIIb_6548[kk]*errextraNIIb
            else:
                for kk in range(0,len(flux_radial.NIIb_6548)):
                    flux_radial_error.NIIb_6548[kk]=flux_radial_error.NIIb_6548[kk]+flux_radial.NIIb_6548[kk]*errextraNIIb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.NIIb_6548,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.NIIb_6548/max_flux
        plt.plot(radial_dis_arc,radial,label="[N II] 6548")
        
        print >>file2,"dist ","[N II] 6548"
        for i in range(0,len(flux_radial.NIIb_6548)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.NIIb_6548[i])
        
        print >>file10, "Distance in arcsec for the peak of [N II]  6548","%.3f" % rad_formax   
    
    index_line=line_names.index("N2_6583s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.NIIc_6584,flux_radial_error.NIIc_6584,radial_dis_arc=rlfs.radialfluxline(flux2D.NIIc_6584,flux2D_error.NIIc_6584,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale) 
        
        flux_radial.NIIc_6584,flux_radial_error.NIIc_6584,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.NIIc_6584,flux_radial_error.NIIc_6584,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("N2_6583e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraNIIc=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.NIIc_6584)):
                    flux_radial_error.NIIc_6584[kk]=flux_radial_error.NIIc_6584[kk]+flux_radial.NIIc_6584[kk]*errextraNIIc
            else:
                for kk in range(0,len(flux_radial.NIIc_6584)):
                    flux_radial_error.NIIc_6584[kk]=flux_radial_error.NIIc_6584[kk]+flux_radial.NIIc_6584[kk]*errextraNIIc
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.NIIc_6584,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.NIIc_6584/max_flux
        plt.plot(radial_dis_arc,radial,label="[N II] 6584")
        
        print >>file2,"dist ","[N II] 6584"
        for i in range(0,len(flux_radial.NIIc_6584)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.NIIc_6584[i])
        
        print >>file10, "Distance in arcsec for the peak of [N II]  6584","%.3f" % rad_formax
        
    index_line=line_names.index("N1_5199s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.NI_5199,flux_radial_error.NI_5199,radial_dis_arc=rlfs.radialfluxline(flux2D.NI_5199,flux2D_error.NI_5199,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.NI_5199,flux_radial_error.NI_5199,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.NI_5199,flux_radial_error.NI_5199,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("N1_5199e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraNI=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.NI_5199)):
                    flux_radial_error.NI_5199[kk]=flux_radial_error.NI_5199[kk]+flux_radial.NI_5199[kk]*errextraNI
            else:
                for kk in range(0,len(flux_radial.NI_5199)):
                    flux_radial_error.NI_5199[kk]=flux_radial_error.NI_5199[kk]+flux_radial.NI_5199[kk]*errextraNI
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.NI_5199,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.NI_5199/max_flux
        plt.plot(radial_dis_arc,radial,label="[N I] 5199")
        
        print >>file2,"dist ","[N I] 5199"
        for i in range(0,len(flux_radial.NI_5199)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.NI_5199[i])
        
        print >>file10, "Distance in arcsec for the peak of [N I]   5199","%.3f" % rad_formax
    
    index_line=line_names.index("O3_4363s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIIIa_4363,flux_radial_error.OIIIa_4363,radial_dis_arc=rlfs.radialfluxline(flux2D.OIIIa_4363,flux2D_error.OIIIa_4363,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIIIa_4363,flux_radial_error.OIIIa_4363,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIIIa_4363,flux_radial_error.OIIIa_4363,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O3_4363e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIIIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIIIa_4363)):
                    flux_radial_error.OIIIa_4363[kk]=flux_radial_error.OIIIa_4363[kk]+flux_radial.OIIIa_4363[kk]*errextraOIIIa
            else:
                for kk in range(0,len(flux_radial.OIIIa_4363)):
                    flux_radial_error.OIIIa_4363[kk]=flux_radial_error.OIIIa_4363[kk]+flux_radial.OIIIa_4363[kk]*errextraOIIIa
                    
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIIIa_4363,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIIIa_4363/max_flux
        plt.plot(radial_dis_arc,radial,label="[O III] 4363")
        
        print >>file2,"dist ","[O III] 4363"
        for i in range(0,len(flux_radial.OIIIa_4363)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIIIa_4363[i])
        
        print >>file10, "Distance in arcsec for the peak of [O III] 4363","%.3f" % rad_formax
        
    index_line=line_names.index("O3_4959s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIIIb_4959,flux_radial_error.OIIIb_4959,radial_dis_arc=rlfs.radialfluxline(flux2D.OIIIb_4959,flux2D_error.OIIIb_4959,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIIIb_4959,flux_radial_error.OIIIb_4959,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIIIb_4959,flux_radial_error.OIIIb_4959,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O3_4959e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIIIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIIIb_4959)):
                    flux_radial_error.OIIIb_4959[kk]=flux_radial_error.OIIIb_4959[kk]+flux_radial.OIIIb_4959[kk]*errextraOIIIb
            else:
                for kk in range(0,len(flux_radial.OIIIb_4959)):
                    flux_radial_error.OIIIb_4959[kk]=flux_radial_error.OIIIb_4959[kk]+flux_radial.OIIIb_4959[kk]*errextraOIIIb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIIIb_4959,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIIIb_4959/max_flux
        plt.plot(radial_dis_arc,radial,label="[O III] 4959")
        
        print >>file2,"dist ","[O III] 4959"
        for i in range(0,len(flux_radial.OIIIb_4959)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIIIb_4959[i])
        
        print >>file10, "Distance in arcsec for the peak of [O III] 4959","%.3f" % rad_formax
        
    index_line=line_names.index("O3_5007s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIIIc_5007,flux_radial_error.OIIIc_5007,radial_dis_arc=rlfs.radialfluxline(flux2D.OIIIc_5007,flux2D_error.OIIIc_5007,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIIIc_5007,flux_radial_error.OIIIc_5007,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIIIc_5007,flux_radial_error.OIIIc_5007,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O3_5007e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIIIc=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIIIc_5007)):
                    flux_radial_error.OIIIc_5007[kk]=flux_radial_error.OIIIc_5007[kk]+flux_radial.OIIIc_5007[kk]*errextraOIIIc
            else:
                for kk in range(0,len(flux_radial.OIIIc_5007)):
                    flux_radial_error.OIIIc_5007[kk]=flux_radial_error.OIIIc_5007[kk]+flux_radial.OIIIc_5007[kk]*errextraOIIIc
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIIIc_5007,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIIIc_5007/max_flux
        plt.plot(radial_dis_arc,radial,label="[O III] 5007")
        
        print >>file2,"dist ","[O III] 5007"
        for i in range(0,len(flux_radial.OIIIc_5007)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIIIc_5007[i])
        
        print >>file10, "Distance in arcsec for the peak of [O III] 5007","%.3f" % rad_formax
    
    index_line=line_names.index("O2_3727s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIIa_3727,flux_radial_error.OIIa_3727,radial_dis_arc=rlfs.radialfluxline(flux2D.OIIa_3727,flux2D_error.OIIa_3727,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIIa_3727,flux_radial_error.OIIa_3727,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIIa_3727,flux_radial_error.OIIa_3727,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O2_3727e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIIa_3727)):
                    flux_radial_error.OIIa_3727[kk]=flux_radial_error.OIIa_3727[kk]+flux_radial.OIIa_3727[kk]*errextraOIIa
            else:
                for kk in range(0,len(flux_radial.OIIa_3727)):
                    flux_radial_error.OIIa_3727[kk]=flux_radial_error.OIIa_3727[kk]+flux_radial.OIIa_3727[kk]*errextraOIIa
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIIa_3727,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIIa_3727/max_flux
        plt.plot(radial_dis_arc,radial,label="[O II] 3727")
        
        print >>file2,"dist ","[O II] 3727"
        for i in range(0,len(flux_radial.OIIa_3727)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIIa_3727[i])
        
        print >>file10, "Distance in arcsec for the peak of [O II]  3727","%.3f" % rad_formax

    index_line=line_names.index("O2_3729s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIIb_3729,flux_radial_error.OIIb_3729,radial_dis_arc=rlfs.radialfluxline(flux2D.OIIb_3729,flux2D_error.OIIb_3729,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIIb_3729,flux_radial_error.OIIb_3729,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIIb_3729,flux_radial_error.OIIb_3729,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O2_3729e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIIb_3729)):
                    flux_radial_error.OIIb_3729[kk]=flux_radial_error.OIIb_3729[kk]+flux_radial.OIIb_3729[kk]*errextraOIIb
            else:
                for kk in range(0,len(flux_radial.OIIb_3729)):
                    flux_radial_error.OIIb_3729[kk]=flux_radial_error.OIIb_3729[kk]+flux_radial.OIIb_3729[kk]*errextraOIIb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIIb_3729,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIIb_3729/max_flux
        plt.plot(radial_dis_arc,radial,label="[O II] 3729")
        
        print >>file2,"dist ","[O II] 3729"
        for i in range(0,len(flux_radial.OIIb_3729)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIIb_3729[i])
        
        print >>file10, "Distance in arcsec for the peak of [O II]  3729","%.3f" % rad_formax
        
    index_line=line_names.index("O2_7320s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIIc_7320,flux_radial_error.OIIc_7320,radial_dis_arc=rlfs.radialfluxline(flux2D.OIIc_7320,flux2D_error.OIIc_7320,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIIc_7320,flux_radial_error.OIIc_7320,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIIc_7320,flux_radial_error.OIIc_7320,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O2_7320e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIIc=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIIc_7320)):
                    flux_radial_error.OIIc_7320[kk]=flux_radial_error.OIIc_7320[kk]+flux_radial.OIIc_7320[kk]*errextraOIIc
            else:
                for kk in range(0,len(flux_radial.OIIc_7320)):
                    flux_radial_error.OIIc_7320[kk]=flux_radial_error.OIIc_7320[kk]+flux_radial.OIIc_7320[kk]*errextraOIIc
        
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIIc_7320,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIIc_7320/max_flux
        plt.plot(radial_dis_arc,radial,label="[O II] 7320")
        
        print >>file2,"dist ","[O II] 7320"
        for i in range(0,len(flux_radial.OIIc_7320)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIIc_7320[i])
        
        print >>file10, "Distance in arcsec for the peak of [O II]  7320","%.3f" % rad_formax
        
    index_line=line_names.index("O2_7330s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIId_7330,flux_radial_error.OIId_7330,radial_dis_arc=rlfs.radialfluxline(flux2D.OIId_7330,flux2D_error.OIId_7330,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIId_7330,flux_radial_error.OIId_7330,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIId_7330,flux_radial_error.OIId_7330,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O2_7330e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIId=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIId_7330)):
                    flux_radial_error.OIId_7330[kk]=flux_radial_error.OIId_7330[kk]+flux_radial.OIId_7330[kk]*errextraOIId
            else:
                for kk in range(0,len(flux_radial.OIId_7330)):
                    flux_radial_error.OIId_7330[kk]=flux_radial_error.OIId_7330[kk]+flux_radial.OIId_7330[kk]*errextraOIId
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIId_7330,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIId_7330/max_flux
        plt.plot(radial_dis_arc,radial,label="[O II] 7330")
        
        print >>file2,"dist ","[O II] 7330"
        for i in range(0,len(flux_radial.OIId_7330)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIId_7330[i])
        
        print >>file10, "Distance in arcsec for the peak of [O II]  7330","%.3f" % rad_formax   
    
    index_line=line_names.index("O1_5577s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIa_5577,flux_radial_error.OIa_5577,radial_dis_arc=rlfs.radialfluxline(flux2D.OIa_5577,flux2D_error.OIa_5577,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIa_5577,flux_radial_error.OIa_5577,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIa_5577,flux_radial_error.OIa_5577,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O1_5577s")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIa_5577)):
                    flux_radial_error.OIa_5577[kk]=flux_radial_error.OIa_5577[kk]+flux_radial.OIa_5577[kk]*errextraOIa
            else:
                for kk in range(0,len(flux_radial.OIa_5577)):
                    flux_radial_error.OIa_5577[kk]=flux_radial_error.OIa_5577[kk]+flux_radial.OIa_5577[kk]*errextraOIa
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIa_5577,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIa_5577/max_flux
        plt.plot(radial_dis_arc,radial,label="[O I] 5577")
        
        print >>file2,"dist ","[O I] 5577"
        for i in range(0,len(flux_radial.OIa_5577)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIa_5577[i])
        
        print >>file10, "Distance in arcsec for the peak of [O I]   5577","%.3f" % rad_formax

    index_line=line_names.index("O1_6300s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIb_6300,flux_radial_error.OIb_6300,radial_dis_arc=rlfs.radialfluxline(flux2D.OIb_6300,flux2D_error.OIb_6300,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIb_6300,flux_radial_error.OIb_6300,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIb_6300,flux_radial_error.OIb_6300,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O1_6300e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIb_6300)):
                    flux_radial_error.OIb_6300[kk]=flux_radial_error.OIb_6300[kk]+flux_radial.OIb_6300[kk]*errextraOIb
            else:
                for kk in range(0,len(flux_radial.OIb_6300)):
                    flux_radial_error.OIb_6300[kk]=flux_radial_error.OIb_6300[kk]+flux_radial.OIb_6300[kk]*errextraOIb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIb_6300,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIb_6300/max_flux
        plt.plot(radial_dis_arc,radial,label="[O I] 6300")
        
        print >>file2,"dist ","[O I] 6300"
        for i in range(0,len(flux_radial.OIb_6300)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIb_6300[i])
        
        print >>file10, "Distance in arcsec for the peak of [O I]   6300","%.3f" % rad_formax
            
    index_line=line_names.index("O1_6363s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.OIc_6363,flux_radial_error.OIc_6363,radial_dis_arc=rlfs.radialfluxline(flux2D.OIc_6363,flux2D_error.OIc_6363,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.OIc_6363,flux_radial_error.OIc_6363,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.OIc_6363,flux_radial_error.OIc_6363,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("O1_6363e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraOIc=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.OIc_6363)):
                    flux_radial_error.OIc_6363[kk]=flux_radial_error.OIc_6363[kk]+flux_radial.OIc_6363[kk]*errextraOIc
            else:
                for kk in range(0,len(flux_radial.OIc_6363)):
                    flux_radial_error.OIc_6363[kk]=flux_radial_error.OIc_6363[kk]+flux_radial.OIc_6363[kk]*errextraOIc
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.OIc_6363,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.OIc_6363/max_flux
        plt.plot(radial_dis_arc,radial,label="[O I] 6363")
        
        print >>file2,"dist ","[O I] 6363"
        for i in range(0,len(flux_radial.OIc_6363)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.OIc_6363[i])
        
        print >>file10, "Distance in arcsec for the peak of [O I]   6363","%.3f" % rad_formax

    index_line=line_names.index("S2_6716s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.SIIa_6716,flux_radial_error.SIIa_6716,radial_dis_arc=rlfs.radialfluxline(flux2D.SIIa_6716,flux2D_error.SIIa_6716,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.SIIa_6716,flux_radial_error.SIIa_6716,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.SIIa_6716,flux_radial_error.SIIa_6716,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("S2_6716e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraSIIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.SIIa_6716)):
                    flux_radial_error.SIIa_6716[kk]=flux_radial_error.SIIa_6716[kk]+flux_radial.SIIa_6716[kk]*errextraSIIa
            else:
                for kk in range(0,len(flux_radial.SIIa_6716)):
                    flux_radial_error.SIIa_6716[kk]=flux_radial_error.SIIa_6716[kk]+flux_radial.SIIa_6716[kk]*errextraSIIa
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.SIIa_6716,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.SIIa_6716/max_flux
        plt.plot(radial_dis_arc,radial,label="[S II] 6716")
        
        print >>file2,"dist ","[S II] 6716"
        for i in range(0,len(flux_radial.SIIa_6716)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.SIIa_6716[i])
        
        print >>file10, "Distance in arcsec for the peak of [S II]  6716","%.3f" % rad_formax
        
    index_line=line_names.index("S2_6731s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.SIIb_6731,flux_radial_error.SIIb_6731,radial_dis_arc=rlfs.radialfluxline(flux2D.SIIb_6731,flux2D_error.SIIb_6731,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.SIIb_6731,flux_radial_error.SIIb_6731,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.SIIb_6731,flux_radial_error.SIIb_6731,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("S2_6731e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraSIIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.SIIb_6731)):
                    flux_radial_error.SIIb_6731[kk]=flux_radial_error.SIIb_6731[kk]+flux_radial.SIIb_6731[kk]*errextraSIIb
            else:
                for kk in range(0,len(flux_radial.SIIb_6731)):
                    flux_radial_error.SIIb_6731[kk]=flux_radial_error.SIIb_6731[kk]+flux_radial.SIIb_6731[kk]*errextraSIIb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.SIIb_6731,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.SIIb_6731/max_flux
        plt.plot(radial_dis_arc,radial,label="[S II] 6731")
        
        print >>file2,"dist ","[S II] 6731"
        for i in range(0,len(flux_radial.SIIb_6731)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.SIIb_6731[i])
        
        print >>file10, "Distance in arcsec for the peak of [S II]  6731","%.3f" % rad_formax

    index_line=line_names.index("S3_6312s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.SIIIa_6312,flux_radial_error.SIIIa_6312,radial_dis_arc=rlfs.radialfluxline(flux2D.SIIIa_6312,flux2D_error.SIIIa_6312,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.SIIIa_6312,flux_radial_error.SIIIa_6312,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.SIIIa_6312,flux_radial_error.SIIIa_6312,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("S3_6312e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraSIIIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.SIIIa_6312)):
                    flux_radial_error.SIIIa_6312[kk]=flux_radial_error.SIIIa_6312[kk]+flux_radial.SIIIa_6312[kk]*errextraSIIIa
            else:
                for kk in range(0,len(flux_radial.SIIIa_6312)):
                    flux_radial_error.SIIIa_6312[kk]=flux_radial_error.SIIIa_6312[kk]+flux_radial.SIIIa_6312[kk]*errextraSIIIa
                    
        max_flux,rad_formax=fms.maxvalue(flux_radial.SIIIa_6312,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.SIIIa_6312/max_flux
        plt.plot(radial_dis_arc,radial,label="[S III] 6312")
        
        print >>file2,"dist ","[S III] 6312"
        for i in range(0,len(flux_radial.SIIIa_6312)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.SIIIa_6312[i])
        
        print >>file10, "Distance in arcsec for the peak of [S III] 6312","%.3f" % rad_formax
        
    index_line=line_names.index("S3_9069s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.SIIIb_9069,flux_radial_error.SIIIb_9069,radial_dis_arc=rlfs.radialfluxline(flux2D.SIIIb_9069,flux2D_error.SIIIb_9069,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.SIIIb_9069,flux_radial_error.SIIIb_9069,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.SIIIb_9069,flux_radial_error.SIIIb_9069,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("S3_9069e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraSIIIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.SIIIb_9069)):
                    flux_radial_error.SIIIb_9069[kk]=flux_radial_error.SIIIb_9069[kk]+flux_radial.SIIIb_9069[kk]*errextraSIIIb
            else:
                for kk in range(0,len(flux_radial.SIIIb_9069)):
                    flux_radial_error.SIIIb_9069[kk]=flux_radial_error.SIIIb_9069[kk]+flux_radial.SIIIb_9069[kk]*errextraSIIIb
                    
        max_flux,rad_formax=fms.maxvalue(flux_radial.SIIIb_9069,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.SIIIb_9069/max_flux
        plt.plot(radial_dis_arc,radial,label="[S III] 9069")
        
        print >>file2,"dist ","[S III] 9069"
        for i in range(0,len(flux_radial.SIIIb_9069)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.SIIIb_9069[i])
        
        print >>file10, "Distance in arcsec for the peak of [S III] 9069","%.3f" % rad_formax

    index_line=line_names.index("Cl3_5517s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.ClIIIa_5517,flux_radial_error.ClIIIa_5517,radial_dis_arc=rlfs.radialfluxline(flux2D.ClIIIa_5517,flux2D_error.ClIIIa_5517,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.ClIIIa_5517,flux_radial_error.ClIIIa_5517,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.ClIIIa_5517,flux_radial_error.ClIIIa_5517,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("Cl3_5517e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraClIIIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.ClIIIa_5517)):
                    flux_radial_error.ClIIIa_5517[kk]=flux_radial_error.ClIIIa_5517[kk]+flux_radial.ClIIIa_5517[kk]*errextraClIIIa
            else:
                for kk in range(0,len(flux_radial.ClIIIa_5517)):
                    flux_radial_error.ClIIIa_5517[kk]=flux_radial_error.ClIIIa_5517[kk]+flux_radial.ClIIIa_5517[kk]*errextraClIIIa
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.ClIIIa_5517,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.ClIIIa_5517/max_flux
        plt.plot(radial_dis_arc,radial,label="[Cl III] 5517")
        
        print >>file2,"dist ","[Cl III] 5517"
        for i in range(0,len(flux_radial.ClIIIa_5517)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.ClIIIa_5517[i])
        
        print >>file10, "Distance in arcsec for the peak of [Cl III] 5517","%.3f" % rad_formax

    index_line=line_names.index("Cl3_5538s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.ClIIIb_5538,flux_radial_error.ClIIIb_5538,radial_dis_arc=rlfs.radialfluxline(flux2D.ClIIIb_5538,flux2D_error.ClIIIb_5538,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.ClIIIb_5538,flux_radial_error.ClIIIb_5538,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.ClIIIb_5538,flux_radial_error.ClIIIb_5538,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("Cl3_5538e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraClIIIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.ClIIIb_5538)):
                    flux_radial_error.ClIIIb_5538[kk]=flux_radial_error.ClIIIb_5538[kk]+flux_radial.ClIIIb_5538[kk]*errextraClIIIb
            else:
                for kk in range(0,len(flux_radial.ClIIIb_5538)):
                    flux_radial_error.ClIIIb_5538[kk]=flux_radial_error.ClIIIb_5538[kk]+flux_radial.ClIIIb_5538[kk]*errextraClIIIb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.ClIIIb_5538,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.ClIIIb_5538/max_flux
        plt.plot(radial_dis_arc,radial,label="[Cl III] 5538")
        
        print >>file2,"dist ","[Cl III] 5538"
        for i in range(0,len(flux_radial.ClIIIb_5538)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.ClIIIb_5538[i])
        
        print >>file10, "Distance in arcsec for the peak of [Cl III] 5538","%.3f" % rad_formax
        
    index_line=line_names.index("Ar3_7136s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.ArIII_7136,flux_radial_error.ArIII_7136,radial_dis_arc=rlfs.radialfluxline(flux2D.ArIII_7136,flux2D_error.ArIII_7136,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.ArIII_7136,flux_radial_error.ArIII_7136,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.ArIII_7136,flux_radial_error.ArIII_7136,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("Ar3_7136e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraArIII=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.ArIII_7136)):
                    flux_radial_error.ArIII_7136[kk]=flux_radial_error.ArIII_7136[kk]+flux_radial.ArIII_7136[kk]*errextraArIII
            else:
                for kk in range(0,len(flux_radial.ArIII_7136)):
                    flux_radial_error.ArIII_7136[kk]=flux_radial_error.ArIII_7136[kk]+flux_radial.ArIII_7136[kk]*errextraArIII
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.ArIII_7136,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.ArIII_7136/max_flux
        plt.plot(radial_dis_arc,radial,label="[Ar III] 7136")
        
        print >>file2,"dist ","[Ar III] 7136"
        for i in range(0,len(flux_radial.ArIII_7136)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.ArIII_7136[i])
        
        print >>file10, "Distance in arcsec for the peak of [Ar III] 7136","%.3f" % rad_formax    
        
    index_line=line_names.index("Ar4_4712s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.ArIVa_4712,flux_radial_error.ArIVa_4712,radial_dis_arc=rlfs.radialfluxline(flux2D.ArIVa_4712,flux2D_error.ArIVa_4712,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.ArIVa_4712,flux_radial_error.ArIVa_4712,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.ArIVa_4712,flux_radial_error.ArIVa_4712,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("Ar4_4712e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraArIVa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.ArIVa_4712)):
                    flux_radial_error.ArIVa_4712[kk]=flux_radial_error.ArIVa_4712[kk]+flux_radial.ArIVa_4712[kk]*errextraArIVa
            else:
                for kk in range(0,len(flux_radial.ArIVa_4712)):
                    flux_radial_error.ArIVa_4712[kk]=flux_radial_error.ArIVa_4712[kk]+flux_radial.ArIVa_4712[kk]*errextraArIVa
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.ArIVa_4712,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.ArIVa_4712/max_flux
        plt.plot(radial_dis_arc,radial,label="[Ar IV] 4712")
        
        print >>file2,"dist ","[Ar IV] 4712"
        for i in range(0,len(flux_radial.ArIVa_4712)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.ArIVa_4712[i]) 
        
        print >>file10, "Distance in arcsec for the peak of [Ar IV] 4712","%.3f" % rad_formax 
        
    index_line=line_names.index("Ar4_4740s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.ArIVb_4740,flux_radial_error.ArIVb_4740,radial_dis_arc=rlfs.radialfluxline(flux2D.ArIVb_4740,flux2D_error.ArIVb_4740,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.ArIVb_4740,flux_radial_error.ArIVb_4740,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.ArIVb_4740,flux_radial_error.ArIVb_4740,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("Ar4_4740e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraArIVb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.ArIVb_4740)):
                    flux_radial_error.ArIVb_4740[kk]=flux_radial_error.ArIVb_4740[kk]+flux_radial.ArIVb_4740[kk]*errextraArIVb
            else:
                for kk in range(0,len(flux_radial.ArIVb_4740)):
                    flux_radial_error.ArIVb_4740[kk]=flux_radial_error.ArIVb_4740[kk]+flux_radial.ArIVb_4740[kk]*errextraArIVb
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.ArIVb_4740,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.ArIVb_4740/max_flux
        plt.plot(radial_dis_arc,radial,label="[Ar IV] 4740")
        
        print >>file2,"dist ","[Ar IV] 4740"
        for i in range(0,len(flux_radial.ArIVb_4740)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.ArIVb_4740[i])
        
        print >>file10, "Distance in arcsec for the peak of [Ar IV] 4740","%.3f" % rad_formax    
        
    index_line=line_names.index("Ne3_3868s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.NeIIIa_3868,flux_radial_error.NeIIIa_3868,radial_dis_arc=rlfs.radialfluxline(flux2D.NeIIIa_3868,flux2D_error.NeIIIa_3868,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.NeIIIa_3868,flux_radial_error.NeIIIa_3868,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.NeIIIa_3868,flux_radial_error.NeIIIa_3868,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("Ne3_3868e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraNeIIIa=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.NeIIIa_3868)):
                    flux_radial_error.NeIIIa_3868[kk]=flux_radial_error.NeIIIa_3868[kk]+flux_radial.NeIIIa_3868[kk]*errextraNeIIIa
            else:
                for kk in range(0,len(flux_radial.NeIIIa_3868)):
                    flux_radial_error.NeIIIa_3868[kk]=flux_radial_error.NeIIIa_3868[kk]+flux_radial.NeIIIa_3868[kk]*errextraNeIIIa
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.NeIIIa_3868,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.NeIIIa_3868/max_flux
        plt.plot(radial_dis_arc,radial,label="[Ne III] 3868")
        
        print >>file2,"dist ","[Ne III] 3868"
        for i in range(0,len(flux_radial.NeIIIa_3868)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.NeIIIa_3868[i])
        
        print >>file10, "Distance in arcsec for the peak of [Ne III] 3868","%.3f" % rad_formax
    
    index_line=line_names.index("Ne3_3967s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.NeIIIb_3967,flux_radial_error.NeIIIb_3967,radial_dis_arc=rlfs.radialfluxline(flux2D.NeIIIb_3967,flux2D_error.NeIIIb_3967,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.NeIIIb_3967,flux_radial_error.NeIIIb_3967,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.NeIIIb_3967,flux_radial_error.NeIIIb_3967,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("Ne3_3967e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraNeIIIb=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.NeIIIb_3967)):
                    flux_radial_error.NeIIIb_3967[kk]=flux_radial_error.NeIIIb_3967[kk]+flux_radial.NeIIIb_3967[kk]*errextraNeIIIb
            else:
                for kk in range(0,len(flux_radial.NeIIIb_3967)):
                    flux_radial_error.NeIIIb_3967[kk]=flux_radial_error.NeIIIb_3967[kk]+flux_radial.NeIIIb_3967[kk]*errextraNeIIIb
        
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.NeIIIb_3967,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.NeIIIb_3967/max_flux
        plt.plot(radial_dis_arc,radial,label="[Ne III] 3967")
        
        print >>file2,"dist ","[Ne III] 3967"
        for i in range(0,len(flux_radial.NeIIIb_3967)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.NeIIIb_3967[i])
        
        print >>file10, "Distance in arcsec for the peak of [Ne III] 3967","%.3f" % rad_formax
     
    index_line=line_names.index("C1_8727s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.CI_8727,flux_radial_error.CI_8727,radial_dis_arc=rlfs.radialfluxline(flux2D.CI_8727,flux2D_error.CI_8727,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.CI_8727,flux_radial_error.CI_8727,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.CI_8727,flux_radial_error.CI_8727,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("C1_8727e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraCI=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.CI_8727)):
                    flux_radial_error.CI_8727[kk]=flux_radial_error.CI_8727[kk]+flux_radial.CI_8727[kk]*errextraCI
            else:
                for kk in range(0,len(flux_radial.CI_8727)):
                    flux_radial_error.CI_8727[kk]=flux_radial_error.CI_8727[kk]+flux_radial.CI_8727[kk]*errextraCI
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.CI_8727,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.CI_8727/max_flux
        plt.plot(radial_dis_arc,radial,label="[C I] 8727")
        
        print >>file2,"dist ","[C I] 8727"
        for i in range(0,len(flux_radial.CI_8727)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.CI_8727[i])
        
        print >>file10, "Distance in arcsec for the peak of [C I]    8727","%.3f" % rad_formax 

    index_line=line_names.index("C2_6461s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        flux_radial.CII_6461,flux_radial_error.CII_6461,radial_dis_arc=rlfs.radialfluxline(flux2D.CII_6461,flux2D_error.CII_6461,angle_for_radial,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize,pixscale)
        
        flux_radial.CII_6461,flux_radial_error.CII_6461,radial_dis_arc=tuas.testing_unequal_arrays(flux_radial.CII_6461,flux_radial_error.CII_6461,radial_dis_arc,pixscale,slitlength)
        
        index_extra_error=line_names.index("C2_6461e")          
        if lines_available[index_extra_error]=="yes":
            if line_ext_error[index_extra_error]!=0:
                errextraCII=float(line_ext_error[index_extra_error])
                for kk in range(0,len(flux_radial.CII_6461)):
                    flux_radial_error.CII_6461[kk]=flux_radial_error.CII_6461[kk]+flux_radial.CII_6461[kk]*errextraCII
            else:
                for kk in range(0,len(flux_radial.CII_6461)):
                    flux_radial_error.CII_6461[kk]=flux_radial_error.CII_6461[kk]+flux_radial.CII_6461[kk]*errextraCII
        
        max_flux,rad_formax=fms.maxvalue(flux_radial.CII_6461,radial_dis_arc,param_mod_name,param_model_values,pixscale)
        radial=flux_radial.CII_6461/max_flux
        plt.plot(radial_dis_arc,radial,label="[C II] 6461")
        
        print >>file2,"dist ","[C II] 6461"
        for i in range(0,len(flux_radial.CII_6461)):
            print >>file2, '{:03.3f}'.format(radial_dis_arc[i]),'{:06.3f}'.format(flux_radial.CII_6461[i])
        
        print >>file10, "Distance in arcsec for the peak of [C II]   6461","%.3f" % rad_formax  
     
    print >>file10, "#######################################################################################"
    
    max_len=max(len(flux_radial.Ha_6563),len(flux_radial.Hb_4861),len(flux_radial.Hg_4340),len(flux_radial.Hd_4101),len(flux_radial.HeIa_5876),len(flux_radial.HeIb_6678),len(flux_radial.HeIIa_4686),len(flux_radial.HeIIb_5412),len(flux_radial.OIa_5577),len(flux_radial.OIb_6300),len(flux_radial.OIc_6363),len(flux_radial.OIIa_3727),len(flux_radial.OIIb_3729),len(flux_radial.OIIc_7320),len(flux_radial.OIId_7330),len(flux_radial.OIIIa_4363),len(flux_radial.OIIIb_4959),len(flux_radial.OIIIc_5007),len(flux_radial.NI_5199),len(flux_radial.NIIa_5755),len(flux_radial.NIIb_6548),len(flux_radial.NIIc_6584),len(flux_radial.SIIa_6716),len(flux_radial.SIIb_6731),len(flux_radial.SIIIa_6312),len(flux_radial.SIIIb_9069),len(flux_radial.NeIIIa_3868),len(flux_radial.NeIIIb_3967),len(flux_radial.ArIII_7136),len(flux_radial.ArIVa_4712),len(flux_radial.ArIVb_4740),len(flux_radial.ClIIIa_5517),len(flux_radial.ClIIIb_5538),len(flux_radial.CI_8727),len(flux_radial.CII_6461))
    
    radial_dis_arc3=[]
    jrad=0    
    for kk in range(0,max_len):
        
        if flux_radial.Hb_4861[jrad]>0 and flux_radial.Ha_6563[jrad]>0:
            
            file3 = open('test.dat','w')
            print >>file3, "LINE test err"

            radial_dis_arc3.append(0)
            radial_dis_arc3[jrad]=radial_dis_arc[jrad]
            
            index_line=line_names.index("HI_6563s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.Ha_6563[jrad],flux_radial_error.Ha_6563[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "H1r_6563A", '{:03.3f}'.format(flux_radial.Ha_6563[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("HI_4861s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "H1r_4861A", '{:03.3f}'.format(flux_radial.Hb_4861[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("HI_4340s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.Hg_4340[jrad],flux_radial_error.Hg_4340[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "H1r_4341A", '{:03.3f}'.format(flux_radial.Hg_4340[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("HI_4101s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.Hd_4101[jrad],flux_radial_error.Hd_4101[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "H1r_4102A", '{:03.3f}'.format(flux_radial.Hg_4340[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("HeII_4686s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.HeIIa_4686[jrad],flux_radial_error.HeIIa_4686[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "He2r_4686A", '{:03.3f}'.format(flux_radial.HeIIa_4686[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("HeII_5412s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.HeIIb_5412[jrad],flux_radial_error.HeIIb_5412[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "He2r_5411A", '{:03.3f}'.format(flux_radial.HeIIb_5412[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("HeI_5876s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.HeIa_5876[jrad],flux_radial_error.HeIa_5876[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "He1r_5876A", '{:03.3f}'.format(flux_radial.HeIa_5876[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("HeI_6678s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.HeIb_6678[jrad],flux_radial_error.HeIb_6678[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "He1r_6678A", '{:03.3f}'.format(flux_radial.HeIb_6678[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("S2_6716s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.SIIa_6716[jrad],flux_radial_error.SIIa_6716[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "S2_6716A", '{:03.3f}'.format(flux_radial.SIIa_6716[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("S2_6731s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.SIIb_6731[jrad],flux_radial_error.SIIb_6731[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "S2_6731A", '{:03.3f}'.format(flux_radial.SIIb_6731[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("S3_6312s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.SIIIa_6312[jrad],flux_radial_error.SIIIa_6312[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "S3_6312A", '{:03.3f}'.format(flux_radial.SIIIa_6312[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("S3_9069s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.SIIIb_9069[jrad],flux_radial_error.SIIIb_9069[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "S3_9069A", '{:03.3f}'.format(flux_radial.SIIIb_9069[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("N1_5199s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.NI_5199[jrad],flux_radial_error.NI_5199[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "N1_5200A", '{:03.3f}'.format(flux_radial.NI_5199[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("N2_5755s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.NIIa_5755[jrad],flux_radial_error.NIIa_5755[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "N2_5755A", '{:03.3f}'.format(flux_radial.NIIa_5755[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("N2_6548s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.NIIb_6548[jrad],flux_radial_error.NIIb_6548[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "N2_6548A", '{:03.3f}'.format(flux_radial.NIIb_6548[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("N2_6583s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.NIIc_6584[jrad],flux_radial_error.NIIc_6584[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "N2_6584A", '{:03.3f}'.format(flux_radial.NIIc_6584[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O1_5577s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.v[jrad],flux_radial_error.OIa_5577[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O1_5577A", '{:03.3f}'.format(flux_radial.OIa_5577[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O1_6300s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIb_6300[jrad],flux_radial_error.OIb_6300[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O1_6300A", '{:03.3f}'.format(flux_radial.OIb_6300[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O1_6363s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIc_6363[jrad],flux_radial_error.OIc_6363[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O1_6364A", '{:03.3f}'.format(flux_radial.OIc_6363[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O2_3727s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIIa_3727[jrad],flux_radial_error.OIIa_3727[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O2_3726A", '{:03.3f}'.format(flux_radial.OIIa_3727[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O2_3729s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIIb_3729[jrad],flux_radial_error.OIIb_3729[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O2_3729A", '{:03.3f}'.format(flux_radial.OIIb_3729[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O2_7320s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIIc_7320[jrad],flux_radial_error.OIIc_7320[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O2_7319A+", '{:03.3f}'.format(flux_radial.OIIc_7320[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O2_7330s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIId_7330[jrad],flux_radial_error.OIId_7330[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O2_7330A+", '{:03.3f}'.format(flux_radial.OIId_7330[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O3_4363s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIIIa_4363[jrad],flux_radial_error.OIIIa_4363[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O3_4363A", '{:03.3f}'.format(flux_radial.OIIIa_4363[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O3_4959s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIIIb_4959[jrad],flux_radial_error.OIIIb_4959[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O3_4959A", '{:03.3f}'.format(flux_radial.OIIIb_4959[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("O3_5007s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.OIIIc_5007[jrad],flux_radial_error.OIIIc_5007[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "O3_5007A", '{:03.3f}'.format(flux_radial.OIIIc_5007[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("Cl3_5517s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.ClIIIa_5517[jrad],flux_radial_error.ClIIIa_5517[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "Cl3_5518A", '{:03.3f}'.format(flux_radial.ClIIIa_5517[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("Cl3_5517s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.ClIIIb_5538[jrad],flux_radial_error.ClIIIb_5538[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "Cl3_5538A", '{:03.3f}'.format(flux_radial.ClIIIb_5538[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)     
            
            index_line=line_names.index("Ne3_3868s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.NeIIIa_3868[jrad],flux_radial_error.NeIIIa_3868[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "Ne3_3869A", '{:03.3f}'.format(flux_radial.NeIIIa_3868[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("Ne3_3967s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.NeIIIb_3967[jrad],flux_radial_error.NeIIIb_3967[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "Ne3_3968A", '{:03.3f}'.format(flux_radial.NeIIIb_3967[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)               
            
            index_line=line_names.index("Ar3_7136s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.ArIII_7136[jrad],flux_radial_error.ArIII_7136[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "Ar3_7136A", '{:03.3f}'.format(flux_radial.ArIII_7136[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("Ar4_4712s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.ArIVa_4712[jrad],flux_radial_error.ArIVa_4712[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "Ar4_4711A", '{:03.3f}'.format(flux_radial.ArIVa_4712[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("Ar4_4740s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.ArIVb_4740[jrad],flux_radial_error.ArIVb_4740[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "Ar4_4740A", '{:03.3f}'.format(flux_radial.ArIVb_4740[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            index_line=line_names.index("C1_8727s")
            if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
                error=nfes.flux_Error(flux_radial.CI_8727[jrad],flux_radial_error.CI_8727[jrad],flux_radial.Hb_4861[jrad],flux_radial_error.Hb_4861[jrad])
                print >>file3, "C1_8728A", '{:03.3f}'.format(flux_radial.CI_8727[jrad]*100/flux_radial.Hb_4861[jrad]),'{:03.3f}'.format(error)
            
            
            file3.close() 
                
            index_extlaw=11
            name_ext=str(param_mod_name[index_extlaw])
            if name_ext=='CCM89': name_ext=str('CCM89')
            if name_ext=='K76': name_ext=str('K76')
            if name_ext=='F99': name_ext=str('F99')
            if name_ext=='S79_H83_CCM89': name_ext=str('S79') + ' ' + str('H83') + ' ' + str('CCM89')
            if name_ext=='CCM89_Bal07': name_ext=str('CCM89') + ' ' + str('Bal07')
            if name_ext=='CCM89_oD94': name_ext=str('CCM89') + ' ' + str('oD94')
            if name_ext=='SM79_Gal': name_ext=str('SM79') + ' ' + str('Gal')
            if name_ext=='G03_LMC': name_ext=str('G03') + ' ' + str('LMC')
            if name_ext=='MCC99_FM90_LMC': name_ext=str('MCC99') + ' ' + str('FM90') + ' ' + str('LMC')
            #if name_ext=='F99_like': name_ext=str('F99') + ' ' + str('-') + ' ' + str('like')
            if name_ext=='F88_F99_LMC': name_ext=str('F88') + ' ' + str('F99') + ' ' + str('LMC')
            if name_ext=='No_correction': name_ext=str('No') + ' ' + str('correction')
            RV_ext=float(param_model_values[index_extlaw])*0.1
        
            pn.log_.level = 3 # set this to 3 to have more details
            # Adopt an extinction law
            extinction_law = name_ext
        
            # Define the data file
            obs_data = 'test.dat'
            
            ### Read and deredden observational data
            # define an Observation object and assign it to name 'obs'
            obs = pn.Observation()
            obs_fake = pn.Observation()
                
            # fill obs with data read from file obs_data, with lines varying across rows and a default percent error on line intensities
            obs.readData(obs_data, fileFormat='lines_in_rows_err_cols', errIsRelative=False)#, err_default=0.05)
            obs_fake.readData(obs_data, fileFormat='lines_in_rows_err_cols', errIsRelative=False)#, err_default=0.05)
        
            index_MC=param_mod_name.index("number_fake_obs_MC")
            obs_fake.addMonteCarloObs(N = param_model_values[index_MC])
        
            obs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
            obs_fake.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        
            obs.extinction.law = extinction_law
            obs.extinction.R_V=RV_ext
            obs.correctData(normWave=4861.)
        
            obs_fake.extinction.law = extinction_law
            obs_fake.extinction.R_V=RV_ext
            obs_fake.correctData(normWave=4861.)
            
            cHb.append(0)
            cHb_err.append(0)
            cHb[jrad]=float(obs.extinction.cHbeta[0])
            cHb_err[jrad]=float(np.std(obs_fake.extinction.cHbeta))    
        
            RC = pn.RedCorr(E_BV = obs.extinction.E_BV[0], law = extinction_law)
            corHa=RC.getCorr(6563, 4861) 
            corHa_er=RC.getErrCorr(6563, np.std(obs_fake.extinction.E_BV),4861)
            corHb=RC.getCorr(4861, 4861) 
            corHb_er=RC.getErrCorr(4861, np.std(obs_fake.extinction.E_BV),4861)
            corHg=RC.getCorr(4340, 4861) 
            corHg_er=RC.getErrCorr(4340, np.std(obs_fake.extinction.E_BV),4861)
            corHd=RC.getCorr(4101, 4861) 
            corHd_er=RC.getErrCorr(4101, np.std(obs_fake.extinction.E_BV),4861)
        
            corHeIa=RC.getCorr(5876, 4861) 
            corHeIa_er=RC.getErrCorr(5876, np.std(obs_fake.extinction.E_BV),4861)
            corHeIb=RC.getCorr(6678, 4861) 
            corHeIb_er=RC.getErrCorr(6678, np.std(obs_fake.extinction.E_BV),4861)
        
            corHeIIa=RC.getCorr(4686, 4861) 
            corHeIIa_er=RC.getErrCorr(4686, np.std(obs_fake.extinction.E_BV),4861)
            corHeIIb=RC.getCorr(5412, 4861) 
            corHeIIb_er=RC.getErrCorr(5412, np.std(obs_fake.extinction.E_BV),4861)
        
            corHeIIa=RC.getCorr(4686, 4861) 
            corHeIIa_er=RC.getErrCorr(4686, np.std(obs_fake.extinction.E_BV),4861)
            corHeIIb=RC.getCorr(5412, 4861) 
            corHeIIb_er=RC.getErrCorr(5412, np.std(obs_fake.extinction.E_BV),4861)
            
            corNI=RC.getCorr(5199, 4861) 
            corNI_er=RC.getErrCorr(5199, np.std(obs_fake.extinction.E_BV),4861)
            corNIIa=RC.getCorr(5755, 4861) 
            corNIIa_er=RC.getErrCorr(5755, np.std(obs_fake.extinction.E_BV),4861)
            corNIIb=RC.getCorr(6548, 4861) 
            corNIIb_er=RC.getErrCorr(6548, np.std(obs_fake.extinction.E_BV),4861)
            corNIIc=RC.getCorr(6584, 4861) 
            corNIIc_er=RC.getErrCorr(6584, np.std(obs_fake.extinction.E_BV),4861)
                
            corOIa=RC.getCorr(5577, 4861) 
            corOIa_er=RC.getErrCorr(5577, np.std(obs_fake.extinction.E_BV),4861)
            corOIb=RC.getCorr(6300, 4861) 
            corOIb_er=RC.getErrCorr(6300, np.std(obs_fake.extinction.E_BV),4861)
            corOIc=RC.getCorr(6363, 4861) 
            corOIc_er=RC.getErrCorr(6363, np.std(obs_fake.extinction.E_BV),4861)
                
            corOIIa=RC.getCorr(3727, 4861) 
            corOIIa_er=RC.getErrCorr(3727, np.std(obs_fake.extinction.E_BV),4861)
            corOIIb=RC.getCorr(3729, 4861) 
            corOIIb_er=RC.getErrCorr(3729, np.std(obs_fake.extinction.E_BV),4861)
            corOIIc=RC.getCorr(7320, 4861) 
            corOIIc_er=RC.getErrCorr(7320, np.std(obs_fake.extinction.E_BV),4861)
            corOIId=RC.getCorr(7330, 4861) 
            corOIId_er=RC.getErrCorr(7330, np.std(obs_fake.extinction.E_BV),4861)
        
            corOIIIa=RC.getCorr(4363, 4861)    
            corOIIIa_er=RC.getErrCorr(4363, np.std(obs_fake.extinction.E_BV),4861)  
            corOIIIb=RC.getCorr(4959, 4861)    
            corOIIIb_er=RC.getErrCorr(4959, np.std(obs_fake.extinction.E_BV),4861)  
            corOIIIc=RC.getCorr(5007, 4861)    
            corOIIIc_er=RC.getErrCorr(5007, np.std(obs_fake.extinction.E_BV),4861)  
                
            corSIIa=RC.getCorr(6716, 4861) 
            corSIIa_er=RC.getErrCorr(6716, np.std(obs_fake.extinction.E_BV),4861)
            corSIIb=RC.getCorr(6731, 4861) 
            corSIIb_er=RC.getErrCorr(6731, np.std(obs_fake.extinction.E_BV),4861)
            corSIIIa=RC.getCorr(6312, 4861) 
            corSIIIa_er=RC.getErrCorr(6312, np.std(obs_fake.extinction.E_BV),4861)
            corSIIIb=RC.getCorr(9069, 4861) 
            corSIIIb_er=RC.getErrCorr(9069, np.std(obs_fake.extinction.E_BV),4861)
                
            corClIIIa=RC.getCorr(5517, 4861)    
            corClIIIa_er=RC.getErrCorr(5517, np.std(obs_fake.extinction.E_BV),4861)  
            corClIIIb=RC.getCorr(5538, 4861)    
            corClIIIb_er=RC.getErrCorr(5538, np.std(obs_fake.extinction.E_BV),4861)  
                
            corNeIIIa=RC.getCorr(3868, 4861)    
            corNeIIIa_er=RC.getErrCorr(3868, np.std(obs_fake.extinction.E_BV),4861)  
            corNeIIIb=RC.getCorr(3967, 4861)    
            corNeIIIb_er=RC.getErrCorr(3967, np.std(obs_fake.extinction.E_BV),4861)
                
            corCI=RC.getCorr(8727, 4861)    
            corCI_er=RC.getErrCorr(8727, np.std(obs_fake.extinction.E_BV),4861)  
            corCII=RC.getCorr(6461, 4861)    
            corCII_er=RC.getErrCorr(6461, np.std(obs_fake.extinction.E_BV),4861)
                
            corArIII=RC.getCorr(7136, 4861)    
            corArIII_er=RC.getErrCorr(7136, np.std(obs_fake.extinction.E_BV),4861)  
            corArIVa=RC.getCorr(4712, 4861)    
            corArIVa_er=RC.getErrCorr(4712, np.std(obs_fake.extinction.E_BV),4861)  
            corArIVb=RC.getCorr(4740, 4861)    
            corArIVb_er=RC.getErrCorr(4740, np.std(obs_fake.extinction.E_BV),4861)  
        
        
            
            index_extra_error=line_names.index("HI_6563s")      
            index_extra_error2=line_names.index("HI_4861s")          
            if lines_available[index_extra_error]=="yes" and lines_available[index_extra_error2]=="yes":
                flux_radial_derred.Ha_6563.append(float(obs.getIntens()['H1r_6563A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['H1r_6563A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHa_er/corHa)
                flux_radial_derred_error.Ha_6563.append(flux_radial_derred.Ha_6563[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                flux_radial_derred_error.Hb_4861.append(0.0)
            else:
                print("Halpha and/or Hbeta emission lines are missing!. Program aborted!")
        
            index_extra_error=line_names.index("HI_4340e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.Hg_4340.append(float(obs.getIntens()['H1r_4341A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['H1r_4341A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHg_er/corHg)
                flux_radial_derred_error.Hg_4340.append(flux_radial_derred.Hg_4340[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("HI_4101e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.Hd_4101.append(float(obs.getIntens()['H1r_4102A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['H1r_4102A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHd_er/corHd)
                flux_radial_derred_error.Hd_4101.append(flux_radial_derred.Hd_4101[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("HeII_4686e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.HeIIa_4686.append(float(obs.getIntens()['He2r_4686A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['He2r_4686A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIIa_er/corHeIIa)
                flux_radial_derred_error.HeIIa_4686.append(flux_radial_derred.HeIIa_4686[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("HeII_5412e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.HeIIb_5412.append(float(obs.getIntens()['He2r_5411A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['He2r_5411A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIIb_er/corHeIIb)
                flux_radial_derred_error.HeIIb_5412.append(flux_radial_derred.HeIIb_5412[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("HeI_5876e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.HeIa_5876.append(float(obs.getIntens()['He1r_5876A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['He1r_5876A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIa_er/corHeIa)
                flux_radial_derred_error.HeIa_5876.append(flux_radial_derred.HeIa_5876[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                    
            index_extra_error=line_names.index("HeI_6678e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.HeIb_6678.append(float(obs.getIntens()['He1r_6678A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['He1r_6678A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIb_er/corHeIb)
                flux_radial_derred_error.HeIb_6678.append(flux_radial_derred.HeIb_6678[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("N2_5755e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.NIIa_5755.append(float(obs.getIntens()['N2_5755A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['N2_5755A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIa_er/corNIIa)
                flux_radial_derred_error.NIIa_5755.append(flux_radial_derred.NIIa_5755[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("N2_6548e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.NIIb_6548.append(float(obs.getIntens()['N2_6548A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['N2_6548A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIb_er/corNIIb)
                flux_radial_derred_error.NIIb_6548.append(flux_radial_derred.NIIb_6548[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("N2_6583e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.NIIc_6584.append(float(obs.getIntens()['N2_6584A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['N2_6584A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIc_er/corNIIc)
                flux_radial_derred_error.NIIc_6584.append(flux_radial_derred.NIIc_6584[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                    
            index_extra_error=line_names.index("N1_5199e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.NI_5199.append(float(obs.getIntens()['N1_5200A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['N1_5200A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNI_er/corNI)
                flux_radial_derred_error.NI_5199.append(flux_radial_derred.NI_5199[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("O3_4363e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIIIa_4363.append(float(obs.getIntens()['O3_4363A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O3_4363A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIa_er/corOIIIa)
                flux_radial_derred_error.OIIIa_4363.append(flux_radial_derred.OIIIa_4363[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("O3_4959e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIIIb_4959.append(float(obs.getIntens()['O3_4959A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O3_4959A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIb_er/corOIIIb)
                flux_radial_derred_error.OIIIb_4959.append(flux_radial_derred.OIIIb_4959[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                    
            index_extra_error=line_names.index("O3_5007e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIIIc_5007.append(float(obs.getIntens()['O3_5007A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O3_5007A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIc_er/corOIIIc)
                flux_radial_derred_error.OIIIc_5007.append(flux_radial_derred.OIIIc_5007[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("O2_3727e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIIa_3727.append(float(obs.getIntens()['O2_3726A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O2_3726A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIa_er/corOIIa)
                flux_radial_derred_error.OIIa_3727.append(flux_radial_derred.OIIa_3727[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("O2_3729e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIIb_3729.append(float(obs.getIntens()['O2_3729A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O2_3729A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIb_er/corOIIb)
                flux_radial_derred_error.OIIb_3729.append(flux_radial_derred.OIIb_3729[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("O2_7320e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIIc_7320.append(float(obs.getIntens()['O2_7319A+']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O2_7319A+'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIc_er/corOIIc)
                flux_radial_derred_error.OIIc_7320.append(flux_radial_derred.OIIc_7320[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("O2_7330e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIId_7330.append(float(obs.getIntens()['O2_7330A+']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O2_7330A+'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIId_er/corOIId)
                flux_radial_derred_error.OIId_7330.append(flux_radial_derred.OIId_7330[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("O1_5577e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIa_5577.append(float(obs.getIntens()['O1_5577A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O1_5577A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIa_er/corOIa)
                flux_radial_derred_error.OIa_5577.append(flux_radial_derred.OIa_5577[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("O1_6300e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIb_6300.append(float(obs.getIntens()['O1_6300A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O1_6300A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIb_er/corOIb)
                flux_radial_derred_error.OIb_6300.append(flux_radial_derred.OIb_6300[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("O1_6363e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.OIc_6363.append(float(obs.getIntens()['O1_6364A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['O1_6364A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIc_er/corOIc)
                flux_radial_derred_error.OIc_6363.append(flux_radial_derred.OIc_6363[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                
            index_extra_error=line_names.index("S2_6716e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.SIIa_6716.append(float(obs.getIntens()['S2_6716A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['S2_6716A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIa_er/corSIIa)
                flux_radial_derred_error.SIIa_6716.append(flux_radial_derred.SIIa_6716[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                
            index_extra_error=line_names.index("S2_6731e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.SIIb_6731.append(float(obs.getIntens()['S2_6731A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['S2_6731A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIb_er/corSIIb)
                flux_radial_derred_error.SIIb_6731.append(flux_radial_derred.SIIb_6731[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                    
            index_extra_error=line_names.index("S3_6312e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.SIIIa_6312.append(float(obs.getIntens()['S3_6312A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['S3_6312A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIIa_er/corSIIIa)
                flux_radial_derred_error.SIIIa_6312.append(flux_radial_derred.SIIIa_6312[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("S3_9069e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.SIIIb_9069.append(float(obs.getIntens()['S3_9069A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['S3_9069A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIIb_er/corSIIIb)
                flux_radial_derred_error.SIIIb_9069.append(flux_radial_derred.SIIIb_9069[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                
            index_extra_error=line_names.index("Cl3_5517e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.ClIIIa_5517.append(float(obs.getIntens()['Cl3_5518A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['Cl3_5518A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corClIIIa_er/corClIIIa)
                flux_radial_derred_error.ClIIIa_5517.append(flux_radial_derred.ClIIIa_5517[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("Cl3_5538e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.ClIIIb_5538.append(float(obs.getIntens()['Cl3_5538A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['Cl3_5538A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corClIIIb_er/corClIIIb)
                flux_radial_derred_error.ClIIIb_5538.append(flux_radial_derred.ClIIIb_5538[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                        
            index_extra_error=line_names.index("Ar3_7136e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.ArIII_7136.append(float(obs.getIntens()['Ar3_7136A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['Ar3_7136A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIII_er/corArIII)
                flux_radial_derred_error.ArIII_7136.append(flux_radial_derred.ArIII_7136[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                    
            index_extra_error=line_names.index("Ar4_4712e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.ArIVa_4712.append(float(obs.getIntens()['Ar4_4711A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['Ar4_4711A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIVa_er/corArIVa)
                flux_radial_derred_error.ArIVa_4712.append(flux_radial_derred.ArIVa_4712[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("Ar4_4740e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.ArIVb_4740.append(float(obs.getIntens()['Ar4_4740A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['Ar4_4740A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIVb_er/corArIVb)
                flux_radial_derred_error.ArIVb_4740.append(flux_radial_derred.ArIVb_4740[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("C1_8727e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.CI_8727.append(float(obs.getIntens()['C1_8728A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['C1_8728A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corCI_er/corCI)
                flux_radial_derred_error.CI_8727.append(flux_radial_derred.CI_8727[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("C2_6461e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.CII_6461.append(float(obs.getIntens()['H1r_6563A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['H1r_6563A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corCII_er/corCII)
                flux_radial_derred_error.CII_6461.append(flux_radial_derred.CII_6461[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            index_extra_error=line_names.index("Ne3_3868e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.NeIIIa_3868.append(float(obs.getIntens()['Ne3_3869A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['Ne3_3869A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNeIIIa_er/corNeIIIa)
                flux_radial_derred_error.NeIIIa_3868.append(flux_radial_derred.NeIIIa_3868[jrad]*np.sqrt(par1**2+par2**2+par3**2))
                
            index_extra_error=line_names.index("Ne3_3967e")          
            if lines_available[index_extra_error]=="yes":
                flux_radial_derred.NeIIIb_3967.append(float(obs.getIntens()['Ne3_3968A']))
                flux_radial_derred.Hb_4861.append(float(obs.getIntens()['H1r_4861A']))
                par1=float(obs.getError()['Ne3_3968A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNeIIIb_er/corNeIIIb)
                flux_radial_derred_error.NeIIIb_3967.append(flux_radial_derred.NeIIIb_3967[jrad]*np.sqrt(par1**2+par2**2+par3**2))
        
            jrad=jrad+1
        
    ##############################################################################
    ##############################################################################
    ##############################################################################

    plt.legend(loc='upper right')
    plt.xlabel("Distance from CS (arcsec)",fontsize=21,labelpad=6)
    plt.ylabel("Normalized flux (reddered) ",fontsize=21,labelpad=6)
    ax3.tick_params(axis='x', labelsize=21, pad=6)
    ax3.tick_params(axis='y', labelsize=21, pad=6)
    plt.show()

    
    for i in range(0,len(flux_radial_derred.Ha_6563)):
        index_line=line_names.index("HI_6563s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.Ha_6563[i]=flux_radial_derred.Ha_6563[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("HI_4861s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.Hb_4861[i]=flux_radial_derred.Hb_4861[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("HI_4340s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.Hg_4340[i]=flux_radial_derred.Hg_4340[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("HI_4101s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.Hd_4101[i]=flux_radial_derred.Hd_4101[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("HeI_5876s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.HeIa_5876[i]=flux_radial_derred.HeIa_5876[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("HeI_6678s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.HeIb_6678[i]=flux_radial_derred.HeIb_6678[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("HeII_4686s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.HeIIa_4686[i]=flux_radial_derred.HeIIa_4686[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("HeII_5412s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.HeIIb_5412[i]=flux_radial_derred.HeIIb_5412[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("N1_5199s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.NI_5199[i]=flux_radial_derred.NI_5199[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("N2_5755s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.NIIa_5755[i]=flux_radial_derred.NIIa_5755[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("N2_6548s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.NIIb_6548[i]=flux_radial_derred.NIIb_6548[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("N2_6583s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.NIIc_6584[i]=flux_radial_derred.NIIc_6584[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O1_5577s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIa_5577[i]=flux_radial_derred.OIa_5577[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O1_6300s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIb_6300[i]=flux_radial_derred.OIb_6300[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O1_6363s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIc_6363[i]=flux_radial_derred.OIc_6363[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O2_3727s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIIa_3727[i]=flux_radial_derred.OIIa_3727[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O2_3729s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIIb_3729[i]=flux_radial_derred.OIIb_3729[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O2_7320s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIIc_7320[i]=flux_radial_derred.OIIc_7320[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O2_7330s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIId_7330[i]=flux_radial_derred.OIId_7330[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O3_4363s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIIIa_4363[i]=flux_radial_derred.OIIIa_4363[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O3_4959s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIIIb_4959[i]=flux_radial_derred.OIIIb_4959[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("O3_5007s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.OIIIc_5007[i]=flux_radial_derred.OIIIc_5007[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("S2_6716s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.SIIa_6716[i]=flux_radial_derred.SIIa_6716[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("S2_6731s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.SIIb_6731[i]=flux_radial_derred.SIIb_6731[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("S3_6312s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.SIIIa_6312[i]=flux_radial_derred.SIIIa_6312[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("S3_9069s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.SIIIb_9069[i]=flux_radial_derred.SIIIb_9069[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("Cl3_5517s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.ClIIIa_5517[i]=flux_radial_derred.ClIIIa_5517[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("Cl3_5538s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.ClIIIb_5538[i]=flux_radial_derred.ClIIIb_5538[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("Ne3_3868s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.NeIIIa_3868[i]=flux_radial_derred.NeIIIa_3868[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("Ne3_3967s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.NeIIIb_3967[i]=flux_radial_derred.NeIIIb_3967[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("Ar3_7136s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.ArIII_7136[i]=flux_radial_derred.ArIII_7136[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("Ar4_4712s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.ArIVa_4712[i]=flux_radial_derred.ArIVa_4712[i]*flux_radial.Hb_4861[i]/100.        
        index_line=line_names.index("Ar4_4740s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.ArIVb_4740[i]=flux_radial_derred.ArIVb_4740[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("C1_8727s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.CI_8727[i]=flux_radial_derred.CI_8727[i]*flux_radial.Hb_4861[i]/100.
        index_line=line_names.index("C2_6461s")
        if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
            flux_radial_derred.CII_6461[i]=flux_radial_derred.CII_6461[i]*flux_radial.Hb_4861[i]/100.

    fig,ax4=plt.subplots(figsize=(12,7))

    max_flux,rad_formax=fms.maxvalue(flux_radial_derred.Ha_6563,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
    radial=flux_radial_derred.Ha_6563/max_flux
    plt.plot(radial_dis_arc3,radial,label="HI 6563")

#    max_flux,rad_formax=fms.maxvalue(flux_radial_derred.Hb_4861,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
#    radial=flux_radial_derred.Hb_4861/max_flux
#    plt.plot(radial_dis_arc3,radial,label="HI 4861")
    
    index_line=line_names.index("HI_4340s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.Hg_4340,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.Hg_4340/max_flux
        plt.plot(radial_dis_arc33,radial,label="HI 4340")
    
    index_line=line_names.index("HI_4101s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.Hd_4101,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.Hd_4101/max_flux
        plt.plot(radial_dis_arc33,radial,label="HI 4101")
    
    index_line=line_names.index("HeI_5876s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.HeIa_5876,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.HeIa_5876/max_flux
        plt.plot(radial_dis_arc3,radial,label="He I 5876")
    
    index_line=line_names.index("HeI_6678s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.HeIb_6678,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.HeIb_6678/max_flux
        plt.plot(radial_dis_arc3,radial,label="He I 6678")
    
    index_line=line_names.index("HeII_4686s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.HeIIa_4686,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.HeIIa_4686/max_flux
        plt.plot(radial_dis_arc3,radial,label="He II 4686")
    
    index_line=line_names.index("HeII_5412s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.HeIIb_5412,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.HeIIb_5412/max_flux
        plt.plot(radial_dis_arc3,radial,label="He II 5412")
    
    index_line=line_names.index("N2_5755s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.NIIa_5755,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.NIIa_5755/max_flux
        plt.plot(radial_dis_arc3,radial,label="[N II] 5755")
    
    index_line=line_names.index("N2_6548s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.NIIb_6548,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.NIIb_6548/max_flux
        plt.plot(radial_dis_arc3,radial,label="[N II] 6548")

    index_line=line_names.index("N2_6583s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.NIIc_6584,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.NIIc_6584/max_flux
        plt.plot(radial_dis_arc3,radial,label="[N II] 6584")
        
    index_line=line_names.index("N1_5199s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.NI_5199,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.NI_5199/max_flux
        plt.plot(radial_dis_arc3,radial,label="[N I] 5199")
    
    index_line=line_names.index("O3_4363s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIIIa_4363,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIIIa_4363/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O III] 4363")
    
    index_line=line_names.index("O3_4959s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIIIb_4959,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIIIb_4959/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O III] 4959")
     
    index_line=line_names.index("O3_5007s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":   
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIIIc_5007,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIIIc_5007/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O III] 5007")
    
    index_line=line_names.index("O2_3727s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIIa_3727,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIIa_3727/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O II] 3727")
       
    index_line=line_names.index("O2_3729s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":     
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIIb_3729,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIIb_3729/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O II] 3729")

    index_line=line_names.index("O2_7320s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":    
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIIc_7320,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIIc_7320/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O II] 7320")
    
    index_line=line_names.index("O2_7330s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIId_7330,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIId_7330/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O II] 7330")
    
    index_line=line_names.index("O1_5577s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIa_5577,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIa_5577/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O I] 5577")
    
    index_line=line_names.index("O1_6300s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIb_6300,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIb_6300/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O I] 6300")
    
    index_line=line_names.index("O1_6363s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.OIc_6363,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.OIc_6363/max_flux
        plt.plot(radial_dis_arc3,radial,label="[O I] 6363")
    
    index_line=line_names.index("S2_6716s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.SIIa_6716,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.SIIa_6716/max_flux
        plt.plot(radial_dis_arc3,radial,label="[S II] 6716")
    
    index_line=line_names.index("S2_6731s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.SIIb_6731,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.SIIb_6731/max_flux
        plt.plot(radial_dis_arc3,radial,label="[S II] 6731")
    
    index_line=line_names.index("S3_6312s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.SIIIa_6312,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.SIIIa_6312/max_flux
        plt.plot(radial_dis_arc3,radial,label="[S III] 6312")
   
    index_line=line_names.index("S3_9069s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.SIIIb_9069,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.SIIIb_9069/max_flux
        plt.plot(radial_dis_arc3,radial,label="[S III] 9069")

    index_line=line_names.index("Cl3_5517s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":    
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.ClIIIa_5517,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.ClIIIa_5517/max_flux
        plt.plot(radial_dis_arc3,radial,label="[Cl III] 5517")
    
    index_line=line_names.index("Cl3_5538s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.ClIIIb_5538,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.ClIIIb_5538/max_flux
        plt.plot(radial_dis_arc3,radial,label="[Cl III] 5538")
    
    index_line=line_names.index("Ar3_7136s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.ArIII_7136,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.ArIII_7136/max_flux
        plt.plot(radial_dis_arc3,radial,label="[Ar III] 7136")
    
    index_line=line_names.index("Ar4_4712s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.ArIVa_4712,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.ArIVa_4712/max_flux
        plt.plot(radial_dis_arc3,radial,label="[Ar IV] 4712")
    
    index_line=line_names.index("Ar4_4740s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.ArIVb_4740,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.ArIVb_4740/max_flux
        plt.plot(radial_dis_arc3,radial,label="[Ar IV] 4740")
    
    index_line=line_names.index("Ne3_3868s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.NeIIIa_3868,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.NeIIIa_3868/max_flux
        plt.plot(radial_dis_arc3,radial,label="[Ne III] 3868")
    
    index_line=line_names.index("Ne3_3967s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.NeIIIb_3967,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.NeIIIb_3967/max_flux
        plt.plot(radial_dis_arc3,radial,label="[Ne III] 3967")
    
    index_line=line_names.index("C1_8727s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.CI_8727,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.CI_8727/max_flux
        plt.plot(radial_dis_arc3,radial,label="[C I] 8727")
    
    index_line=line_names.index("C2_6461s")
    if lines_available[index_line]=="yes" and lines_radial[index_line]=="radial_yes":
        max_flux,rad_formax=fms.maxvalue(flux_radial_derred.CII_6461,radial_dis_arc3,param_mod_name,param_model_values,pixscale)
        radial=flux_radial_derred.CII_6461/max_flux
        plt.plot(radial_dis_arc3,radial,label="[C II] 6461")
    
    print >>file22,"dist ","Ha"
    for i in range(0,len(flux_radial_derred.Ha_6563)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.Ha_6563[i])
                    
#    print >>file22,"dist ","Hb"
#    for i in range(0,len(flux_radial_derred.Hb_4861)):
#        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.Hb_4861[i])
            
    print >>file22,"dist ","Hc"
    for i in range(0,len(flux_radial_derred.Hg_4340)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.Hg_4340[i])
    
    print >>file22,"dist ","Hd"
    for i in range(0,len(flux_radial_derred.Hd_4101)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.Hd_4101[i])

    print >>file22,"dist ","HeI 5876"
    for i in range(0,len(flux_radial_derred.HeIa_5876)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.HeIa_5876[i])
            
    print >>file22,"dist ","HeI 6678"
    for i in range(0,len(flux_radial_derred.HeIb_6678)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.HeIb_6678[i])
    
    print >>file22,"dist ","HeII 4686"
    for i in range(0,len(flux_radial_derred.HeIIa_4686)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.HeIIa_4686[i])
        
    print >>file22,"dist ","HeII 5412"
    for i in range(0,len(flux_radial_derred.HeIIb_5412)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.HeIIb_5412[i])
    
    print >>file22,"dist ","[N II] 5755"
    for i in range(0,len(flux_radial_derred.NIIa_5755)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.NIIa_5755[i])

    print >>file22,"dist ","[N II] 6548"
    for i in range(0,len(flux_radial_derred.NIIb_6548)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.NIIb_6548[i])
        
    print >>file22,"dist ","[N II] 6584"
    for i in range(0,len(flux_radial_derred.NIIc_6584)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.NIIc_6584[i])
    
    print >>file22,"dist ","[N I] 5199"
    for i in range(0,len(flux_radial_derred.NI_5199)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.NI_5199[i])
        
    print >>file22,"dist ","[O III] 4363"
    for i in range(0,len(flux_radial_derred.OIIIa_4363)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIIIa_4363[i])
       
    print >>file22,"dist ","[O III] 4959"
    for i in range(0,len(flux_radial_derred.OIIIb_4959)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIIIb_4959[i])
    
    print >>file22,"dist ","[O III] 5007"
    for i in range(0,len(flux_radial_derred.OIIIc_5007)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIIIc_5007[i])
    
    print >>file22,"dist ","[O II] 3727"
    for i in range(0,len(flux_radial_derred.OIIa_3727)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIIa_3727[i])
        
    print >>file22,"dist ","[O II] 3729"
    for i in range(0,len(flux_radial_derred.OIIb_3729)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIIb_3729[i])
        
    print >>file22,"dist ","[O II] 7320"
    for i in range(0,len(flux_radial_derred.OIIc_7320)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIIc_7320[i])
        
    print >>file22,"dist ","[O II] 7330"
    for i in range(0,len(flux_radial_derred.OIId_7330)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIId_7330[i])
        
    print >>file22,"dist ","[O I] 5577"
    for i in range(0,len(flux_radial_derred.OIa_5577)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIa_5577[i])
        
    print >>file22,"dist ","[O I] 6300"
    for i in range(0,len(flux_radial_derred.OIb_6300)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIb_6300[i])
        
    print >>file22,"dist ","[O I] 6363"
    for i in range(0,len(flux_radial_derred.OIc_6363)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.OIc_6363[i])
        
    print >>file22,"dist ","[S II] 6716"
    for i in range(0,len(flux_radial_derred.SIIa_6716)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.SIIa_6716[i])
   
    print >>file22,"dist ","[S II] 6731"
    for i in range(0,len(flux_radial_derred.SIIb_6731)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.SIIb_6731[i])
    
    print >>file22,"dist ","[S III] 6312"
    for i in range(0,len(flux_radial_derred.SIIIa_6312)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.SIIIa_6312[i])
        
    print >>file22,"dist ","[S III] 9069"
    for i in range(0,len(flux_radial_derred.SIIIb_9069)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.SIIIb_9069[i])
        
    print >>file22,"dist ","[Cl III] 5517"
    for i in range(0,len(flux_radial_derred.ClIIIa_5517)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.ClIIIa_5517[i])
        
    print >>file22,"dist ","[Cl III] 5538"
    for i in range(0,len(flux_radial_derred.ClIIIb_5538)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.ClIIIb_5538[i])
        
    print >>file22,"dist ","[Ar III] 7136"
    for i in range(0,len(flux_radial_derred.ArIII_7136)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.ArIII_7136[i])
        
    print >>file22,"dist ","[Ar IV] 4712"
    for i in range(0,len(flux_radial_derred.ArIVa_4712)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.ArIVa_4712[i]) 

    print >>file22,"dist ","[Ar IV] 4740"
    for i in range(0,len(flux_radial_derred.ArIVb_4740)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.ArIVb_4740[i])
    
    print >>file22,"dist ","[Ne III] 3868"
    for i in range(0,len(flux_radial_derred.NeIIIa_3868)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.NeIIIa_3868[i])

    print >>file22,"dist ","[Ne III] 3967"
    for i in range(0,len(flux_radial_derred.NeIIIb_3967)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.NeIIIb_3967[i])
       
    print >>file22,"dist ","[C I] 8727"
    for i in range(0,len(flux_radial_derred.CI_8727)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.CI_8727[i])
            
    print >>file22,"dist ","[C II] 6461"
    for i in range(0,len(flux_radial_derred.CII_6461)):
        print >>file22, '{:03.3f}'.format(radial_dis_arc3[i]),'{:06.3f}'.format(flux_radial_derred.CII_6461[i])
    
    
    plt.legend(loc='upper right')
    plt.xlabel("Distance from CS (arcsec)",fontsize=21,labelpad=6)
    plt.ylabel("Normalized flux (dereddered) ",fontsize=21,labelpad=6)
    ax4.tick_params(axis='x', labelsize=21, pad=6)
    ax4.tick_params(axis='y', labelsize=21, pad=6)
    plt.show()
    
    # plot a specific rage of radius from 0 up to 11 row/column
    cHb2=[]
    cHb2_err=[]
    for i in range(0,len(cHb)):
        cHb2.append(0)
        cHb2_err.append(0)
        radial_dis_arc2.append(0)
        radial_dis_arc2[i]=radial_dis_arc[i]
        cHb2[i]=cHb[i]
        cHb2_err[i]=cHb_err[i]
    
    index_plotname=par_plotname.index("cHbeta")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    fig,ax4=plt.subplots(figsize=(12,7))
    plt.errorbar(radial_dis_arc2,cHb2,yerr=cHb2_err,label=r'c(H$\beta$)')
    plt.xlabel("Distance from CS (arcsec)",fontsize=21,labelpad=6)
    plt.ylabel(r'c(H$\beta$)',fontsize=21,labelpad=6)
    ax4.tick_params(axis='x', labelsize=21, pad=6)
    ax4.tick_params(axis='y', labelsize=21, pad=6)
    plt.ylim(ymin,ymax)
    fig.savefig("output_radial_plots/c(Hb)_radial.png")
    fig.savefig("output_radial_plots/c(Hb)_radial.pdf")


    return flux_radial, flux_radial_error, radial_dis_arc
     
  
