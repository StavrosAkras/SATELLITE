# rotate_line_fluxes_script.py:
# Calculates the emission line fluxes, line intensities (relative to Hbeta), extincion c(Hb) and emission line ratio defined by the user
# in the output.txt file, for a number of radial pseudo-slits with PA form 0 to 360 degrees.
# The extincion coefficient is computed using the PyNeb package.
# The minimum, maximum and step angles are given by the user as well as the width and the length of the slits.
# The results are also saved in two external files, so the user can performe any further analysis. 
# The name of the files are : 'output_linesintensities_per_angles.txt' and 'output_lineratios_per_angles.txt'
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



from satellite import slit_line_flux_script as slfs
from satellite import print_results_script as prs
from satellite import norm_flux_error_script as nfes



def rotate_line_fluxes(flux2D,flux2D_error,line_names,line_ext_error,lines_available,param_estimated,param_requered,param_mod_name,param_model_values,par_plotname,par_plotymin,par_plotymax):
    
    
    class ratio_angles:
        HeIa5876_Ha=[]
        HeIb6678_Ha=[]
        HeIIa4686_Hb=[]
        HeIIb5412_Hb=[]
        HeIa5876_HeIIa4686=[]
        HeIa5876_HeIIb5412=[]
        NIIc6583_Ha=[]
        NIIbc654884_Ha=[]
        NIIbc654884_NIIa5755=[]
        NIIbc654884_OIIIbc_495907=[]
        NI5199_Hb=[]  
        SIIab671631_Ha=[]
        SIIa6716_SIIb6731=[]
        SIIab671631_SIIIab_631269=[] 
        OIbc630063_Ha=[]     
        OIbc630063_OIa5577=[]   
        OIbc630063_OIIIbc495907=[]
        OIbc630063_OIIcd732030=[]
        OIbc630063_OIIab372729=[]
        OIb6300_Ha=[]
        OIIIc5007_Hb=[]
        OIIIbc495907_Hb=[]       
        OIIIbc45907_OIIIa4363=[]   
        OIIab372729_Hb=[]
        OIIab372729_OIIIbc495907=[] 
        OIIcd732030_OIIIbc495907=[]
        ArIVa4712_ArIVb4740=[]
        ArIVab471240_Hb=[]
        NeIIIab386867_Hb=[]   
        ClIIIa5517_ClIIIb5538=[]
        ClIIIab551738_Hb=[]
        CI8727_Ha=[]
        CII6461_Ha=[]
 
    class ratio_angles_error:
        HeIa5876_Ha=[]
        HeIb6678_Ha=[]
        HeIIa4686_Hb=[]
        HeIIb5412_Hb=[]
        HeIa5876_HeIIa4686=[]
        HeIa5876_HeIIb5412=[]
        NIIc6583_Ha=[]
        NIIbc654884_Ha=[]
        NIIbc654884_NIIa5755=[]
        NIIbc654884_OIIIbc_495907=[]
        NI5199_Hb=[]  
        SIIab671631_Ha=[]
        SIIa6716_SIIb6731=[]
        SIIab671631_SIIIab_631269=[] 
        OIbc630063_Ha=[]     
        OIbc630063_OIa5577=[]   
        OIbc630063_OIIIbc495907=[]
        OIbc630063_OIIcd732030=[]
        OIbc630063_OIIab372729=[]
        OIb6300_Ha=[]
        OIIIc5007_Hb=[]
        OIIIbc495907_Hb=[]       
        OIIIbc45907_OIIIa4363=[]   
        OIIab372729_Hb=[]
        OIIab372729_OIIIbc495907=[] 
        OIIcd732030_OIIIbc495907=[]
        ArIVa4712_ArIVb4740=[]
        ArIVab471240_Hb=[]
        NeIIIab386867_Hb=[]   
        ClIIIa5517_ClIIIb5538=[]
        ClIIIab551738_Hb=[]
        CI8727_Ha=[]
        CII6461_Ha=[]


    class flux_angles:
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

    class flux_angles_error:
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



    class flux_angles_norm:
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

    class flux_angles_norm_error:
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

 

    
    index_step_ang=param_mod_name.index("position_angle_step")
    index_star_ang=param_mod_name.index("minimum_position_angle")
    index_end_ang=param_mod_name.index("maximum_position_angle")
    start_angle=param_model_values[index_star_ang]
    end_angle=param_model_values[index_end_ang]
    step_angle=param_model_values[index_step_ang]

    index_tot_pix=param_mod_name.index("total_num_pixels_horiz")
    index_slit_wid=param_mod_name.index("slit_width_in_spaxels")
    index_slit_len=param_mod_name.index("slit_length_in_spaxels")
    slitwidth=param_model_values[index_slit_wid]
    slitlength=param_model_values[index_slit_len]
    maxsize=param_model_values[index_tot_pix]
    
    dataHanew_err,dataHbnew_err,dataHgnew_err,dataHdnew_err,dataHeIanew_err,dataHeIbnew_err,dataHeIIanew_err,dataHeIIbnew_err,dataNIIanew_err,dataNIIbnew_err,dataNIIcnew_err,dataNInew_err,dataOIIIanew_err,dataOIIIbnew_err,dataOIIIcnew_err,dataOIIanew_err,dataOIIbnew_err,dataOIIcnew_err,dataOIIdnew_err,dataOIanew_err,dataOIbnew_err,dataOIcnew_err,dataSIIanew_err,dataSIIbnew_err,dataSIIIanew_err,dataSIIIbnew_err,dataClIIIanew_err,dataClIIIbnew_err,dataArIIInew_err,dataArIVanew_err,dataArIVbnew_err,dataCInew_err,dataCIInew_err,dataNeIIIanew_err,dataNeIIIbnew_err=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    
    dataHanew,dataHbnew,dataHgnew,dataHdnew,dataHeIanew,dataHeIbnew,dataHeIIanew,dataHeIIbnew,dataNIIanew,dataNIIbnew,dataNIIcnew,dataNInew,dataOIIIanew,dataOIIIbnew,dataOIIIcnew,dataOIIanew,dataOIIbnew,dataOIIcnew,dataOIIdnew,dataOIanew,dataOIbnew,dataOIcnew,dataSIIanew,dataSIIbnew,dataSIIIanew,dataSIIIbnew,dataClIIIanew,dataClIIIbnew,dataArIIInew,dataArIVanew,dataArIVbnew,dataCInew,dataCIInew,dataNeIIIanew,dataNeIIIbnew=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    ang=[]
    Ha,Hb,Hg,Hd,HeIa,HeIb,HeIIa,HeIIb,NIIa,NIIb,NIIc,NI,OIIIa,OIIIb,OIIIc,OIIa,OIIb,OIIc,OIId,OIa,OIb,OIc,SIIa,SIIb,SIIIa,SIIIb,ClIIIa,ClIIIb,ArIII,ArIVa,ArIVb,CI,CII,NeIIIa,NeIIIb=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    Ha_err,Hb_err,Hg_err,Hd_err,HeIa_err,HeIb_err,HeIIa_err,HeIIb_err,NIIa_err,NIIb_err,NIIc_err,NI_err,OIIIa_err,OIIIb_err,OIIIc_err,OIIa_err,OIIb_err,OIIc_err,OIId_err,OIa_err,OIb_err,OIc_err,SIIa_err,SIIb_err,SIIIa_err,SIIIb_err,ClIIIa_err,ClIIIb_err,ArIII_err,ArIVa_err,ArIVb_err,CI_err,CII_err,NeIIIa_err,NeIIIb_err=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    cHb,cHg,cHd=[],[],[] 
    cHb_err,cHg_err,cHd_err=[],[],[]

    ratio_HeIa_Ha,ratio_HeIb_Ha,ratio_HeIIa_Hb,ratio_HeIIb_Hb,ratio_HeIa_HeIIa,ratio_HeIa_HeIIb,ratio_NIIc_Ha,ratio_NIIbc_Ha,ratio_NIIbc_NIIa,ratio_NIIbc_OIIIbc,ratio_NI_Hb,ratio_SIIab_Ha,ratio_SIIa_SIIb,ratio_SIIab_SIIIab,ratio_OIbc_Ha,ratio_OIbc_OIa,ratio_OIbc_OIIIbc,ratio_OIbc_OIIcd,ratio_OIbc_OIIab,ratio_OIb_Ha,ratio_OIIIc_Hb,ratio_OIIIbc_Hb,ratio_OIIIbc_OIIIa,ratio_OIIab_Hb,ratio_OIIab_OIIIbc,ratio_OIIcd_OIIIbc,ratio_ArIVa_ArIVb,ratio_ArIVab_Hb,ratio_NeIIIab_Hb,ratio_ClIIIa_ClIIIb,ratio_ClIIIab_Hb,ratio_CI_Ha,ratio_CII_Ha=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    ratio_HeIa_Ha_err,ratio_HeIb_Ha_err,ratio_HeIIa_Hb_err,ratio_HeIIb_Hb_err,ratio_HeIa_HeIIa_err,ratio_HeIa_HeIIb_err,ratio_NIIc_Ha_err,ratio_NIIbc_Ha_err,ratio_NIIbc_NIIa_err,ratio_NIIbc_OIIIbc_err,ratio_NI_Hb_err,ratio_SIIab_Ha_err,ratio_SIIa_SIIb_err,ratio_SIIab_SIIIab_err,ratio_OIbc_Ha_err,ratio_OIbc_OIa_err,ratio_OIbc_OIIIbc_err,ratio_OIbc_OIIcd_err,ratio_OIbc_OIIab_err,ratio_OIb_Ha_err,ratio_OIIIc_Hb_err,ratio_OIIIbc_Hb_err,ratio_OIIIbc_OIIIa_err,ratio_OIIab_Hb_err,ratio_OIIab_OIIIbc_err,ratio_OIIcd_OIIIbc_err,ratio_ArIVa_ArIVb_err,ratio_ArIVab_Hb_err,ratio_NeIIIab_Hb_err,ratio_ClIIIa_ClIIIb_err,ratio_ClIIIab_Hb_err,ratio_CI_Ha_err,ratio_CII_Ha_err=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    file1 = open('output_linesintensities_per_angles.txt','w')
    print >>file1,"ang ","cHb ","OII3727","OII3729","NeIII3868", "NeIII3967","Hd    ", "Hg    ","OIII4363","HeII4686","ArIV4712","ArIV4740","Hb     ","OIII4959","OIII5007","NI5200","HeI5412","ClIII5517","ClIII5538","OI5577","NII5755","HeI5876","OI6300","SIII6312","OI6363","CII6461","NII6548","Ha     ","NII6584","HeI6678","SII6716","SII6731","ArIII7136","OII7320","OII7330","CI8727","SIII9069",'Total f(Hb)'

    file2 = open('output_lineratios_per_angles.txt','w')
    print >>file2,"ang","c_Hb ","c_Hg ","c_Hd ","log(HeI_5876s/HI_6563s)","log(HeI_6678s/HI_6563s)","log(HeII_5412s/HI_4861s)","log(HeII_4686s/HI_4861s)","log(HeI_5876s/HeII_4686s)","log(HeI_5876s/HeII_5412s)","log(N2_6583s/HI_6563s)","log(N2_6548s+N2_6583s)/HI_6563s)","log((N2_6548s+N2_6583s)/N2_5755s)","log((N2_6548s+N2_6583s)/(O3_4959s+O3_5007s))","log(N1_5199s/HI_4861s)","log((S2_6716s+S2_6731s)/HI_6563s)","log(S2_6716s/S2_6731s)","log((S2_6716s+S2_6731s)/(S3_6312s+S3_9069s))","log((O1_6300s+O1_6363s)/HI_6563s)","log((O1_6300s+O1_6363s)/O1_5577s)","log((O1_6300s+O1_6363s)/(O3_4959s+O3_5007s))","log((O1_6300s+O1_6363s)/(O2_7320s+O2_7330s))","log((O1_6300s+O1_6363s)/(O2_3727s+O2_3729s))","log(O1_6300s/HI_6563s)","log(O3_5007s/HI_4861s)","log((O3_4959s+O3_5007s)/HI_4861s)","log((O3_4959s+O3_5007s)/O3_4363s)","log((O2_3727s+O2_3729s)/HI_4861s)","log((O2_3727s+O2_3729s)/(O3_4959s+O3_5007s))","log((O2_7320s+O2_7330s)/(O3_4959s+O3_5007s))","log(Ar4_4712s/Ar4_4740s)","log((Ar4_4712s+Ar4_4740s)/HI_4861s)","log((Ne3_3868s+Ne3_3967s)/HI_4861s)","log(Cl3_5517s/Cl3_5538s)","log((Cl3_5517s+Cl3_5538s)/HI_4861s)","log(C1_8727s/HI_6563s)","log(C2_6461s/HI_6563s)"   

    pn.log_.open_file('my_log.txt')
    
    #In this "for", the fluxes of various emission lines requested by the used  
    # are calculated for a specific PA range, minmum and maximum provided by the used
    # and for a specific angle step. The flxues for all the different PAs are stored in tables.
    kkangle=0
    for i in range(start_angle,end_angle,step_angle): # -90 to 90 degree slit, 0 to 360 for the entire nebula
        dataHanew.append(0),dataHbnew.append(0),dataHgnew.append(0),dataHdnew.append(0),dataHeIanew.append(0),dataHeIbnew.append(0),dataHeIIanew.append(0),dataHeIIbnew.append(0),dataNIIanew.append(0),dataNIIbnew.append(0),dataNIIcnew.append(0),dataNInew.append(0),dataOIIIanew.append(0),dataOIIIbnew.append(0),dataOIIIcnew.append(0),dataOIIanew.append(0),dataOIIbnew.append(0),dataOIIcnew.append(0),dataOIIdnew.append(0),dataOIanew.append(0),dataOIbnew.append(0),dataOIcnew.append(0),dataSIIanew.append(0),dataSIIbnew.append(0),dataSIIIanew.append(0),dataSIIIbnew.append(0),dataClIIIanew.append(0),dataClIIIbnew.append(0),dataArIIInew.append(0),dataArIVanew.append(0),dataArIVbnew.append(0),dataCInew.append(0),dataCIInew.append(0),dataNeIIIanew.append(0),dataNeIIIbnew.append(0)
        ang.append(0)
    
        cHb.append(0),cHg.append(0),cHd.append(0)
        cHb_err.append(0),cHg_err.append(0),cHd_err.append(0)

        Ha.append(0),Hb.append(0),Hg.append(0),Hd.append(0),HeIa.append(0),HeIb.append(0),HeIIa.append(0),HeIIb.append(0),NIIa.append(0),NIIb.append(0),NIIc.append(0),NI.append(0),OIIIa.append(0),OIIIb.append(0),OIIIc.append(0),OIIa.append(0),OIIb.append(0),OIIc.append(0),OIId.append(0),OIa.append(0),OIb.append(0),OIc.append(0),SIIa.append(0),SIIb.append(0),SIIIa.append(0),SIIIb.append(0),ClIIIa.append(0),ClIIIb.append(0),ArIII.append(0),ArIVa.append(0),ArIVb.append(0),CI.append(0),CII.append(0),NeIIIa.append(0),NeIIIb.append(0)
        
        Ha_err.append(0),Hb_err.append(0),Hg_err.append(0),Hd_err.append(0),HeIa_err.append(0),HeIb_err.append(0),HeIIa_err.append(0),HeIIb_err.append(0),NIIa_err.append(0),NIIb_err.append(0),NIIc_err.append(0),NI_err.append(0),OIIIa_err.append(0),OIIIb_err.append(0),OIIIc_err.append(0),OIIa_err.append(0),OIIb_err.append(0),OIIc_err.append(0),OIId_err.append(0),OIa_err.append(0),OIb_err.append(0),OIc_err.append(0),SIIa_err.append(0),SIIb_err.append(0),SIIIa_err.append(0),SIIIb_err.append(0),ClIIIa_err.append(0),ClIIIb_err.append(0),ArIII_err.append(0),ArIVa_err.append(0),ArIVb_err.append(0),CI_err.append(0),CII_err.append(0),NeIIIa_err.append(0),NeIIIb_err.append(0)
   
        ratio_HeIa_Ha.append(-10e-10),ratio_HeIb_Ha.append(-10e-10),ratio_HeIIa_Hb.append(-10e-10),ratio_HeIIb_Hb.append(-10e-10),ratio_HeIa_HeIIa.append(-10e-10),ratio_HeIa_HeIIb.append(-10e-10),ratio_NIIc_Ha.append(-10e-10),ratio_NIIbc_Ha.append(-10e-10),ratio_NIIbc_NIIa.append(-10e-10),ratio_NIIbc_OIIIbc.append(-10e-10),ratio_NI_Hb.append(-10e-10),ratio_SIIab_Ha.append(-10e-10),ratio_SIIa_SIIb.append(-10e-10),ratio_SIIab_SIIIab.append(-10e-10),ratio_OIbc_Ha.append(-10e-10),ratio_OIbc_OIa.append(-10e-10),ratio_OIbc_OIIIbc.append(-10e-10),ratio_OIbc_OIIcd.append(-10e-10),ratio_OIbc_OIIab.append(-10e-10),ratio_OIb_Ha.append(-10e-10),ratio_OIIIc_Hb.append(-10e-10),ratio_OIIIbc_Hb.append(-10e-10),ratio_OIIIbc_OIIIa.append(-10e-10),ratio_OIIab_Hb.append(-10e-10),ratio_OIIab_OIIIbc.append(-10e-10),ratio_OIIcd_OIIIbc.append(-10e-10),ratio_ArIVa_ArIVb.append(-10e-10),ratio_ArIVab_Hb.append(-10e-10),ratio_NeIIIab_Hb.append(-10e-10),ratio_ClIIIa_ClIIIb.append(-10e-10),ratio_ClIIIab_Hb.append(-10e-10),ratio_CI_Ha.append(-10e-10),ratio_CII_Ha.append(-10e-10)

        ratio_HeIa_Ha_err.append(-10e-10),ratio_HeIb_Ha_err.append(-10e-10),ratio_HeIIa_Hb_err.append(-10e-10),ratio_HeIIb_Hb_err.append(-10e-10),ratio_HeIa_HeIIa_err.append(-10e-10),ratio_HeIa_HeIIb_err.append(-10e-10),ratio_NIIc_Ha_err.append(-10e-10),ratio_NIIbc_Ha_err.append(-10e-10),ratio_NIIbc_NIIa_err.append(-10e-10),ratio_NIIbc_OIIIbc_err.append(-10e-10),ratio_NI_Hb_err.append(-10e-10),ratio_SIIab_Ha_err.append(-10e-10),ratio_SIIa_SIIb_err.append(-10e-10),ratio_SIIab_SIIIab_err.append(-10e-10),ratio_OIbc_Ha_err.append(-10e-10),ratio_OIbc_OIa_err.append(-10e-10),ratio_OIbc_OIIIbc_err.append(-10e-10),ratio_OIbc_OIIcd_err.append(-10e-10),ratio_OIbc_OIIab_err.append(-10e-10),ratio_OIb_Ha_err.append(-10e-10),ratio_OIIIc_Hb_err.append(-10e-10),ratio_OIIIbc_Hb_err.append(-10e-10),ratio_OIIIbc_OIIIa_err.append(-10e-10),ratio_OIIab_Hb_err.append(-10e-10),ratio_OIIab_OIIIbc_err.append(-10e-10),ratio_OIIcd_OIIIbc_err.append(-10e-10),ratio_ArIVa_ArIVb_err.append(-10e-10),ratio_ArIVab_Hb_err.append(-10e-10),ratio_NeIIIab_Hb_err.append(-10e-10),ratio_ClIIIa_ClIIIb_err.append(-10e-10),ratio_ClIIIab_Hb_err.append(-10e-10),ratio_CI_Ha_err.append(-10e-10),ratio_CII_Ha_err.append(-10e-10)

        dataHanew_err.append(0),dataHbnew_err.append(0),dataHgnew_err.append(0),dataHdnew_err.append(0),dataHeIanew_err.append(0),dataHeIbnew_err.append(0),dataHeIIanew_err.append(0),dataHeIIbnew_err.append(0),dataNIIanew_err.append(0),dataNIIbnew_err.append(0),dataNIIcnew_err.append(0),dataNInew_err.append(0),dataOIIIanew_err.append(0),dataOIIIbnew_err.append(0),dataOIIIcnew_err.append(0),dataOIIanew_err.append(0),dataOIIbnew_err.append(0),dataOIIcnew_err.append(0),dataOIIdnew_err.append(0),dataOIanew_err.append(0),dataOIbnew_err.append(0),dataOIcnew_err.append(0),dataSIIanew_err.append(0),dataSIIbnew_err.append(0),dataSIIIanew_err.append(0),dataSIIIbnew_err.append(0),dataClIIIanew_err.append(0),dataClIIIbnew_err.append(0),dataArIIInew_err.append(0),dataArIVanew_err.append(0),dataArIVbnew_err.append(0),dataCInew_err.append(0),dataCIInew_err.append(0),dataNeIIIanew_err.append(0),dataNeIIIbnew_err.append(0)

#   dataHanew[] array with the flxues of Ha line for different angles. the "a" parameter 
#   corresponds to the number of pixels used for the flux and angle.  
#   correct fluxes for interstellar extintion and normalized to Hbeta=100 for all the different angles.

        index_extra_error_line=line_names.index("HI_6563s")
        if lines_available[index_extra_error_line]=="yes":
            dataHanew[kkangle],dataHanew_err[kkangle],num_pix_Ha=slfs.slitfluxline(flux2D.Ha_6563,flux2D_error.Ha_6563,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
        else: 
            print("Halpha emission line is missing!. Program aborted!")
        
        index_extra_error_line=line_names.index("HI_4861s")
        if lines_available[index_extra_error_line]=="yes":
            dataHbnew[kkangle],dataHbnew_err[kkangle],num_pix_Hb=slfs.slitfluxline(flux2D.Hb_4861,flux2D_error.Hb_4861,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
        else: 
            print("Hbeta emission line is missing!. Program aborted!")
        
        
        if dataHanew[kkangle]==0 or dataHbnew[kkangle]==0:
            print("Halpha and/or Hbeta emission line is missing!. Program aborted!")
            
        if dataHanew[kkangle]>0 and dataHbnew[kkangle] > 0:
            
            index_extra_error=line_names.index("HI_6563e")          
            if lines_available[index_extra_error]=="yes":
                if line_ext_error[index_extra_error]!=0:
                    errextraHa=float(line_ext_error[index_extra_error])
                    dataHanew_err[kkangle]=dataHanew_err[kkangle]+dataHanew[kkangle]*errextraHa
                else:
                    errextraHa=float(line_ext_error[index_extra_error])
                    dataHanew_err[kkangle]=dataHanew[kkangle]*errextraHa
            
            index_extra_error=line_names.index("HI_4861e")          
            if lines_available[index_extra_error]=="yes":
                if line_ext_error[index_extra_error]!=0:
                    errextraHb=float(line_ext_error[index_extra_error])
                    dataHbnew_err[kkangle]=dataHbnew_err[kkangle]+dataHbnew[kkangle]*errextraHb
                else:
                    errextraHb=float(line_ext_error[index_extra_error])
                    dataHbnew_err[kkangle]=dataHbnew[kkangle]*errextraHb
            
            
            index_line=line_names.index("HI_4340s")
            if lines_available[index_line]=="yes":
                dataHgnew[kkangle],dataHgnew_err[kkangle],num_pix_Hg=slfs.slitfluxline(flux2D.Hg_4340,flux2D_error.Hg_4340,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("HI_4340e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHg=float(line_ext_error[index_extra_error])
                        dataHgnew_err[kkangle]=dataHgnew_err[kkangle]+dataHgnew[kkangle]*errextraHg
                    else:
                        errextraHg=float(line_ext_error[index_extra_error])
                        dataHgnew_err[kkangle]=dataHgnew[kkangle]*errextraHg    
                
            index_line=line_names.index("HI_4101s")
            if lines_available[index_line]=="yes":
                dataHdnew[kkangle],dataHdnew_err[kkangle],num_pix_Hd=slfs.slitfluxline(flux2D.Hd_4101,flux2D_error.Hd_4101,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("HI_4101e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHd=float(line_ext_error[index_extra_error])
                        dataHdnew_err[kkangle]=dataHdnew_err[kkangle]+dataHdnew[kkangle]*errextraHd
                    else:
                        errextraHd=float(line_ext_error[index_extra_error])
                        dataHdnew_err[kkangle]=dataHdnew[kkangle]*errextraHd    
           
            index_line=line_names.index("HeI_5876s")
            if lines_available[index_line]=="yes":
                dataHeIanew[kkangle],dataHeIanew_err[kkangle],num_pix_HeIa=slfs.slitfluxline(flux2D.HeIa_5876,flux2D_error.HeIa_5876,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("HeI_5876e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIa=float(line_ext_error[index_extra_error])
                        dataHeIanew_err[kkangle]=dataHeIanew_err[kkangle]+dataHeIanew[kkangle]*errextraHeIa
                    else:
                        errextraHeIa=float(line_ext_error[index_extra_error])
                        dataHeIanew_err[kkangle]=dataHeIanew[kkangle]*errextraHeIa    
            
            index_line=line_names.index("HeI_6678s")
            if lines_available[index_line]=="yes":
                dataHeIbnew[kkangle],dataHeIbnew_err[kkangle],num_pix_HeIb=slfs.slitfluxline(flux2D.HeIb_6678,flux2D_error.HeIb_6678,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("HeI_6678e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIb=float(line_ext_error[index_extra_error])
                        dataHeIbnew_err[kkangle]=dataHeIbnew_err[kkangle]+dataHeIbnew[kkangle]*errextraHeIb
                    else:
                        errextraHeIb=float(line_ext_error[index_extra_error])
                        dataHeIbnew_err[kkangle]=dataHeIbnew[kkangle]*errextraHeIb  
            
            index_line=line_names.index("HeII_4686s")
            if lines_available[index_line]=="yes":
                dataHeIIanew[kkangle],dataHeIIanew_err[kkangle],num_pix_HeIIa=slfs.slitfluxline(flux2D.HeIIa_4686,flux2D_error.HeIIa_4686,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize) 
                
                index_extra_error=line_names.index("HeII_4686e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIIa=float(line_ext_error[index_extra_error])
                        dataHeIIanew_err[kkangle]=dataHeIIanew_err[kkangle]+dataHeIIanew[kkangle]*errextraHeIIa
                    else:
                        errextraHeIIa=float(line_ext_error[index_extra_error])
                        dataHeIIanew_err[kkangle]=dataHeIIanew[kkangle]*errextraHeIIa 
            
            index_line=line_names.index("HeII_5412s")
            if lines_available[index_line]=="yes":
                dataHeIIbnew[kkangle],dataHeIIbnew_err[kkangle],num_pix_HeIIb=slfs.slitfluxline(flux2D.HeIIb_5412,flux2D_error.HeIIb_5412,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("HeII_5412e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIIb=float(line_ext_error[index_extra_error])
                        dataHeIIbnew_err[kkangle]=dataHeIIbnew_err[kkangle]+dataHeIIbnew[kkangle]*errextraHeIIb
                    else:
                        errextraHeIIb=float(line_ext_error[index_extra_error])
                        dataHeIIbnew_err[kkangle]=dataHeIIbnew[kkangle]*errextraHeIIb
           
            index_line=line_names.index("N2_5755s")
            if lines_available[index_line]=="yes":
                dataNIIanew[kkangle],dataNIIanew_err[kkangle],num_pix_NIIa=slfs.slitfluxline(flux2D.NIIa_5755,flux2D_error.NIIa_5755,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("N2_5755e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIa=float(line_ext_error[index_extra_error])
                        dataNIIanew_err[kkangle]=dataNIIanew_err[kkangle]+dataNIIanew[kkangle]*errextraNIIa
                    else:
                        errextraNIIa=float(line_ext_error[index_extra_error])
                        dataNIIanew_err[kkangle]=dataNIIanew[kkangle]*errextraNIIa
            
            index_line=line_names.index("N2_6548s")
            if lines_available[index_line]=="yes":
                dataNIIbnew[kkangle],dataNIIbnew_err[kkangle],num_pix_NIIb=slfs.slitfluxline(flux2D.NIIb_6548,flux2D_error.NIIb_6548,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("N2_6548e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIb=float(line_ext_error[index_extra_error])
                        dataNIIbnew_err[kkangle]=dataNIIbnew_err[kkangle]+dataNIIbnew[kkangle]*errextraNIIb
                    else:
                        errextraNIIb=float(line_ext_error[index_extra_error])
                        dataNIIbnew_err[kkangle]=dataNIIbnew[kkangle]*errextraNIIb 

            index_line=line_names.index("N2_6583s")
            if lines_available[index_line]=="yes":
                dataNIIcnew[kkangle],dataNIIcnew_err[kkangle],num_pix_NIIc=slfs.slitfluxline(flux2D.NIIc_6584,flux2D_error.NIIc_6584,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("N2_6583e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIc=float(line_ext_error[index_extra_error])
                        dataNIIcnew_err[kkangle]=dataNIIcnew_err[kkangle]+dataNIIcnew[kkangle]*errextraNIIc
                    else:
                        errextraNIIc=float(line_ext_error[index_extra_error])
                        dataNIIcnew_err[kkangle]=dataNIIcnew[kkangle]*errextraNIIc 
            
            index_line=line_names.index("N1_5199s")
            if lines_available[index_line]=="yes":
                dataNInew[kkangle],dataNInew_err[kkangle],num_pix_NI=slfs.slitfluxline(flux2D.NI_5199,flux2D_error.NI_5199,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("N1_5199e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNI=float(line_ext_error[index_extra_error])
                        dataNInew_err[kkangle]=dataNInew_err[kkangle]+dataNInew[kkangle]*errextraNI
                    else:
                        errextraNI=float(line_ext_error[index_extra_error])
                        dataNInew_err[kkangle]=dataNInew[kkangle]*errextraNI
            
           
            index_line=line_names.index("O3_4363s")
            if lines_available[index_line]=="yes":
                dataOIIIanew[kkangle],dataOIIIanew_err[kkangle],num_pix_OIIIa=slfs.slitfluxline(flux2D.OIIIa_4363,flux2D_error.OIIIa_4363,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("O3_4363e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIa=float(line_ext_error[index_extra_error])
                        dataOIIIanew_err[kkangle]=dataOIIIanew_err[kkangle]+dataOIIIanew[kkangle]*errextraOIIIa
                    else:
                        errextraOIIIa=float(line_ext_error[index_extra_error])
                        dataOIIIanew_err[kkangle]=dataOIIIanew[kkangle]*errextraOIIIa
    
            index_line=line_names.index("O3_4959s")
            if lines_available[index_line]=="yes":
                dataOIIIbnew[kkangle],dataOIIIbnew_err[kkangle],num_pix_OIIIb=slfs.slitfluxline(flux2D.OIIIb_4959,flux2D_error.OIIIb_4959,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("O3_4959e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIb=float(line_ext_error[index_extra_error])
                        dataOIIIbnew_err[kkangle]=dataOIIIbnew_err[kkangle]+dataOIIIbnew[kkangle]*errextraOIIIb
                    else:
                        errextraOIIIb=float(line_ext_error[index_extra_error])
                        dataOIIIbnew_err[kkangle]=dataOIIIbnew[kkangle]*errextraOIIIb

            index_line=line_names.index("O3_5007s")
            if lines_available[index_line]=="yes":
                dataOIIIcnew[kkangle],dataOIIIcnew_err[kkangle],num_pix_OIIIc=slfs.slitfluxline(flux2D.OIIIc_5007,flux2D_error.OIIIc_5007,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("O3_5007e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIc=float(line_ext_error[index_extra_error])
                        dataOIIIcnew_err[kkangle]=dataOIIIcnew_err[kkangle]+dataOIIIcnew[kkangle]*errextraOIIIc
                    else:
                        errextraOIIIc=float(line_ext_error[index_extra_error])
                        dataOIIIcnew_err[kkangle]=dataOIIIcnew[kkangle]*errextraOIIIc
            
            index_line=line_names.index("O2_3727s")
            if lines_available[index_line]=="yes":
                dataOIIanew[kkangle],dataOIIanew_err[kkangle],num_pix_OIIa=slfs.slitfluxline(flux2D.OIIa_3727,flux2D_error.OIIa_3727,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("O2_3727e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIa=float(line_ext_error[index_extra_error])
                        dataOIIanew_err[kkangle]=dataOIIanew_err[kkangle]+dataOIIanew[kkangle]*errextraOIIa
                    else:
                        errextraOIIa=float(line_ext_error[index_extra_error])
                        dataOIIanew_err[kkangle]=dataOIIanew[kkangle]*errextraOIIa
            
            index_line=line_names.index("O2_3729s")
            if lines_available[index_line]=="yes":
                dataOIIbnew[kkangle],dataOIIbnew_err[kkangle],num_pix_OIIb=slfs.slitfluxline(flux2D.OIIb_3729,flux2D_error.OIIb_3729,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("O2_3729e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIb=float(line_ext_error[index_extra_error])
                        dataOIIbnew_err[kkangle]=dataOIIbnew_err[kkangle]+dataOIIbnew[kkangle]*errextraOIIb
                    else:
                        errextraOIIb=float(line_ext_error[index_extra_error])
                        dataOIIbnew_err[kkangle]=dataOIIbnew[kkangle]*errextraOIIb
           
            index_line=line_names.index("O2_7320s")
            if lines_available[index_line]=="yes":
                dataOIIcnew[kkangle],dataOIIcnew_err[kkangle],num_pix_OIIc=slfs.slitfluxline(flux2D.OIIc_7320,flux2D_error.OIIc_7320,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
            
                index_extra_error=line_names.index("O2_7320e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIc=float(line_ext_error[index_extra_error])
                        dataOIIcnew_err[kkangle]=dataOIIcnew_err[kkangle]+dataOIIcnew[kkangle]*errextraOIIc
                    else:
                        errextraOIIc=float(line_ext_error[index_extra_error])
                        dataOIIcnew_err[kkangle]=dataOIIcnew[kkangle]*errextraOIIc
                        
            index_line=line_names.index("O2_7330s")
            if lines_available[index_line]=="yes":
                dataOIIdnew[kkangle],dataOIIdnew_err[kkangle],num_pix_OIId=slfs.slitfluxline(flux2D.OIId_7330,flux2D_error.OIId_7330,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("O2_7330e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIId=float(line_ext_error[index_extra_error])
                        dataOIIdnew_err[kkangle]=dataOIIdnew_err[kkangle]+dataOIIdnew[kkangle]*errextraOIId
                    else:
                        errextraOIId=float(line_ext_error[index_extra_error])
                        dataOIIdnew_err[kkangle]=dataOIIdnew[kkangle]*errextraOIId

            index_line=line_names.index("O1_5577s")
            if lines_available[index_line]=="yes":
                dataOIanew[kkangle],dataOIanew_err[kkangle],num_pix_OIa=slfs.slitfluxline(flux2D.OIa_5577,flux2D_error.OIa_5577,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("O1_5577e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIa=float(line_ext_error[index_extra_error])
                        dataOIanew_err[kkangle]=dataOIanew_err[kkangle]+dataOIanew[kkangle]*errextraOIa
                    else:
                        errextraOIa=float(line_ext_error[index_extra_error])
                        dataOIanew_err[kkangle]=dataOIanew[kkangle]*errextraOIa
                
            index_line=line_names.index("O1_6300s")
            if lines_available[index_line]=="yes":
                dataOIbnew[kkangle],dataOIbnew_err[kkangle],num_pix_OIb=slfs.slitfluxline(flux2D.OIb_6300,flux2D_error.OIb_6300,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("O1_6300e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIb=float(line_ext_error[index_extra_error])
                        dataOIbnew_err[kkangle]=dataOIbnew_err[kkangle]+dataOIbnew[kkangle]*errextraOIb
                    else:
                        errextraOIb=float(line_ext_error[index_extra_error])
                        dataOIbnew_err[kkangle]=dataOIbnew[kkangle]*errextraOIb
                    
            index_line=line_names.index("O1_6363s")
            if lines_available[index_line]=="yes":
                dataOIcnew[kkangle],dataOIcnew_err[kkangle],num_pix_OIc=slfs.slitfluxline(flux2D.OIc_6363,flux2D_error.OIc_6363,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
            
                index_extra_error=line_names.index("O1_6363e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIc=float(line_ext_error[index_extra_error])
                        dataOIcnew_err[kkangle]=dataOIcnew_err[kkangle]+dataOIcnew[kkangle]*errextraOIc
                    else:
                        errextraOIc=float(line_ext_error[index_extra_error])
                        dataOIcnew_err[kkangle]=dataOIcnew[kkangle]*errextraOIc
            
            index_line=line_names.index("S2_6716s")
            if lines_available[index_line]=="yes":
                dataSIIanew[kkangle],dataSIIanew_err[kkangle],num_pix_SIIa=slfs.slitfluxline(flux2D.SIIa_6716,flux2D_error.SIIa_6716,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("S2_6716e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIa=float(line_ext_error[index_extra_error])
                        dataSIIanew_err[kkangle]=dataSIIanew_err[kkangle]+dataSIIanew[kkangle]*errextraSIIa
                    else:
                        errextraSIIa=float(line_ext_error[index_extra_error])
                        dataSIIanew_err[kkangle]=dataSIIanew[kkangle]*errextraSIIa
        
            index_line=line_names.index("S2_6731s")
            if lines_available[index_line]=="yes":
                dataSIIbnew[kkangle],dataSIIbnew_err[kkangle],num_pix_SIIb=slfs.slitfluxline(flux2D.SIIb_6731,flux2D_error.SIIb_6731,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("S2_6731e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIb=float(line_ext_error[index_extra_error])
                        dataSIIbnew_err[kkangle]=dataSIIbnew_err[kkangle]+dataSIIbnew[kkangle]*errextraSIIb
                    else:
                        errextraSIIb=float(line_ext_error[index_extra_error])
                        dataSIIbnew_err[kkangle]=dataSIIbnew[kkangle]*errextraSIIb
            
            index_line=line_names.index("S3_6312s")
            if lines_available[index_line]=="yes":
                dataSIIIanew[kkangle],dataSIIIanew_err[kkangle],num_pix_SIIIa=slfs.slitfluxline(flux2D.SIIIa_6312,flux2D_error.SIIIa_6312,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
            
                index_extra_error=line_names.index("S3_6312e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIIa=float(line_ext_error[index_extra_error])
                        dataSIIIanew_err[kkangle]=dataSIIIanew_err[kkangle]+dataSIIIanew[kkangle]*errextraSIIIa
                    else:
                        errextraSIIIa=float(line_ext_error[index_extra_error])
                        dataSIIIanew_err[kkangle]=dataSIIIanew[kkangle]*errextraSIIIa
                        
            index_line=line_names.index("S3_9069s")
            if lines_available[index_line]=="yes":
                dataSIIIbnew[kkangle],dataSIIIbnew_err[kkangle],num_pix_SIIIb=slfs.slitfluxline(flux2D.SIIIb_9069,flux2D_error.SIIIb_9069,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("S3_9069e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIIb=float(line_ext_error[index_extra_error])
                        dataSIIIbnew_err[kkangle]=dataSIIIbnew_err[kkangle]+dataSIIIbnew[kkangle]*errextraSIIIb
                    else:
                        errextraSIIIb=float(line_ext_error[index_extra_error])
                        dataSIIIbnew_err[kkangle]=dataSIIIbnew[kkangle]*errextraSIIIb
                        
            index_line=line_names.index("Cl3_5517s")
            if lines_available[index_line]=="yes":
                dataClIIIanew[kkangle],dataClIIIanew_err[kkangle],num_pix_ClIIIa=slfs.slitfluxline(flux2D.ClIIIa_5517,flux2D_error.ClIIIa_5517,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("Cl3_5517e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraClIIIa=float(line_ext_error[index_extra_error])
                        dataClIIIanew_err[kkangle]=dataClIIIanew_err[kkangle]+dataClIIIanew[kkangle]*errextraClIIIa
                    else:
                        errextraClIIIa=float(line_ext_error[index_extra_error])
                        dataClIIIanew_err[kkangle]=dataClIIIanew[kkangle]*errextraClIIIa

            index_line=line_names.index("Cl3_5538s")
            if lines_available[index_line]=="yes":
                dataClIIIbnew[kkangle],dataClIIIbnew_err[kkangle],num_pix_ClIIIb=slfs.slitfluxline(flux2D.ClIIIb_5538,flux2D_error.ClIIIb_5538,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)

                index_extra_error=line_names.index("Cl3_5538e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraClIIIb=float(line_ext_error[index_extra_error])
                        dataClIIIbnew_err[kkangle]=dataClIIIbnew_err[kkangle]+dataClIIIbnew[kkangle]*errextraClIIIb
                    else:
                        errextraClIIIb=float(line_ext_error[index_extra_error])
                        dataClIIIbnew_err[kkangle]=dataClIIIbnew[kkangle]*errextraClIIIb

            index_line=line_names.index("Ar3_7136s")
            if lines_available[index_line]=="yes":
                dataArIIInew[kkangle],dataArIIInew_err[kkangle],num_pix_ArIII=slfs.slitfluxline(flux2D.ArIII_7136,flux2D_error.ArIII_7136,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
            
                index_extra_error=line_names.index("Ar3_7136e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIII=float(line_ext_error[index_extra_error])
                        dataArIIInew_err[kkangle]=dataArIIInew_err[kkangle]+dataArIIInew[kkangle]*errextraArIII
                    else:
                        errextraArIII=float(line_ext_error[index_extra_error])
                        dataArIIInew_err[kkangle]=dataArIIInew[kkangle]*errextraArIII
            
            index_line=line_names.index("Ar4_4712s")
            if lines_available[index_line]=="yes":
                dataArIVanew[kkangle],dataArIVanew_err[kkangle],num_pix_ArIVa=slfs.slitfluxline(flux2D.ArIVa_4712,flux2D_error.ArIVa_4712,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
            
                index_extra_error=line_names.index("Ar4_4712e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIVa=float(line_ext_error[index_extra_error])
                        dataArIVanew_err[kkangle]=dataArIVanew_err[kkangle]+dataArIVanew[kkangle]*errextraArIVa
                    else:
                        errextraArIVa=float(line_ext_error[index_extra_error])
                        dataArIVanew_err[kkangle]=dataArIVanew[kkangle]*errextraArIVa
                        
            index_line=line_names.index("Ar4_4740s")
            if lines_available[index_line]=="yes":
                dataArIVbnew[kkangle],dataArIVbnew_err[kkangle],num_pix_ArIVb=slfs.slitfluxline(flux2D.ArIVb_4740,flux2D_error.ArIVb_4740,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("Ar4_4740e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIVb=float(line_ext_error[index_extra_error])
                        dataArIVbnew_err[kkangle]=dataArIVbnew_err[kkangle]+dataArIVbnew[kkangle]*errextraArIVb
                    else:
                        errextraArIVb=float(line_ext_error[index_extra_error])
                        dataArIVbnew_err[kkangle]=dataArIVbnew[kkangle]*errextraArIVb
                
            index_line=line_names.index("C1_8727s")
            if lines_available[index_line]=="yes":
                dataCInew[kkangle],dataCInew_err[kkangle],num_pix_CI=slfs.slitfluxline(flux2D.CI_8727,flux2D_error.CI_8727,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
           
                index_extra_error=line_names.index("C1_8727e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraCI=float(line_ext_error[index_extra_error])
                        dataCInew_err[kkangle]=dataCInew_err[kkangle]+dataCInew[kkangle]*errextraCI
                    else:
                        errextraCI=float(line_ext_error[index_extra_error])
                        dataCInew_err[kkangle]=dataCInew[kkangle]*errextraCI
           
            index_line=line_names.index("C2_6461s")
            if lines_available[index_line]=="yes":
                dataCIInew[kkangle],dataCIInew_err[kkangle],num_pix_CII=slfs.slitfluxline(flux2D.CII_6461,flux2D_error.CII_6461,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
                
                index_extra_error=line_names.index("C2_6461e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraCII=float(line_ext_error[index_extra_error])
                        dataCIInew_err[kkangle]=dataCIInew_err[kkangle]+dataCIInew[kkangle]*errextraCII
                    else:
                        errextraCII=float(line_ext_error[index_extra_error])
                        dataCIInew_err[kkangle]=dataCIInew[kkangle]*errextraCII
                
            index_line=line_names.index("Ne3_3868s")
            if lines_available[index_line]=="yes":
                dataNeIIIanew[kkangle],dataNeIIIanew_err[kkangle],num_pix_NeIIIa=slfs.slitfluxline(flux2D.NeIIIa_3868,flux2D_error.NeIIIa_3868,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
    
                index_extra_error=line_names.index("Ne3_3868e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNeIIIa=float(line_ext_error[index_extra_error])
                        dataNeIIIanew_err[kkangle]=dataNeIIIanew_err[kkangle]+dataNeIIIanew[kkangle]*errextraNeIIIa
                    else:
                        errextraNeIIIa=float(line_ext_error[index_extra_error])
                        dataNeIIIanew_err[kkangle]=dataNeIIIanew[kkangle]*errextraNeIIIa
    
            index_line=line_names.index("Ne3_3967s")
            if lines_available[index_line]=="yes":
                dataNeIIIbnew[kkangle],dataNeIIIbnew_err[kkangle],num_pix_NeIIIb=slfs.slitfluxline(flux2D.NeIIIb_3967,flux2D_error.NeIIIb_3967,i,flux2D.Ha_6563,flux2D.Hb_4861,slitwidth,slitlength,maxsize)
            
                index_extra_error=line_names.index("Ne3_3967e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNeIIIb=float(line_ext_error[index_extra_error])
                        dataNeIIIbnew_err[kkangle]=dataNeIIIbnew_err[kkangle]+dataNeIIIbnew[kkangle]*errextraNeIIIb
                    else:
                        errextraNeIIIb=float(line_ext_error[index_extra_error])
                        dataNeIIIbnew_err[kkangle]=dataNeIIIbnew[kkangle]*errextraNeIIIb
            
            file3 = open('test.dat','w')
            print >>file3, "LINE test err"

            error=nfes.flux_Error(dataHanew[kkangle],dataHanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "H1r_6563A", '{:03.3f}'.format(dataHanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
            
            error=nfes.flux_Error(dataHbnew[kkangle],dataHbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "H1r_4861A", '{:03.3f}'.format(dataHbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
                
            error=nfes.flux_Error(dataHgnew[kkangle],dataHgnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "H1r_4341A", '{:03.3f}'.format(dataHgnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHdnew[kkangle],dataHdnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "H1r_4102A", '{:03.3f}'.format(dataHdnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHeIIanew[kkangle],dataHeIIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "He2r_4686A", '{:03.3f}'.format(dataHeIIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHeIIbnew[kkangle],dataHeIIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "He2r_5411A", '{:03.3f}'.format(dataHeIIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHeIanew[kkangle],dataHeIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "He1r_5876A", '{:03.3f}'.format(dataHeIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHeIbnew[kkangle],dataHeIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "He1r_6678A", '{:03.3f}'.format(dataHeIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataSIIanew[kkangle],dataSIIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "S2_6716A", '{:03.3f}'.format(dataSIIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataSIIbnew[kkangle],dataSIIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "S2_6731A", '{:03.3f}'.format(dataSIIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataSIIIanew[kkangle],dataSIIIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "S3_6312A", '{:03.3f}'.format(dataSIIIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataSIIIbnew[kkangle],dataSIIIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "S3_9069A", '{:03.3f}'.format(dataSIIIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNInew[kkangle],dataNInew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "N1_5200A", '{:03.3f}'.format(dataNInew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNIIanew[kkangle],dataNIIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "N2_5755A", '{:03.3f}'.format(dataNIIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNIIbnew[kkangle],dataNIIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "N2_6548A", '{:03.3f}'.format(dataNIIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNIIcnew[kkangle],dataNIIcnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "N2_6584A", '{:03.3f}'.format(dataNIIcnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIanew[kkangle],dataOIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O1_5577A", '{:03.3f}'.format(dataOIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIbnew[kkangle],dataOIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O1_6300A", '{:03.3f}'.format(dataOIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIcnew[kkangle],dataOIcnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O1_6364A", '{:03.3f}'.format(dataOIcnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIanew[kkangle],dataOIIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O2_3726A", '{:03.3f}'.format(dataOIIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIbnew[kkangle],dataOIIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O2_3729A", '{:03.3f}'.format(dataOIIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIcnew[kkangle],dataOIIcnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O2_7319A+", '{:03.3f}'.format(dataOIIcnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIdnew[kkangle],dataOIIdnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O2_7330A+", '{:03.3f}'.format(dataOIIdnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIIanew[kkangle],dataOIIIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O3_4363A", '{:03.3f}'.format(dataOIIIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIIbnew[kkangle],dataOIIIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O3_4959A", '{:03.3f}'.format(dataOIIIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIIcnew[kkangle],dataOIIIcnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "O3_5007A", '{:03.3f}'.format(dataOIIIcnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataClIIIanew[kkangle],dataClIIIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])        
            print >>file3, "Cl3_5518A", '{:03.3f}'.format(dataClIIIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)

            error=nfes.flux_Error(dataClIIIbnew[kkangle],dataClIIIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "Cl3_5538A", '{:03.3f}'.format(dataClIIIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)        
            
            error=nfes.flux_Error(dataNeIIIanew[kkangle],dataNeIIIanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "Ne3_3869A", '{:03.3f}'.format(dataNeIIIanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNeIIIbnew[kkangle],dataNeIIIbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "Ne3_3968A", '{:03.3f}'.format(dataNeIIIbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)                
        
            error=nfes.flux_Error(dataArIIInew[kkangle],dataArIIInew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "Ar3_7136A", '{:03.3f}'.format(dataArIIInew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataArIVanew[kkangle],dataArIVanew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "Ar4_4711A", '{:03.3f}'.format(dataArIVanew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataArIVbnew[kkangle],dataArIVbnew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "Ar4_4740A", '{:03.3f}'.format(dataArIVbnew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataCInew[kkangle],dataCInew_err[kkangle],dataHbnew[kkangle],dataHbnew_err[kkangle])
            print >>file3, "C1_8728A", '{:03.3f}'.format(dataCInew[kkangle]*100/dataHbnew[kkangle]),'{:03.3f}'.format(error)

            ang[kkangle]=i

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
           #if name_ext=='F99_like': name_ext=str('F99') + str('-like')
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
            
            cHb[kkangle]=float(obs.extinction.cHbeta[0])
            cHb_err[kkangle]=float(np.std(obs_fake.extinction.cHbeta))
            
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
            
            Ha[kkangle]=float(obs.getIntens()['H1r_6563A'])
            Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
            par1=float(obs.getError()['H1r_6563A'])
            par2=float(obs.getError()['H1r_4861A'])
            par3=float(corHa_er/corHa)
            Ha_err[kkangle]=Ha[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
            Hb_err[kkangle]=0.0

            index_extra_error=line_names.index("HI_4340e")          
            if lines_available[index_extra_error]=="yes":
                Hg[kkangle]=float(obs.getIntens()['H1r_4341A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['H1r_4341A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHg_er/corHg)
                Hg_err[kkangle]=Hg[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("HI_4101e")          
            if lines_available[index_extra_error]=="yes":
                Hd[kkangle]=float(obs.getIntens()['H1r_4102A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['H1r_4102A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHd_er/corHd)
                Hd_err[kkangle]=Hd[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("HeII_4686e")          
            if lines_available[index_extra_error]=="yes":
                HeIIa[kkangle]=float(obs.getIntens()['He2r_4686A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['He2r_4686A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIIa_er/corHeIIa)
                HeIIa_err[kkangle]=HeIIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("HeII_5412e")          
            if lines_available[index_extra_error]=="yes":
                HeIIb[kkangle]=float(obs.getIntens()['He2r_5411A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['He2r_5411A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIIb_er/corHeIIb)
                HeIIb_err[kkangle]=HeIIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("HeI_5876e")          
            if lines_available[index_extra_error]=="yes":
                HeIa[kkangle]=float(obs.getIntens()['He1r_5876A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['He1r_5876A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIa_er/corHeIa)
                HeIa_err[kkangle]=HeIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("HeI_6678e")          
            if lines_available[index_extra_error]=="yes":
                HeIb[kkangle]=float(obs.getIntens()['He1r_6678A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['He1r_6678A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIb_er/corHeIb)
                HeIb_err[kkangle]=HeIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("N2_5755e")          
            if lines_available[index_extra_error]=="yes":
                NIIa[kkangle]=float(obs.getIntens()['N2_5755A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['N2_5755A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIa_er/corNIIa)
                NIIa_err[kkangle]=NIIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("N2_6548e")          
            if lines_available[index_extra_error]=="yes":
                NIIb[kkangle]=float(obs.getIntens()['N2_6548A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['N2_6548A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIb_er/corNIIb)
                NIIb_err[kkangle]=NIIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("N2_6583e")          
            if lines_available[index_extra_error]=="yes":
                NIIc[kkangle]=float(obs.getIntens()['N2_6584A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['N2_6584A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIc_er/corNIIc)
                NIIc_err[kkangle]=NIIc[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                
            index_extra_error=line_names.index("N1_5199e")          
            if lines_available[index_extra_error]=="yes":
                NI[kkangle]=float(obs.getIntens()['N1_5200A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['N1_5200A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNI_er/corNI)
                NI_err[kkangle]=NI[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("O3_4363e")          
            if lines_available[index_extra_error]=="yes":
                OIIIa[kkangle]=float(obs.getIntens()['O3_4363A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O3_4363A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIa_er/corOIIIa)
                OIIIa_err[kkangle]=OIIIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("O3_4959e")          
            if lines_available[index_extra_error]=="yes":
                OIIIb[kkangle]=float(obs.getIntens()['O3_4959A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O3_4959A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIb_er/corOIIIb)
                OIIIb_err[kkangle]=OIIIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("O3_5007e")          
            if lines_available[index_extra_error]=="yes":
                OIIIc[kkangle]=float(obs.getIntens()['O3_5007A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O3_5007A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIc_er/corOIIIc)
                OIIIc_err[kkangle]=OIIIc[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
     
            index_extra_error=line_names.index("O2_3727e")          
            if lines_available[index_extra_error]=="yes":
                OIIa[kkangle]=float(obs.getIntens()['O2_3726A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O2_3726A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIa_er/corOIIa)
                OIIa_err[kkangle]=OIIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("O2_3729e")          
            if lines_available[index_extra_error]=="yes":
                OIIb[kkangle]=float(obs.getIntens()['O2_3729A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O2_3729A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIb_er/corOIIb)
                OIIb_err[kkangle]=OIIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("O2_7320e")          
            if lines_available[index_extra_error]=="yes":
                OIIc[kkangle]=float(obs.getIntens()['O2_7319A+'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O2_7319A+'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIc_er/corOIIc)
                OIIc_err[kkangle]=OIIc[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
     
            index_extra_error=line_names.index("O2_7330e")          
            if lines_available[index_extra_error]=="yes":
                OIId[kkangle]=float(obs.getIntens()['O2_7330A+'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O2_7330A+'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIId_er/corOIId)
                OIId_err[kkangle]=OIId[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
     
            index_extra_error=line_names.index("O1_5577e")          
            if lines_available[index_extra_error]=="yes":
                OIa[kkangle]=float(obs.getIntens()['O1_5577A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O1_5577A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIa_er/corOIa)
                OIa_err[kkangle]=OIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("O1_6300e")          
            if lines_available[index_extra_error]=="yes":
                OIb[kkangle]=float(obs.getIntens()['O1_6300A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O1_6300A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIb_er/corOIb)
                OIb_err[kkangle]=OIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("O1_6363e")          
            if lines_available[index_extra_error]=="yes":
                OIc[kkangle]=float(obs.getIntens()['O1_6364A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O1_6364A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIc_er/corOIc)
                OIc_err[kkangle]=OIc[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                
            index_extra_error=line_names.index("S2_6716e")          
            if lines_available[index_extra_error]=="yes":
                SIIa[kkangle]=float(obs.getIntens()['S2_6716A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['S2_6716A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIa_er/corSIIa)
                SIIa_err[kkangle]=SIIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
            
            index_extra_error=line_names.index("S2_6731e")          
            if lines_available[index_extra_error]=="yes":
                SIIb[kkangle]=float(obs.getIntens()['S2_6731A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['S2_6731A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIb_er/corSIIb)
                SIIb_err[kkangle]=SIIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("S3_6312e")          
            if lines_available[index_extra_error]=="yes":
                SIIIa[kkangle]=float(obs.getIntens()['S3_6312A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['S3_6312A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIIa_er/corSIIIa)
                SIIIa_err[kkangle]=SIIIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("S3_9069e")          
            if lines_available[index_extra_error]=="yes":
                SIIIb[kkangle]=float(obs.getIntens()['S3_9069A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['S3_9069A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIIb_er/corSIIIb)
                SIIIb_err[kkangle]=SIIIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
            
            index_extra_error=line_names.index("Cl3_5517e")          
            if lines_available[index_extra_error]=="yes":
                ClIIIa[kkangle]=float(obs.getIntens()['Cl3_5518A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Cl3_5518A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corClIIIa_er/corClIIIa)
                ClIIIa_err[kkangle]=ClIIIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("Cl3_5538e")          
            if lines_available[index_extra_error]=="yes":
                ClIIIb[kkangle]=float(obs.getIntens()['Cl3_5538A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Cl3_5538A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corClIIIb_er/corClIIIb)
                ClIIIb_err[kkangle]=ClIIIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("Ar3_7136e")          
            if lines_available[index_extra_error]=="yes":
                ArIII[kkangle]=float(obs.getIntens()['Ar3_7136A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ar3_7136A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIII_er/corArIII)
                ArIII_err[kkangle]=ArIII[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("Ar4_4712e")          
            if lines_available[index_extra_error]=="yes":
                ArIVa[kkangle]=float(obs.getIntens()['Ar4_4711A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ar4_4711A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIVa_er/corArIVa)
                ArIVa_err[kkangle]=ArIVa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("Ar4_4740e")          
            if lines_available[index_extra_error]=="yes":
                ArIVb[kkangle]=float(obs.getIntens()['Ar4_4740A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ar4_4740A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIVb_er/corArIVb)
                ArIVb_err[kkangle]=ArIVb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("C1_8727e")          
            if lines_available[index_extra_error]=="yes":
                CI[kkangle]=float(obs.getIntens()['C1_8728A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['C1_8728A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corCI_er/corCI)
                CI_err[kkangle]=CI[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("C2_6461e")          
            if lines_available[index_extra_error]=="yes":
               CII[kkangle]=float(obs.getIntens()['H1r_6563A'])
               Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
               par1=float(obs.getError()['H1r_6563A'])
               par2=float(obs.getError()['H1r_4861A'])
               par3=float(corCII_er/corCII)
               CII_err[kkangle]=CII[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("Ne3_3868e")          
            if lines_available[index_extra_error]=="yes":
                NeIIIa[kkangle]=float(obs.getIntens()['Ne3_3869A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ne3_3869A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNeIIIa_er/corNeIIIa)
                NeIIIa_err[kkangle]=NeIIIa[kkangle]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("Ne3_3967e")          
            if lines_available[index_extra_error]=="yes":
                NeIIIb[kkangle]=float(obs.getIntens()['Ne3_3968A'])
                Hb[kkangle]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ne3_3968A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNeIIIb_er/corNeIIIb)
                NeIIIb_err[kkangle]=NeIIIb[kkangle]*np.sqrt(par1**2+par2**2+par3**2)
            
            
            flux_angles.Ha_6563.append(dataHanew[kkangle])
            flux_angles.Hb_4861.append(dataHbnew[kkangle])
            flux_angles.Hg_4340.append(dataHgnew[kkangle])
            flux_angles.Hd_4101.append(dataHdnew[kkangle])
            flux_angles.HeIa_5876.append(dataHeIanew[kkangle])
            flux_angles.HeIb_6678.append(dataHeIbnew[kkangle])
            flux_angles.HeIIa_4686.append(dataHeIIanew[kkangle])
            flux_angles.HeIIb_5412.append(dataHeIIbnew[kkangle])
            flux_angles.NIIa_5755.append(dataNIIanew[kkangle])
            flux_angles.NIIb_6548.append(dataNIIbnew[kkangle])
            flux_angles.NIIc_6584.append(dataNIIcnew[kkangle])
            flux_angles.NI_5199.append(dataNInew[kkangle])
            flux_angles.OIIIa_4363.append(dataOIIIanew[kkangle])
            flux_angles.OIIIb_4959.append(dataOIIIbnew[kkangle])
            flux_angles.OIIIc_5007.append(dataOIIIcnew[kkangle])
            flux_angles.OIIa_3727.append(dataOIIanew[kkangle])
            flux_angles.OIIb_3729.append(dataOIIbnew[kkangle])
            flux_angles.OIIc_7320.append(dataOIIcnew[kkangle])
            flux_angles.OIId_7330.append(dataOIIdnew[kkangle])
            flux_angles.OIa_5577.append(dataOIanew[kkangle])
            flux_angles.OIb_6300.append(dataOIbnew[kkangle])
            flux_angles.OIc_6363.append(dataOIcnew[kkangle])
            flux_angles.SIIa_6716.append(dataSIIanew[kkangle])
            flux_angles.SIIb_6731.append(dataSIIbnew[kkangle])
            flux_angles.SIIIa_6312.append(dataSIIIanew[kkangle])
            flux_angles.SIIIb_9069.append(dataSIIIbnew[kkangle])
            flux_angles.ClIIIa_5517.append(dataClIIIanew[kkangle])
            flux_angles.ClIIIb_5538.append(dataClIIIbnew[kkangle])
            flux_angles.ArIII_7136.append(dataArIIInew[kkangle])
            flux_angles.ArIVa_4712.append(dataArIVanew[kkangle])
            flux_angles.ArIVb_4740.append(dataArIVbnew[kkangle])
            flux_angles.CI_8727.append(dataCInew[kkangle])
            flux_angles.CII_6461.append(dataCIInew[kkangle])
            flux_angles.NeIIIa_3868.append(dataNeIIIanew[kkangle])
            flux_angles.NeIIIb_3967.append(dataNeIIIbnew[kkangle])

            flux_angles_error.Ha_6563.append(dataHanew_err[kkangle])
            flux_angles_error.Hb_4861.append(dataHbnew_err[kkangle])
            flux_angles_error.Hg_4340.append(dataHgnew_err[kkangle])
            flux_angles_error.Hd_4101.append(dataHdnew_err[kkangle])
            flux_angles_error.HeIa_5876.append(dataHeIanew_err[kkangle])
            flux_angles_error.HeIb_6678.append(dataHeIbnew_err[kkangle])
            flux_angles_error.HeIIa_4686.append(dataHeIIanew_err[kkangle])
            flux_angles_error.HeIIb_5412.append(dataHeIIbnew_err[kkangle])
            flux_angles_error.NIIa_5755.append(dataNIIanew_err[kkangle])
            flux_angles_error.NIIb_6548.append(dataNIIbnew_err[kkangle])
            flux_angles_error.NIIc_6584.append(dataNIIcnew_err[kkangle])
            flux_angles_error.NI_5199.append(dataNInew_err[kkangle])
            flux_angles_error.OIIIa_4363.append(dataOIIIanew_err[kkangle])
            flux_angles_error.OIIIb_4959.append(dataOIIIbnew_err[kkangle])
            flux_angles_error.OIIIc_5007.append(dataOIIIcnew_err[kkangle])
            flux_angles_error.OIIa_3727.append(dataOIIanew_err[kkangle])
            flux_angles_error.OIIb_3729.append(dataOIIbnew_err[kkangle])
            flux_angles_error.OIIc_7320.append(dataOIIcnew_err[kkangle])
            flux_angles_error.OIId_7330.append(dataOIIdnew_err[kkangle])
            flux_angles_error.OIa_5577.append(dataOIanew_err[kkangle])
            flux_angles_error.OIb_6300.append(dataOIbnew_err[kkangle])
            flux_angles_error.OIc_6363.append(dataOIcnew_err[kkangle])
            flux_angles_error.SIIa_6716.append(dataSIIanew_err[kkangle])
            flux_angles_error.SIIb_6731.append(dataSIIbnew_err[kkangle])
            flux_angles_error.SIIIa_6312.append(dataSIIIanew_err[kkangle])
            flux_angles_error.SIIIb_9069.append(dataSIIIbnew_err[kkangle])
            flux_angles_error.ClIIIa_5517.append(dataClIIIanew_err[kkangle])
            flux_angles_error.ClIIIb_5538.append(dataClIIIbnew_err[kkangle])
            flux_angles_error.ArIII_7136.append(dataArIIInew_err[kkangle])
            flux_angles_error.ArIVa_4712.append(dataArIVanew_err[kkangle])
            flux_angles_error.ArIVb_4740.append(dataArIVbnew_err[kkangle])
            flux_angles_error.CI_8727.append(dataCInew_err[kkangle])
            flux_angles_error.CII_6461.append(dataCIInew_err[kkangle])
            flux_angles_error.NeIIIa_3868.append(dataNeIIIanew_err[kkangle])
            flux_angles_error.NeIIIb_3967.append(dataNeIIIbnew_err[kkangle])
    
            flux_angles_norm.Ha_6563.append(Ha[kkangle])
            flux_angles_norm.Hb_4861.append(Hb[kkangle])
            flux_angles_norm.Hg_4340.append(Hg[kkangle])
            flux_angles_norm.Hd_4101.append(Hd[kkangle])
            flux_angles_norm.HeIa_5876.append(HeIa[kkangle])
            flux_angles_norm.HeIb_6678.append(HeIb[kkangle])
            flux_angles_norm.HeIIa_4686.append(HeIIa[kkangle])
            flux_angles_norm.HeIIb_5412.append(HeIIb[kkangle])
            flux_angles_norm.NIIa_5755.append(NIIa[kkangle])
            flux_angles_norm.NIIb_6548.append(NIIb[kkangle])
            flux_angles_norm.NIIc_6584.append(NIIc[kkangle])
            flux_angles_norm.NI_5199.append(NI[kkangle])
            flux_angles_norm.OIIIa_4363.append(OIIIa[kkangle])
            flux_angles_norm.OIIIb_4959.append(OIIIb[kkangle])
            flux_angles_norm.OIIIc_5007.append(OIIIc[kkangle])
            flux_angles_norm.OIIa_3727.append(OIIa[kkangle])
            flux_angles_norm.OIIb_3729.append(OIIb[kkangle])
            flux_angles_norm.OIIc_7320.append(OIIc[kkangle])
            flux_angles_norm.OIId_7330.append(OIId[kkangle])
            flux_angles_norm.OIa_5577.append(OIa[kkangle])
            flux_angles_norm.OIb_6300.append(OIb[kkangle])
            flux_angles_norm.OIc_6363.append(OIc[kkangle])
            flux_angles_norm.SIIa_6716.append(SIIa[kkangle])
            flux_angles_norm.SIIb_6731.append(SIIb[kkangle])
            flux_angles_norm.SIIIa_6312.append(SIIIa[kkangle])
            flux_angles_norm.SIIIb_9069.append(SIIIb[kkangle])
            flux_angles_norm.ClIIIa_5517.append(ClIIIa[kkangle])
            flux_angles_norm.ClIIIb_5538.append(ClIIIb[kkangle])
            flux_angles_norm.ArIII_7136.append(ArIII[kkangle])
            flux_angles_norm.ArIVa_4712.append(ArIVa[kkangle])
            flux_angles_norm.ArIVb_4740.append(ArIVb[kkangle])
            flux_angles_norm.CI_8727.append(CI[kkangle])
            flux_angles_norm.CII_6461.append(CII[kkangle])
            flux_angles_norm.NeIIIa_3868.append(NeIIIa[kkangle])
            flux_angles_norm.NeIIIb_3967.append(NeIIIb[kkangle])
        
            flux_angles_norm_error.Ha_6563.append(Ha_err[kkangle])
            flux_angles_norm_error.Hb_4861.append(Hb_err[kkangle])
            flux_angles_norm_error.Hg_4340.append(Hg_err[kkangle])
            flux_angles_norm_error.Hd_4101.append(Hd_err[kkangle])
            flux_angles_norm_error.HeIa_5876.append(HeIa_err[kkangle])
            flux_angles_norm_error.HeIb_6678.append(HeIb_err[kkangle])
            flux_angles_norm_error.HeIIa_4686.append(HeIIa_err[kkangle])
            flux_angles_norm_error.HeIIb_5412.append(HeIIb_err[kkangle])
            flux_angles_norm_error.NIIa_5755.append(NIIa_err[kkangle])
            flux_angles_norm_error.NIIb_6548.append(NIIb_err[kkangle])
            flux_angles_norm_error.NIIc_6584.append(NIIc_err[kkangle])
            flux_angles_norm_error.NI_5199.append(NI_err[kkangle])
            flux_angles_norm_error.OIIIa_4363.append(OIIIa_err[kkangle])
            flux_angles_norm_error.OIIIb_4959.append(OIIIb_err[kkangle])
            flux_angles_norm_error.OIIIc_5007.append(OIIIc_err[kkangle])
            flux_angles_norm_error.OIIa_3727.append(OIIa_err[kkangle])
            flux_angles_norm_error.OIIb_3729.append(OIIb_err[kkangle])
            flux_angles_norm_error.OIIc_7320.append(OIIc_err[kkangle])
            flux_angles_norm_error.OIId_7330.append(OIId_err[kkangle])
            flux_angles_norm_error.OIa_5577.append(OIa_err[kkangle])
            flux_angles_norm_error.OIb_6300.append(OIb_err[kkangle])
            flux_angles_norm_error.OIc_6363.append(OIc_err[kkangle])
            flux_angles_norm_error.SIIa_6716.append(SIIa_err[kkangle])
            flux_angles_norm_error.SIIb_6731.append(SIIb_err[kkangle])
            flux_angles_norm_error.SIIIa_6312.append(SIIIa_err[kkangle])
            flux_angles_norm_error.SIIIb_9069.append(SIIIb_err[kkangle])
            flux_angles_norm_error.ClIIIa_5517.append(ClIIIa_err[kkangle])
            flux_angles_norm_error.ClIIIb_5538.append(ClIIIb_err[kkangle])
            flux_angles_norm_error.ArIII_7136.append(ArIII_err[kkangle])
            flux_angles_norm_error.ArIVa_4712.append(ArIVa_err[kkangle])
            flux_angles_norm_error.ArIVb_4740.append(ArIVb_err[kkangle])
            flux_angles_norm_error.CI_8727.append(CI_err[kkangle])
            flux_angles_norm_error.CII_6461.append(CII_err[kkangle])
            flux_angles_norm_error.NeIIIa_3868.append(NeIIIa_err[kkangle])
            flux_angles_norm_error.NeIIIb_3967.append(NeIIIb_err[kkangle])
    
            index_par=param_estimated.index("log(HeI_5876s_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_HeIa_Ha[kkangle]=np.log10(HeIa[kkangle]/Ha[kkangle])
                par1= (HeIa_err[kkangle]/(HeIa[kkangle]*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_HeIa_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log(HeI_6678s_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_HeIb_Ha[kkangle]=np.log10(HeIb[kkangle]/Ha[kkangle])
                par1= (HeIb_err[kkangle]/(HeIb[kkangle]*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_HeIb_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log(HeII_4686s_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_HeIIa_Hb[kkangle]=np.log10(HeIIa[kkangle]/Hb[kkangle])
                par1= (HeIIa_err[kkangle]/(HeIIa[kkangle]*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_HeIIa_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log(HeII_5412s_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_HeIIb_Hb[kkangle]=np.log10(HeIIb[kkangle]/Hb[kkangle])    
                par1= (HeIIb_err[kkangle]/(HeIIb[kkangle]*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_HeIIb_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log(HeI_5876s_HeII_4686s)")
            if param_requered[index_par]=="yes":
                ratio_HeIa_HeIIa[kkangle]=np.log10(HeIa[kkangle]/HeIIa[kkangle])    
                par1= (HeIa_err[kkangle]/(HeIa[kkangle]*np.log(10)))
                par2= (HeIIa_err[kkangle]/(HeIIa[kkangle]*np.log(10)))    
                ratio_HeIa_HeIIa_err[kkangle]=np.sqrt(par1**2+par2**2)
    
            index_par=param_estimated.index("log(HeI_5876s_HeII_5412s)")
            if param_requered[index_par]=="yes":
                ratio_HeIa_HeIIb[kkangle]=np.log10(HeIa[kkangle]/HeIIb[kkangle])   
                par1= (HeIa_err[kkangle]/(HeIa[kkangle]*np.log(10)))
                par2= (HeIIb_err[kkangle]/(HeIIb[kkangle]*np.log(10)))    
                ratio_HeIa_HeIIb_err[kkangle]=np.sqrt(par1**2+par2**2)
    
            index_par=param_estimated.index("log(N2_6583s_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_NIIc_Ha[kkangle]=np.log10(NIIc[kkangle]/Ha[kkangle])   
                par1= (NIIc_err[kkangle]/(NIIc[kkangle]*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_NIIc_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)    
            
            index_par=param_estimated.index("log(N2_6548s+N2_6583s)_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_NIIbc_Ha[kkangle]=np.log10((NIIb[kkangle]+NIIc[kkangle])/Ha[kkangle])
                par1= ((NIIb_err[kkangle]+NIIc_err[kkangle])/((NIIb[kkangle]+NIIc[kkangle])*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_NIIbc_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)
        
            index_par=param_estimated.index("log((N2_6548s+N2_6583s)_N2_5755s)")
            if param_requered[index_par]=="yes":
                ratio_NIIbc_NIIa[kkangle]=np.log10((NIIb[kkangle]+NIIc[kkangle])/NIIa[kkangle])    
                par1= ((NIIb_err[kkangle]+NIIc_err[kkangle])/((NIIb[kkangle]+NIIc[kkangle])*np.log(10)))
                par2= (NIIa_err[kkangle]/(NIIa[kkangle]*np.log(10)))    
                ratio_NIIbc_NIIa_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((N2_6548s+N2_6583s)_(O3_4959s+O3_5007s))")
            if param_requered[index_par]=="yes":
                ratio_NIIbc_OIIIbc[kkangle]=np.log10((NIIb[kkangle]+NIIc[kkangle])/(OIIIb[kkangle]+OIIIc[kkangle]))     
                par1= ((NIIb_err[kkangle]+NIIc_err[kkangle])/((NIIb[kkangle]+NIIc[kkangle])*np.log(10)))
                par2= ((OIIIb_err[kkangle]+OIIIc_err[kkangle])/((OIIIb[kkangle]+OIIIc[kkangle])*np.log(10)))    
                ratio_NIIbc_OIIIbc_err[kkangle]=np.sqrt(par1**2+par2**2)
    
            index_par=param_estimated.index("log(N1_5199s_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_NI_Hb[kkangle]=np.log10(NI[kkangle]/Hb[kkangle])  
                par1= (NI_err[kkangle]/(NI[kkangle]*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_NI_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)  
        
            index_par=param_estimated.index("log((S2_6716s+S2_6731s)_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_SIIab_Ha[kkangle]=np.log10((SIIa[kkangle]+SIIb[kkangle])/Ha[kkangle]) 
                par1= ((SIIa_err[kkangle]+SIIb_err[kkangle])/((SIIa[kkangle]+SIIb[kkangle])*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_SIIab_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log(S2_6716s_S2_6731s)")
            if param_requered[index_par]=="yes":
                ratio_SIIa_SIIb[kkangle]=np.log10(SIIa[kkangle]/SIIb[kkangle]) 
                par1= (SIIa_err[kkangle]/(SIIa[kkangle]*np.log(10)))
                par2= (SIIb_err[kkangle]/(SIIb[kkangle]*np.log(10)))    
                ratio_SIIa_SIIb_err[kkangle]=np.sqrt(par1**2+par2**2)  
            
            index_par=param_estimated.index("log((S2_6716s+S2_6731s)_(S3_6312s+S3_9069s))")
            if param_requered[index_par]=="yes":
                ratio_SIIab_SIIIab[kkangle]=np.log10((SIIa[kkangle]+SIIb[kkangle])/(SIIIa[kkangle]+SIIIb[kkangle]))   
                par1= ((SIIa_err[kkangle]+SIIb_err[kkangle])/((SIIa[kkangle]+SIIb[kkangle])*np.log(10)))
                par2= ((SIIIa_err[kkangle]+SIIIb_err[kkangle])/((SIIIa[kkangle]+SIIIb[kkangle])*np.log(10)))    
                ratio_SIIab_SIIIab_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((O1_6300s+O1_6363s)_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_OIbc_Ha[kkangle]=np.log10((OIb[kkangle]+OIc[kkangle])/Ha[kkangle])     
                par1= ((OIb_err[kkangle]+OIc_err[kkangle])/((OIb[kkangle]+OIc[kkangle])*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_OIbc_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((O1_6300s+O1_6363s)_O1_5577s)")
            if param_requered[index_par]=="yes":
                ratio_OIbc_OIa[kkangle]=np.log10((OIb[kkangle]+OIc[kkangle])/OIa[kkangle])     
                par1= ((OIb_err[kkangle]+OIc_err[kkangle])/((OIb[kkangle]+OIc[kkangle])*np.log(10)))
                par2= (OIa_err[kkangle]/(OIa[kkangle]*np.log(10)))    
                ratio_OIbc_OIa_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((O1_6300s+O1_6363s)_(O3_4959s+O3_5007s))")
            if param_requered[index_par]=="yes":
                ratio_OIbc_OIIIbc[kkangle]=np.log10((OIb[kkangle]+OIc[kkangle])/(OIIIb[kkangle]+OIIIc[kkangle]))  
                par1= ((OIb_err[kkangle]+OIc_err[kkangle])/((OIb[kkangle]+OIc[kkangle])*np.log(10)))
                par2= ((OIIIb_err[kkangle]+OIIIc_err[kkangle])/((OIIIb[kkangle]+OIIIc[kkangle])*np.log(10)))    
                ratio_OIbc_OIIIbc_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((O1_6300s+O1_6363s)_(O2_7320s+O2_7330s))")
            if param_requered[index_par]=="yes":
                ratio_OIbc_OIIcd[kkangle]=np.log10((OIb[kkangle]+OIc[kkangle])/(OIIc[kkangle]+OIId[kkangle]))
                par1= ((OIb_err[kkangle]+OIc_err[kkangle])/((OIb[kkangle]+OIc[kkangle])*np.log(10)))
                par2= ((OIIc_err[kkangle]+OIId_err[kkangle])/((OIIc[kkangle]+OIId[kkangle])*np.log(10)))    
                ratio_OIbc_OIIcd_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((O1_6300s+O1_6363s)_(O2_3727s+O2_3729s))")
            if param_requered[index_par]=="yes":
                ratio_OIbc_OIIab[kkangle]=np.log10((OIb[kkangle]+OIc[kkangle])/(OIIa[kkangle]+OIIb[kkangle]))
                par1= ((OIb_err[kkangle]+OIc_err[kkangle])/((OIb[kkangle]+OIc[kkangle])*np.log(10)))
                par2= ((OIIa_err[kkangle]+OIIb_err[kkangle])/((OIIa[kkangle]+OIIb[kkangle])*np.log(10)))    
                ratio_OIbc_OIIab_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log(O1_6300s_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_OIb_Ha[kkangle]=np.log10(OIb[kkangle]/Ha[kkangle])
                par1= (OIb_err[kkangle]/(OIb[kkangle]*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_OIb_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)      
  
            index_par=param_estimated.index("log(O3_5007s_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_OIIIc_Hb[kkangle]=np.log10(OIIIc[kkangle]/Hb[kkangle])
                par1= (OIIIc_err[kkangle]/(OIIIc[kkangle]*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_OIIIc_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)      
  
            index_par=param_estimated.index("log((O3_4959s+O3_5007s)_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_OIIIbc_Hb[kkangle]=np.log10((OIIIb[kkangle]+OIIIc[kkangle])/Hb[kkangle])    
                par1= ((OIIIb_err[kkangle]+OIIIc_err[kkangle])/((OIIIb[kkangle]+OIIIc[kkangle])*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_OIIIbc_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((O3_4959s+O3_5007s)_O3_4363s)")
            if param_requered[index_par]=="yes":
                ratio_OIIIbc_OIIIa[kkangle]=np.log10((OIIIb[kkangle]+OIIIc[kkangle])/OIIIa[kkangle])     
                par1= ((OIIIb_err[kkangle]+OIIIc_err[kkangle])/((OIIIb[kkangle]+OIIIc[kkangle])*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_OIIIbc_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)
        
            index_par=param_estimated.index("log((O2_3727s+O2_3729s)_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_OIIab_Hb[kkangle]=np.log10((OIIa[kkangle]+OIIb[kkangle])/Hb[kkangle])       
                par1= ((OIIa_err[kkangle]+OIIb_err[kkangle])/((OIIa[kkangle]+OIIb[kkangle])*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_OIIab_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((O2_3727s+O2_3729s)_(O3_4959s+O3_5007s))")
            if param_requered[index_par]=="yes":
                ratio_OIIab_OIIIbc[kkangle]=np.log10((OIIa[kkangle]+OIIb[kkangle])/(OIIIb[kkangle]+OIIIc[kkangle]))  
                par1= ((OIIa_err[kkangle]+OIIb_err[kkangle])/((OIIa[kkangle]+OIIb[kkangle])*np.log(10)))
                par2= ((OIIIb_err[kkangle]+OIIIc_err[kkangle])/((OIIIb[kkangle]+OIIIc[kkangle])*np.log(10)))    
                ratio_OIIab_OIIIbc_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((O2_7320s+O2_7330s)_(O3_4959s+O3_5007s))")
            if param_requered[index_par]=="yes":
                ratio_OIIcd_OIIIbc[kkangle]=np.log10((OIIc[kkangle]+OIId[kkangle])/(OIIIb[kkangle]+OIIIc[kkangle]))  
                par1= ((OIIc_err[kkangle]+OIId_err[kkangle])/((OIIc[kkangle]+OIId[kkangle])*np.log(10)))
                par2= ((OIIIb_err[kkangle]+OIIIc_err[kkangle])/((OIIIb[kkangle]+OIIIc[kkangle])*np.log(10)))    
                ratio_OIIcd_OIIIbc_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log(Ar4_4712s_Ar4_4740s)")
            if param_requered[index_par]=="yes":
                ratio_ArIVa_ArIVb[kkangle]=np.log10(ArIVa[kkangle]/ArIVb[kkangle])  
                par1= (ArIVa_err[kkangle]/(ArIVa[kkangle]*np.log(10)))
                par2= (ArIVb_err[kkangle]/(ArIVb[kkangle]*np.log(10)))    
                ratio_ArIVa_ArIVb_err[kkangle]=np.sqrt(par1**2+par2**2)      
        
            index_par=param_estimated.index("log((Ar4_4712s+Ar4_4740s)_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_ArIVab_Hb[kkangle]=np.log10((ArIVa[kkangle]+ArIVb[kkangle])/Hb[kkangle])     
                par1= ((ArIVa_err[kkangle]+ArIVb_err[kkangle])/((ArIVa[kkangle]+ArIVb[kkangle])*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_ArIVab_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)
            
            index_par=param_estimated.index("log((Ne3_3868s+Ne3_3967s)_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_NeIIIab_Hb[kkangle]=np.log10((NeIIIa[kkangle]+NeIIIb[kkangle])/Hb[kkangle])     
                par1= ((NeIIIa_err[kkangle]+NeIIIb_err[kkangle])/((NeIIIa[kkangle]+NeIIIb[kkangle])*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_NeIIIab_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)
    
            index_par=param_estimated.index("log(Cl3_5517s_Cl3_5538s)")
            if param_requered[index_par]=="yes":
                ratio_ClIIIa_ClIIIb[kkangle]=np.log10(ClIIIa[kkangle]/ClIIIb[kkangle]) 
                par1= (ClIIIa_err[kkangle]/(ClIIIa[kkangle]*np.log(10)))
                par2= (ClIIIb_err[kkangle]/(ClIIIb[kkangle]*np.log(10)))    
                ratio_ClIIIa_ClIIIb_err[kkangle]=np.sqrt(par1**2+par2**2)      
  
            index_par=param_estimated.index("log((Cl3_5517s+Cl3_5538s)_HI_4861s)")
            if param_requered[index_par]=="yes":
                ratio_ClIIIab_Hb[kkangle]=np.log10((ClIIIa[kkangle]+ClIIIb[kkangle])/Hb[kkangle])  
                par1= ((ClIIIa_err[kkangle]+ClIIIb_err[kkangle])/((ClIIIa[kkangle]+ClIIIb[kkangle])*np.log(10)))
                par2= (Hb_err[kkangle]/(Hb[kkangle]*np.log(10)))    
                ratio_ClIIIab_Hb_err[kkangle]=np.sqrt(par1**2+par2**2)    
            
            index_par=param_estimated.index("log(C1_8727s_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_CI_Ha[kkangle]=np.log10(CI[kkangle]/Ha[kkangle]) 
                par1= (CI_err[kkangle]/(CI[kkangle]*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_CI_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)      

            index_par=param_estimated.index("log(C2_6461s_HI_6563s)")
            if param_requered[index_par]=="yes":
                ratio_CII_Ha[kkangle]=np.log10(CII[kkangle]/Ha[kkangle]) 
                par1= (CII_err[kkangle]/(CII[kkangle]*np.log(10)))
                par2= (Ha_err[kkangle]/(Ha[kkangle]*np.log(10)))    
                ratio_CII_Ha_err[kkangle]=np.sqrt(par1**2+par2**2)
        
            ratio_angles.HeIa5876_Ha.append(ratio_HeIa_Ha[kkangle])
            ratio_angles.HeIb6678_Ha.append(ratio_HeIb_Ha[kkangle])
            ratio_angles.HeIIa4686_Hb.append(ratio_HeIIa_Hb[kkangle])
            ratio_angles.HeIIb5412_Hb.append(ratio_HeIIb_Hb[kkangle])
            ratio_angles.HeIa5876_HeIIa4686.append(ratio_HeIa_HeIIa[kkangle])
            ratio_angles.HeIa5876_HeIIb5412.append(ratio_HeIa_HeIIb[kkangle])
            ratio_angles.NIIc6583_Ha.append(ratio_NIIc_Ha[kkangle])
            ratio_angles.NIIbc654884_Ha.append(ratio_NIIbc_Ha[kkangle])
            ratio_angles.NIIbc654884_NIIa5755.append(ratio_NIIbc_NIIa[kkangle])
            ratio_angles.NIIbc654884_OIIIbc_495907.append(ratio_NIIbc_OIIIbc[kkangle])
            ratio_angles.NI5199_Hb.append(ratio_NI_Hb[kkangle])  
            ratio_angles.SIIab671631_Ha.append(ratio_SIIab_Ha[kkangle])
            ratio_angles.SIIa6716_SIIb6731.append(ratio_SIIa_SIIb[kkangle])
            ratio_angles.SIIab671631_SIIIab_631269.append(ratio_SIIab_SIIIab[kkangle]) 
            ratio_angles.OIbc630063_Ha.append(ratio_OIbc_Ha[kkangle])     
            ratio_angles.OIbc630063_OIa5577.append(ratio_OIbc_OIa[kkangle])   
            ratio_angles.OIbc630063_OIIIbc495907.append(ratio_OIbc_OIIIbc[kkangle])
            ratio_angles.OIbc630063_OIIcd732030.append(ratio_OIbc_OIIcd[kkangle])
            ratio_angles.OIbc630063_OIIab372729.append(ratio_OIbc_OIIab[kkangle])
            ratio_angles.OIb6300_Ha.append(ratio_OIb_Ha[kkangle])
            ratio_angles.OIIIc5007_Hb.append(ratio_OIIIc_Hb[kkangle])
            ratio_angles.OIIIbc495907_Hb.append(ratio_OIIIbc_Hb[kkangle])       
            ratio_angles.OIIIbc45907_OIIIa4363.append(ratio_OIIIbc_OIIIa[kkangle])   
            ratio_angles.OIIab372729_Hb.append(ratio_OIIab_Hb[kkangle])
            ratio_angles.OIIab372729_OIIIbc495907.append(ratio_OIIab_OIIIbc[kkangle]) 
            ratio_angles.OIIcd732030_OIIIbc495907.append(ratio_OIIcd_OIIIbc[kkangle])
            ratio_angles.ArIVa4712_ArIVb4740.append(ratio_ArIVa_ArIVb[kkangle])
            ratio_angles.ArIVab471240_Hb.append(ratio_ArIVab_Hb[kkangle])
            ratio_angles.NeIIIab386867_Hb.append(ratio_NeIIIab_Hb[kkangle])   
            ratio_angles.ClIIIa5517_ClIIIb5538.append(ratio_ClIIIa_ClIIIb[kkangle])
            ratio_angles.ClIIIab551738_Hb.append(ratio_ClIIIab_Hb[kkangle])
            ratio_angles.CI8727_Ha.append(ratio_CI_Ha[kkangle])
            ratio_angles.CII6461_Ha.append(ratio_CII_Ha[kkangle])
 
            ratio_angles_error.HeIa5876_Ha.append(ratio_HeIa_Ha_err[kkangle])
            ratio_angles_error.HeIb6678_Ha.append(ratio_HeIb_Ha_err[kkangle])
            ratio_angles_error.HeIIa4686_Hb.append(ratio_HeIIa_Hb_err[kkangle])
            ratio_angles_error.HeIIb5412_Hb.append(ratio_HeIIb_Hb_err[kkangle])
            ratio_angles_error.HeIa5876_HeIIa4686.append(ratio_HeIa_HeIIa_err[kkangle])
            ratio_angles_error.HeIa5876_HeIIb5412.append(ratio_HeIa_HeIIb_err[kkangle])
            ratio_angles_error.NIIc6583_Ha.append(ratio_NIIc_Ha_err[kkangle])
            ratio_angles_error.NIIbc654884_Ha.append(ratio_NIIbc_Ha_err[kkangle])
            ratio_angles_error.NIIbc654884_NIIa5755.append(ratio_NIIbc_NIIa_err[kkangle])
            ratio_angles_error.NIIbc654884_OIIIbc_495907.append(ratio_NIIbc_OIIIbc_err[kkangle])
            ratio_angles_error.NI5199_Hb.append(ratio_NI_Hb_err[kkangle])  
            ratio_angles_error.SIIab671631_Ha.append(ratio_SIIab_Ha_err[kkangle])
            ratio_angles_error.SIIa6716_SIIb6731.append(ratio_SIIa_SIIb_err[kkangle])
            ratio_angles_error.SIIab671631_SIIIab_631269.append(ratio_SIIab_SIIIab_err[kkangle]) 
            ratio_angles_error.OIbc630063_Ha.append(ratio_OIbc_Ha_err[kkangle])     
            ratio_angles_error.OIbc630063_OIa5577.append(ratio_OIbc_OIa_err[kkangle])   
            ratio_angles_error.OIbc630063_OIIIbc495907.append(ratio_OIbc_OIIIbc_err[kkangle])
            ratio_angles_error.OIbc630063_OIIcd732030.append(ratio_OIbc_OIIcd_err[kkangle])
            ratio_angles_error.OIbc630063_OIIab372729.append(ratio_OIbc_OIIab_err[kkangle])
            ratio_angles_error.OIb6300_Ha.append(ratio_OIb_Ha_err[kkangle])
            ratio_angles_error.OIIIc5007_Hb.append(ratio_OIIIc_Hb_err[kkangle])
            ratio_angles_error.OIIIbc495907_Hb.append(ratio_OIIIbc_Hb_err[kkangle])       
            ratio_angles_error.OIIIbc45907_OIIIa4363.append(ratio_OIIIbc_OIIIa_err[kkangle])   
            ratio_angles_error.OIIab372729_Hb.append(ratio_OIIab_Hb_err[kkangle])
            ratio_angles_error.OIIab372729_OIIIbc495907.append(ratio_OIIab_OIIIbc_err[kkangle]) 
            ratio_angles_error.OIIcd732030_OIIIbc495907.append(ratio_OIIcd_OIIIbc_err[kkangle])
            ratio_angles_error.ArIVa4712_ArIVb4740.append(ratio_ArIVa_ArIVb_err[kkangle])
            ratio_angles_error.ArIVab471240_Hb.append(ratio_ArIVab_Hb_err[kkangle])
            ratio_angles_error.NeIIIab386867_Hb.append(ratio_NeIIIab_Hb_err[kkangle])   
            ratio_angles_error.ClIIIa5517_ClIIIb5538.append(ratio_ClIIIa_ClIIIb_err[kkangle])
            ratio_angles_error.ClIIIab551738_Hb.append(ratio_ClIIIab_Hb_err[kkangle])
            ratio_angles_error.CI8727_Ha.append(ratio_CI_Ha_err[kkangle])
            ratio_angles_error.CII6461_Ha.append(ratio_CII_Ha_err[kkangle])
        
            print >>file1, "######################################"," slit PA  ", ang[kkangle] ," ##############################"
    
            
            print >>file1, '{:03.0f}'.format(ang[kkangle]),'{:04.3f}'.format(cHb[kkangle]),'{:07.3f}'.format(OIIa[kkangle]), '{:07.3f}'.format(OIIb[kkangle]),'{:09.3f}'.format(NeIIIa[kkangle]),'{:09.3f}'.format(NeIIIb[kkangle]),'{:06.2f}'.format(Hd[kkangle]),'{:06.2f}'.format(Hg[kkangle]),'{:08.3f}'.format(OIIIa[kkangle]),'{:08.3f}'.format(HeIIa[kkangle]),'{:08.3f}'.format(ArIVa[kkangle]),'{:08.3f}'.format(ArIVb[kkangle]),'{:06.2f}'.format(Hb[kkangle]),'{:08.3f}'.format(OIIIb[kkangle]),'{:08.3f}'.format(OIIIc[kkangle]),'{:06.3f}'.format(NI[kkangle]),'{:07.3f}'.format(HeIIb[kkangle]),'{:09.3f}'.format(ClIIIa[kkangle]),'{:09.3f}'.format(ClIIIb[kkangle]),'{:06.3f}'.format(OIa[kkangle]),'{:07.3f}'.format(NIIa[kkangle]),'{:07.3f}'.format(HeIa[kkangle]),'{:06.3f}'.format(OIb[kkangle]),'{:08.3f}'.format(SIIIa[kkangle]), '{:06.3f}'.format(OIc[kkangle]),'{:07.3f}'.format(CII[kkangle]),'{:07.3f}'.format(NIIb[kkangle]),'{:06.3f}'.format(Ha[kkangle]),'{:07.3f}'.format(NIIc[kkangle]),'{:07.3f}'.format(HeIb[kkangle]),'{:07.3f}'.format(SIIa[kkangle]),'{:07.3f}'.format(SIIb[kkangle]),'{:09.3f}'.format(ArIII[kkangle]),'{:07.3f}'.format(OIIc[kkangle]),'{:07.3f}'.format(OIId[kkangle]),'{:06.3f}'.format(CI[kkangle]),'{:08.3f}'.format(SIIIb[kkangle]),'{:05.4e}'.format(flux_angles.Hb_4861[kkangle]*10**(-20))
        
            print >>file1, '{:03.0f}'.format(ang[kkangle]),'{:04.3f}'.format(cHb_err[kkangle]),'{:07.3f}'.format(OIIa_err[kkangle]), '{:07.3f}'.format(OIIb_err[kkangle]),'{:09.3f}'.format(NeIIIa_err[kkangle]),'{:09.3f}'.format(NeIIIb_err[kkangle]),'{:06.2f}'.format(Hd_err[kkangle]),'{:06.2f}'.format(Hg_err[kkangle]),'{:08.3f}'.format(OIIIa_err[kkangle]),'{:08.3f}'.format(HeIIa_err[kkangle]),'{:08.3f}'.format(ArIVa_err[kkangle]),'{:08.3f}'.format(ArIVb_err[kkangle]),'{:06.2f}'.format(Hb_err[kkangle]),'{:08.3f}'.format(OIIIb_err[kkangle]),'{:08.3f}'.format(OIIIc_err[kkangle]),'{:06.3f}'.format(NI_err[kkangle]),'{:07.3f}'.format(HeIIb_err[kkangle]),'{:09.3f}'.format(ClIIIa_err[kkangle]),'{:09.3f}'.format(ClIIIb_err[kkangle]),'{:06.3f}'.format(OIa_err[kkangle]),'{:07.3f}'.format(NIIa_err[kkangle]),'{:07.3f}'.format(HeIa_err[kkangle]),'{:06.3f}'.format(OIb_err[kkangle]),'{:08.3f}'.format(SIIIa_err[kkangle]), '{:06.3f}'.format(OIc_err[kkangle]),'{:07.3f}'.format(CII_err[kkangle]),'{:07.3f}'.format(NIIb_err[kkangle]),'{:06.3f}'.format(Ha_err[kkangle]),'{:07.3f}'.format(NIIc_err[kkangle]),'{:07.3f}'.format(HeIb_err[kkangle]),'{:07.3f}'.format(SIIa_err[kkangle]),'{:07.3f}'.format(SIIb_err[kkangle]),'{:09.3f}'.format(ArIII_err[kkangle]),'{:07.3f}'.format(OIIc_err[kkangle]),'{:07.3f}'.format(OIId_err[kkangle]),'{:06.3f}'.format(CI_err[kkangle]),'{:08.3f}'.format(SIIIb_err[kkangle]),'{:05.4e}'.format(flux_angles_error.Hb_4861[kkangle]*10**(-20))
        
#fluxes in ergs/s/cm^2 are only the sum from the images.  
            print >>file2, "######################################"," slit PA  ", ang[kkangle] ," ##############################"

            print >>file2,'{:03.0f}'.format(ang[kkangle]),'{:04.3f}'.format(cHb[kkangle]),'{:04.3f}'.format(cHg[kkangle]),'{:04.3f}'.format(cHd[kkangle]),'{:05.3f}'.format(ratio_HeIa_Ha[kkangle]),'{:05.3f}'.format(ratio_HeIb_Ha[kkangle]),'{:05.3f}'.format(ratio_HeIIa_Hb[kkangle]),'{:05.3f}'.format(ratio_HeIIb_Hb[kkangle]),'{:05.3f}'.format(ratio_HeIa_HeIIa[kkangle]),'{:05.3f}'.format(ratio_HeIa_HeIIb[kkangle]),'{:05.3f}'.format(ratio_NIIc_Ha[kkangle]),'{:05.3f}'.format(ratio_NIIbc_Ha[kkangle]),'{:05.3f}'.format(ratio_NIIbc_NIIa[kkangle]),'{:05.3f}'.format(ratio_NIIbc_OIIIbc[kkangle]),'{:05.3f}'.format(ratio_NI_Hb[kkangle]),'{:05.3f}'.format(ratio_SIIab_Ha[kkangle]),'{:05.3f}'.format(ratio_SIIa_SIIb[kkangle]),'{:05.3f}'.format(ratio_SIIab_SIIIab[kkangle]),'{:05.3f}'.format(ratio_OIbc_Ha[kkangle]),'{:05.3f}'.format(ratio_OIbc_OIa[kkangle]),'{:05.3f}'.format(ratio_OIbc_OIIIbc[kkangle]),'{:05.3f}'.format(ratio_OIbc_OIIcd[kkangle]),'{:05.3f}'.format(ratio_OIbc_OIIab[kkangle]),'{:05.3f}'.format(ratio_OIb_Ha[kkangle]),'{:05.3f}'.format(ratio_OIIIc_Hb[kkangle]),'{:05.3f}'.format(ratio_OIIIbc_Hb[kkangle]),'{:05.3f}'.format(ratio_OIIIbc_OIIIa[kkangle]),'{:05.3f}'.format(ratio_OIIab_Hb[kkangle]),'{:05.3f}'.format(ratio_OIIab_OIIIbc[kkangle]),'{:05.3f}'.format(ratio_OIIcd_OIIIbc[kkangle]),'{:05.3f}'.format(ratio_ArIVa_ArIVb[kkangle]),'{:05.3f}'.format(ratio_ArIVab_Hb[kkangle]),'{:05.3f}'.format(ratio_NeIIIab_Hb[kkangle]),'{:05.3f}'.format(ratio_ClIIIa_ClIIIb[kkangle]),'{:05.3f}'.format(ratio_ClIIIab_Hb[kkangle]),'{:05.3f}'.format(ratio_CI_Ha[kkangle]),'{:05.3f}'.format(ratio_CII_Ha[kkangle])
        
            print >>file2,'{:03.0f}'.format(ang[kkangle]),'{:04.3f}'.format(cHb_err[kkangle]),'{:04.3f}'.format(cHg_err[kkangle]),'{:04.3f}'.format(cHd_err[kkangle]),'{:05.3f}'.format(ratio_HeIa_Ha_err[kkangle]),'{:05.3f}'.format(ratio_HeIb_Ha_err[kkangle]),'{:05.3f}'.format(ratio_HeIIa_Hb_err[kkangle]),'{:05.3f}'.format(ratio_HeIIb_Hb_err[kkangle]),'{:05.3f}'.format(ratio_HeIa_HeIIa_err[kkangle]),'{:05.3f}'.format(ratio_HeIa_HeIIb_err[kkangle]),'{:05.3f}'.format(ratio_NIIc_Ha_err[kkangle]),'{:05.3f}'.format(ratio_NIIbc_Ha_err[kkangle]),'{:05.3f}'.format(ratio_NIIbc_NIIa_err[kkangle]),'{:05.3f}'.format(ratio_NIIbc_OIIIbc_err[kkangle]),'{:05.3f}'.format(ratio_NI_Hb_err[kkangle]),'{:05.3f}'.format(ratio_SIIab_Ha_err[kkangle]),'{:05.3f}'.format(ratio_SIIa_SIIb_err[kkangle]),'{:05.3f}'.format(ratio_SIIab_SIIIab_err[kkangle]),'{:05.3f}'.format(ratio_OIbc_Ha_err[kkangle]),'{:05.3f}'.format(ratio_OIbc_OIa_err[kkangle]),'{:05.3f}'.format(ratio_OIbc_OIIIbc_err[kkangle]),'{:05.3f}'.format(ratio_OIbc_OIIcd_err[kkangle]),'{:05.3f}'.format(ratio_OIbc_OIIab_err[kkangle]),'{:05.3f}'.format(ratio_OIb_Ha_err[kkangle]),'{:05.3f}'.format(ratio_OIIIc_Hb_err[kkangle]),'{:05.3f}'.format(ratio_OIIIbc_Hb_err[kkangle]),'{:05.3f}'.format(ratio_OIIIbc_OIIIa_err[kkangle]),'{:05.3f}'.format(ratio_OIIab_Hb_err[kkangle]),'{:05.3f}'.format(ratio_OIIab_OIIIbc_err[kkangle]),'{:05.3f}'.format(ratio_OIIcd_OIIIbc_err[kkangle]),'{:05.3f}'.format(ratio_ArIVa_ArIVb_err[kkangle]),'{:05.3f}'.format(ratio_ArIVab_Hb_err[kkangle]),'{:05.3f}'.format(ratio_NeIIIab_Hb_err[kkangle]),'{:05.3f}'.format(ratio_ClIIIa_ClIIIb_err[kkangle]),'{:05.3f}'.format(ratio_ClIIIab_Hb_err[kkangle]),'{:05.3f}'.format(ratio_CI_Ha_err[kkangle]),'{:05.3f}'.format(ratio_CII_Ha_err[kkangle])
    

        kkangle=kkangle+1

    file1.close() 
    file2.close() 
    prs.print_output()   
    
    
    index_line1=param_estimated.index("c_Hb")
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("cHbeta")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    if param_requered[index_line1]=="yes" :
        fig,ax1=plt.subplots(figsize=(12,7))
        plt.errorbar(ang,cHb,yerr=cHb_err,label=r'c(H$\beta$)')
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel(r'c(H$\beta$)',fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_cHbeta_angles.png")
        fig.savefig("output_angles_plots/fig_cHbeta_angles.pdf")

    return flux_angles,flux_angles_error,flux_angles_norm,flux_angles_norm_error,ratio_angles,ratio_angles_error,ang
