# specificPA_line_fluxes_script.py:
# Calculates the line intensities (corrected for interstellar extiction and relative to Hb=100) 
# and the emission line ratios for the 10 specific pseudo-slits.
# The rourine returns the results in two ascii files (results output_linesintensities_specific_slit.txt 
# and output_lineratios_specific_slit.txt)
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
from satellite import specific_line_analysis_script as slas
from satellite import norm_flux_error_script as nfes



def specficPA_line_fluxes(flux2D,flux2D_error,line_names,line_ext_error,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values):
    
    
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

    index_pixel_scale=param_mod_name.index("pixel_scale")
    pixscale=float(param_model_values[index_pixel_scale])*0.01
    
    index_energy=param_mod_name.index("energy_parameter")
    energy=10**(param_model_values[index_energy])
    energy_conversion=float(energy*pixscale*pixscale)
    
    index_PA_slit=param_mod_name.index("PA_for_specific_slit_1")
    index_width_slit=param_mod_name.index("width_for_specific_slit_1")
    index_length_slit=param_mod_name.index("length_for_specific_slit_1")
    index_xcoor_slit=param_mod_name.index("x_coor_of_spec_slit_1")
    index_ycoor_slit=param_mod_name.index("y_coor_of_spec_slit_1")
    index_xcoor_CS=param_mod_name.index("x_coor_of_CS")
    index_ycoor_CS=param_mod_name.index("y_coor_of_CS")    

    dataHanew_err,dataHbnew_err,dataHgnew_err,dataHdnew_err,dataHeIanew_err,dataHeIbnew_err,dataHeIIanew_err,dataHeIIbnew_err,dataNIIanew_err,dataNIIbnew_err,dataNIIcnew_err,dataNInew_err,dataOIIIanew_err,dataOIIIbnew_err,dataOIIIcnew_err,dataOIIanew_err,dataOIIbnew_err,dataOIIcnew_err,dataOIIdnew_err,dataOIanew_err,dataOIbnew_err,dataOIcnew_err,dataSIIanew_err,dataSIIbnew_err,dataSIIIanew_err,dataSIIIbnew_err,dataClIIIanew_err,dataClIIIbnew_err,dataArIIInew_err,dataArIVanew_err,dataArIVbnew_err,dataCInew_err,dataCIInew_err,dataNeIIIanew_err,dataNeIIIbnew_err=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    
    dataHanew,dataHbnew,dataHgnew,dataHdnew,dataHeIanew,dataHeIbnew,dataHeIIanew,dataHeIIbnew,dataNIIanew,dataNIIbnew,dataNIIcnew,dataNInew,dataOIIIanew,dataOIIIbnew,dataOIIIcnew,dataOIIanew,dataOIIbnew,dataOIIcnew,dataOIIdnew,dataOIanew,dataOIbnew,dataOIcnew,dataSIIanew,dataSIIbnew,dataSIIIanew,dataSIIIbnew,dataClIIIanew,dataClIIIbnew,dataArIIInew,dataArIVanew,dataArIVbnew,dataCInew,dataCIInew,dataNeIIIanew,dataNeIIIbnew=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    ang=[]
    Ha,Hb,Hg,Hd,HeIa,HeIb,HeIIa,HeIIb,NIIa,NIIb,NIIc,NI,OIIIa,OIIIb,OIIIc,OIIa,OIIb,OIIc,OIId,OIa,OIb,OIc,SIIa,SIIb,SIIIa,SIIIb,ClIIIa,ClIIIb,ArIII,ArIVa,ArIVb,CI,CII,NeIIIa,NeIIIb=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    Ha_err,Hb_err,Hg_err,Hd_err,HeIa_err,HeIb_err,HeIIa_err,HeIIb_err,NIIa_err,NIIb_err,NIIc_err,NI_err,OIIIa_err,OIIIb_err,OIIIc_err,OIIa_err,OIIb_err,OIIc_err,OIId_err,OIa_err,OIb_err,OIc_err,SIIa_err,SIIb_err,SIIIa_err,SIIIb_err,ClIIIa_err,ClIIIb_err,ArIII_err,ArIVa_err,ArIVb_err,CI_err,CII_err,NeIIIa_err,NeIIIb_err=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
    cHb,cHg,cHd=[],[],[] 
    cHb_err,cHg_err,cHd_err=[],[],[]

    ratio_HeIa_Ha,ratio_HeIb_Ha,ratio_HeIIa_Hb,ratio_HeIIb_Hb,ratio_HeIa_HeIIa,ratio_HeIa_HeIIb,ratio_NIIc_Ha,ratio_NIIbc_Ha,ratio_NIIbc_NIIa,ratio_NIIbc_OIIIbc,ratio_NI_Hb,ratio_SIIab_Ha,ratio_SIIa_SIIb,ratio_SIIab_SIIIab,ratio_OIbc_Ha,ratio_OIbc_OIa,ratio_OIbc_OIIIbc,ratio_OIbc_OIIcd,ratio_OIbc_OIIab,ratio_OIb_Ha,ratio_OIIIc_Hb,ratio_OIIIbc_Hb,ratio_OIIIbc_OIIIa,ratio_OIIab_Hb,ratio_OIIab_OIIIbc,ratio_OIIcd_OIIIbc,ratio_ArIVa_ArIVb,ratio_ArIVab_Hb,ratio_NeIIIab_Hb,ratio_ClIIIa_ClIIIb,ratio_ClIIIab_Hb,ratio_CI_Ha,ratio_CII_Ha=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    ratio_HeIa_Ha_err,ratio_HeIb_Ha_err,ratio_HeIIa_Hb_err,ratio_HeIIb_Hb_err,ratio_HeIa_HeIIa_err,ratio_HeIa_HeIIb_err,ratio_NIIc_Ha_err,ratio_NIIbc_Ha_err,ratio_NIIbc_NIIa_err,ratio_NIIbc_OIIIbc_err,ratio_NI_Hb_err,ratio_SIIab_Ha_err,ratio_SIIa_SIIb_err,ratio_SIIab_SIIIab_err,ratio_OIbc_Ha_err,ratio_OIbc_OIa_err,ratio_OIbc_OIIIbc_err,ratio_OIbc_OIIcd_err,ratio_OIbc_OIIab_err,ratio_OIb_Ha_err,ratio_OIIIc_Hb_err,ratio_OIIIbc_Hb_err,ratio_OIIIbc_OIIIa_err,ratio_OIIab_Hb_err,ratio_OIIab_OIIIbc_err,ratio_OIIcd_OIIIbc_err,ratio_ArIVa_ArIVb_err,ratio_ArIVab_Hb_err,ratio_NeIIIab_Hb_err,ratio_ClIIIa_ClIIIb_err,ratio_ClIIIab_Hb_err,ratio_CI_Ha_err,ratio_CII_Ha_err=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

    file1 = open('output_linesintensities_specific_slit.txt','w')
    print >>file1,"ang ","cHb ","OII3727","OII3729","NeIII3868", "NeIII3967","Hd    ", "Hg    ","OIII4363","HeII4686","ArIV4712","ArIV4740","Hb    ","OIII4959","OIII5007","NI5200","HeI5412","ClIII5517","ClIII5538","OI5577","NII5755","HeI5876","OI6300","SIII6312","OI6363","CII6461","NII6548","Ha     ","NII6584","HeI6678","SII6716","SII6731","ArIII7136","OII7320","OII7330","CI8727","SIII9069",'Total f(Hb)'

    file2 = open('output_lineratios_specific_slit.txt','w')
    print >>file2,"ang","c_Hb ","c_Hg ","c_Hd ","log(HeI_5876s/HI_6563s)","log(HeI_6678s/HI_6563s)","log(HeII_5412s/HI_4861s)","log(HeII_4686s/HI_4861s)","log(HeI_5876s/HeII_4686s)","log(HeI_5876s/HeII_5412s)","log(N2_6583s/HI_6563s)","log(N2_6548s+N2_6583s)/HI_6563s)","log((N2_6548s+N2_6583s)/N2_5755s)","log((N2_6548s+N2_6583s)/(O3_4959s+O3_5007s))","log(N1_5199s/HI_4861s)","log((S2_6716s+S2_6731s)/HI_6563s)","log(S2_6716s/S2_6731s)","log((S2_6716s+S2_6731s)/(S3_6312s+S3_9069s))","log((O1_6300s+O1_6363s)/HI_6563s)","log((O1_6300s+O1_6363s)/O1_5577s)","log((O1_6300s+O1_6363s)/(O3_4959s+O3_5007s))","log((O1_6300s+O1_6363s)/(O2_7320s+O2_7330s))","log((O1_6300s+O1_6363s)/(O2_3727s+O2_3729s))","log(O1_6300s/HI_6563s)","log(O3_5007s/HI_4861s)","log((O3_4959s+O3_5007s)/HI_4861s)","log((O3_4959s+O3_5007s)/O3_4363s)","log((O2_3727s+O2_3729s)/HI_4861s)","log((O2_3727s+O2_3729s)/(O3_4959s+O3_5007s))","log((O2_7320s+O2_7330s)/(O3_4959s+O3_5007s))","log(Ar4_4712s/Ar4_4740s)","log((Ar4_4712s+Ar4_4740s)/HI_4861s)","log((Ne3_3868s+Ne3_3967s)/HI_4861s)","log(Cl3_5517s/Cl3_5538s)","log((Cl3_5517s+Cl3_5538s)/HI_4861s)","log(C1_8727s/HI_6563s)","log(C2_6461s/HI_6563s)"   


    xcoor_CS=param_model_values[index_xcoor_CS]
    ycoor_CS=param_model_values[index_ycoor_CS]
    # The fluxes of various emission lines requested by the used  
    # are calculated for specific pseudo-slits provided by the used
    # and specific slit characteristics. 
    # The flxues for all the slits are saved in arrays and files.

    slit_number=0
    for i in range(0,10,1): # for the 10 specific pseudo-slits defined by the user.
        PA_slit=param_model_values[index_PA_slit]
        width_slit=param_model_values[index_width_slit]
        length_slit=param_model_values[index_length_slit]
        xcoor_slit=param_model_values[index_xcoor_slit]
        ycoor_slit=param_model_values[index_ycoor_slit]
  
  
        dataHanew.append(0),dataHbnew.append(0),dataHgnew.append(0),dataHdnew.append(0),dataHeIanew.append(0),dataHeIbnew.append(0),dataHeIIanew.append(0),dataHeIIbnew.append(0),dataNIIanew.append(0),dataNIIbnew.append(0),dataNIIcnew.append(0),dataNInew.append(0),dataOIIIanew.append(0),dataOIIIbnew.append(0),dataOIIIcnew.append(0),dataOIIanew.append(0),dataOIIbnew.append(0),dataOIIcnew.append(0),dataOIIdnew.append(0),dataOIanew.append(0),dataOIbnew.append(0),dataOIcnew.append(0),dataSIIanew.append(0),dataSIIbnew.append(0),dataSIIIanew.append(0),dataSIIIbnew.append(0),dataClIIIanew.append(0),dataClIIIbnew.append(0),dataArIIInew.append(0),dataArIVanew.append(0),dataArIVbnew.append(0),dataCInew.append(0),dataCIInew.append(0),dataNeIIIanew.append(0),dataNeIIIbnew.append(0)
        dataHanew_err.append(0),dataHbnew_err.append(0),dataHgnew_err.append(0),dataHdnew_err.append(0),dataHeIanew_err.append(0),dataHeIbnew_err.append(0),dataHeIIanew_err.append(0),dataHeIIbnew_err.append(0),dataNIIanew_err.append(0),dataNIIbnew_err.append(0),dataNIIcnew_err.append(0),dataNInew_err.append(0),dataOIIIanew_err.append(0),dataOIIIbnew_err.append(0),dataOIIIcnew_err.append(0),dataOIIanew_err.append(0),dataOIIbnew_err.append(0),dataOIIcnew_err.append(0),dataOIIdnew_err.append(0),dataOIanew_err.append(0),dataOIbnew_err.append(0),dataOIcnew_err.append(0),dataSIIanew_err.append(0),dataSIIbnew_err.append(0),dataSIIIanew_err.append(0),dataSIIIbnew_err.append(0),dataClIIIanew_err.append(0),dataClIIIbnew_err.append(0),dataArIIInew_err.append(0),dataArIVanew_err.append(0),dataArIVbnew_err.append(0),dataCInew_err.append(0),dataCIInew_err.append(0),dataNeIIIanew_err.append(0),dataNeIIIbnew_err.append(0)
        
        ang.append(0)
    
        cHb.append(0),cHg.append(0),cHd.append(0)
        cHb_err.append(0),cHg_err.append(0),cHd_err.append(0)

        Ha.append(0),Hb.append(0),Hg.append(0),Hd.append(0),HeIa.append(0),HeIb.append(0),HeIIa.append(0),HeIIb.append(0),NIIa.append(0),NIIb.append(0),NIIc.append(0),NI.append(0),OIIIa.append(0),OIIIb.append(0),OIIIc.append(0),OIIa.append(0),OIIb.append(0),OIIc.append(0),OIId.append(0),OIa.append(0),OIb.append(0),OIc.append(0),SIIa.append(0),SIIb.append(0),SIIIa.append(0),SIIIb.append(0),ClIIIa.append(0),ClIIIb.append(0),ArIII.append(0),ArIVa.append(0),ArIVb.append(0),CI.append(0),CII.append(0),NeIIIa.append(0),NeIIIb.append(0)
    
        Ha_err.append(0),Hb_err.append(0),Hg_err.append(0),Hd_err.append(0),HeIa_err.append(0),HeIb_err.append(0),HeIIa_err.append(0),HeIIb_err.append(0),NIIa_err.append(0),NIIb_err.append(0),NIIc_err.append(0),NI_err.append(0),OIIIa_err.append(0),OIIIb_err.append(0),OIIIc_err.append(0),OIIa_err.append(0),OIIb_err.append(0),OIIc_err.append(0),OIId_err.append(0),OIa_err.append(0),OIb_err.append(0),OIc_err.append(0),SIIa_err.append(0),SIIb_err.append(0),SIIIa_err.append(0),SIIIb_err.append(0),ClIIIa_err.append(0),ClIIIb_err.append(0),ArIII_err.append(0),ArIVa_err.append(0),ArIVb_err.append(0),CI_err.append(0),CII_err.append(0),NeIIIa_err.append(0),NeIIIb_err.append(0)
   
        ratio_HeIa_Ha.append(-10e-10),ratio_HeIb_Ha.append(-10e-10),ratio_HeIIa_Hb.append(-10e-10),ratio_HeIIb_Hb.append(-10e-10),ratio_HeIa_HeIIa.append(-10e-10),ratio_HeIa_HeIIb.append(-10e-10),ratio_NIIc_Ha.append(-10e-10),ratio_NIIbc_Ha.append(-10e-10),ratio_NIIbc_NIIa.append(-10e-10),ratio_NIIbc_OIIIbc.append(-10e-10),ratio_NI_Hb.append(-10e-10),ratio_SIIab_Ha.append(-10e-10),ratio_SIIa_SIIb.append(-10e-10),ratio_SIIab_SIIIab.append(-10e-10),ratio_OIbc_Ha.append(-10e-10),ratio_OIbc_OIa.append(-10e-10),ratio_OIbc_OIIIbc.append(-10e-10),ratio_OIbc_OIIcd.append(-10e-10),ratio_OIbc_OIIab.append(-10e-10),ratio_OIb_Ha.append(-10e-10),ratio_OIIIc_Hb.append(-10e-10),ratio_OIIIbc_Hb.append(-10e-10),ratio_OIIIbc_OIIIa.append(-10e-10),ratio_OIIab_Hb.append(-10e-10),ratio_OIIab_OIIIbc.append(-10e-10),ratio_OIIcd_OIIIbc.append(-10e-10),ratio_ArIVa_ArIVb.append(-10e-10),ratio_ArIVab_Hb.append(-10e-10),ratio_NeIIIab_Hb.append(-10e-10),ratio_ClIIIa_ClIIIb.append(-10e-10),ratio_ClIIIab_Hb.append(-10e-10),ratio_CI_Ha.append(-10e-10),ratio_CII_Ha.append(-10e-10)

        ratio_HeIa_Ha_err.append(-10e-10),ratio_HeIb_Ha_err.append(-10e-10),ratio_HeIIa_Hb_err.append(-10e-10),ratio_HeIIb_Hb_err.append(-10e-10),ratio_HeIa_HeIIa_err.append(-10e-10),ratio_HeIa_HeIIb_err.append(-10e-10),ratio_NIIc_Ha_err.append(-10e-10),ratio_NIIbc_Ha_err.append(-10e-10),ratio_NIIbc_NIIa_err.append(-10e-10),ratio_NIIbc_OIIIbc_err.append(-10e-10),ratio_NI_Hb_err.append(-10e-10),ratio_SIIab_Ha_err.append(-10e-10),ratio_SIIa_SIIb_err.append(-10e-10),ratio_SIIab_SIIIab_err.append(-10e-10),ratio_OIbc_Ha_err.append(-10e-10),ratio_OIbc_OIa_err.append(-10e-10),ratio_OIbc_OIIIbc_err.append(-10e-10),ratio_OIbc_OIIcd_err.append(-10e-10),ratio_OIbc_OIIab_err.append(-10e-10),ratio_OIb_Ha_err.append(-10e-10),ratio_OIIIc_Hb_err.append(-10e-10),ratio_OIIIbc_Hb_err.append(-10e-10),ratio_OIIIbc_OIIIa_err.append(-10e-10),ratio_OIIab_Hb_err.append(-10e-10),ratio_OIIab_OIIIbc_err.append(-10e-10),ratio_OIIcd_OIIIbc_err.append(-10e-10),ratio_ArIVa_ArIVb_err.append(-10e-10),ratio_ArIVab_Hb_err.append(-10e-10),ratio_NeIIIab_Hb_err.append(-10e-10),ratio_ClIIIa_ClIIIb_err.append(-10e-10),ratio_ClIIIab_Hb_err.append(-10e-10),ratio_CI_Ha_err.append(-10e-10),ratio_CII_Ha_err.append(-10e-10)

#        dataHanew[slit_number],dataHanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.Ha_6563,flux2D_error.Ha_6563,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
#        dataHbnew[slit_number],dataHbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.Hb_4861,flux2D_error.Hb_4861,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)


        index_extra_error_line=line_names.index("HI_6563s")
        if lines_available[index_extra_error_line]=="yes":
            dataHanew[slit_number],dataHanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.Ha_6563,flux2D_error.Ha_6563,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
        else: 
            print("Halpha emission line is missing!. Program aborted!")
        
        index_extra_error_line=line_names.index("HI_4861s")
        if lines_available[index_extra_error_line]=="yes":
            dataHbnew[slit_number],dataHbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.Hb_4861,flux2D_error.Hb_4861,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)        
        else: 
            print("Hbeta emission line is missing!. Program aborted!")


        if dataHanew[slit_number]>0 and dataHbnew[slit_number] > 0:
 
            index_extra_error=line_names.index("HI_6563e")          
            if lines_available[index_extra_error]=="yes":
                if line_ext_error[index_extra_error]!=0:
                    errextraHa=float(line_ext_error[index_extra_error])
                    dataHanew_err[slit_number]=dataHanew_err[slit_number]+dataHanew[slit_number]*errextraHa
                else:
                    errextraHa=float(line_ext_error[index_extra_error])
                    dataHanew_err[slit_number]=dataHanew[slit_number]*errextraHa
 
            index_extra_error=line_names.index("HI_4861e")          
            if lines_available[index_extra_error]=="yes":
                if line_ext_error[index_extra_error]!=0:
                    errextraHb=float(line_ext_error[index_extra_error])
                    dataHbnew_err[slit_number]=dataHbnew_err[slit_number]+dataHbnew[slit_number]*errextraHb
                else:
                    errextraHb=float(line_ext_error[index_extra_error])
                    dataHbnew_err[slit_number]=dataHbnew[slit_number]*errextraHb
            
            
            index_line=line_names.index("HI_4340s")
            if lines_available[index_line]=="yes":
                dataHgnew[slit_number],dataHgnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.Hg_4340,flux2D_error.Hg_4340,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)    
                
                index_extra_error=line_names.index("HI_4340e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHg=float(line_ext_error[index_extra_error])
                        dataHgnew_err[slit_number]=dataHgnew_err[slit_number]+dataHgnew[slit_number]*errextraHg
                    else:
                        errextraHg=float(line_ext_error[index_extra_error])
                        dataHgnew_err[slit_number]=dataHgnew[slit_number]*errextraHg    
                
            index_line=line_names.index("HI_4101s")
            if lines_available[index_line]=="yes":
                dataHdnew[slit_number],dataHdnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.Hd_4101,flux2D_error.Hd_4101,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)    
                
                index_extra_error=line_names.index("HI_4101e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHd=float(line_ext_error[index_extra_error])
                        dataHdnew_err[slit_number]=dataHdnew_err[slit_number]+dataHdnew[slit_number]*errextraHd
                    else:
                        errextraHd=float(line_ext_error[index_extra_error])
                        dataHdnew_err[slit_number]=dataHdnew[slit_number]*errextraHd    
           
            index_line=line_names.index("HeI_5876s")
            if lines_available[index_line]=="yes":
                dataHeIanew[slit_number],dataHeIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.HeIa_5876,flux2D_error.HeIa_5876,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)                  
                
                index_extra_error=line_names.index("HeI_5876e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIa=float(line_ext_error[index_extra_error])
                        dataHeIanew_err[slit_number]=dataHeIanew_err[slit_number]+dataHeIanew[slit_number]*errextraHeIa
                    else:
                        errextraHeIa=float(line_ext_error[index_extra_error])
                        dataHeIanew_err[slit_number]=dataHeIanew[slit_number]*errextraHeIa    
            
            index_line=line_names.index("HeI_6678s")
            if lines_available[index_line]=="yes":
                dataHeIbnew[slit_number],dataHeIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.HeIb_6678,flux2D_error.HeIb_6678,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)                  
                
                index_extra_error=line_names.index("HeI_6678e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIb=float(line_ext_error[index_extra_error])
                        dataHeIbnew_err[slit_number]=dataHeIbnew_err[slit_number]+dataHeIbnew[slit_number]*errextraHeIb
                    else:
                        errextraHeIb=float(line_ext_error[index_extra_error])
                        dataHeIbnew_err[slit_number]=dataHeIbnew[slit_number]*errextraHeIb  
            
            index_line=line_names.index("HeII_4686s")
            if lines_available[index_line]=="yes":
                dataHeIIanew[slit_number],dataHeIIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.HeIIa_4686,flux2D_error.HeIIa_4686,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)                  
                
                index_extra_error=line_names.index("HeII_4686e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIIa=float(line_ext_error[index_extra_error])
                        dataHeIIanew_err[slit_number]=dataHeIIanew_err[slit_number]+dataHeIIanew[slit_number]*errextraHeIIa
                    else:
                        errextraHeIIa=float(line_ext_error[index_extra_error])
                        dataHeIIanew_err[slit_number]=dataHeIIanew[slit_number]*errextraHeIIa 
            
            index_line=line_names.index("HeII_5412s")
            if lines_available[index_line]=="yes":
                dataHeIIbnew[slit_number],dataHeIIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.HeIIb_5412,flux2D_error.HeIIb_5412,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)                  
                
                index_extra_error=line_names.index("HeII_5412e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIIb=float(line_ext_error[index_extra_error])
                        dataHeIIbnew_err[slit_number]=dataHeIIbnew_err[slit_number]+dataHeIIbnew[slit_number]*errextraHeIIb
                    else:
                        errextraHeIIb=float(line_ext_error[index_extra_error])
                        dataHeIIbnew_err[slit_number]=dataHeIIbnew[slit_number]*errextraHeIIb
           
            index_line=line_names.index("N2_5755s")
            if lines_available[index_line]=="yes":
                dataNIIanew[slit_number],dataNIIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.NIIa_5755,flux2D_error.NIIa_5755,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)      
                
                index_extra_error=line_names.index("N2_5755e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIa=float(line_ext_error[index_extra_error])
                        dataNIIanew_err[slit_number]=dataNIIanew_err[slit_number]+dataNIIanew[slit_number]*errextraNIIa
                    else:
                        errextraNIIa=float(line_ext_error[index_extra_error])
                        dataNIIanew_err[slit_number]=dataNIIanew[slit_number]*errextraNIIa
            
            index_line=line_names.index("N2_6548s")
            if lines_available[index_line]=="yes":
                dataNIIbnew[slit_number],dataNIIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.NIIb_6548,flux2D_error.NIIb_6548,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("N2_6548e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIb=float(line_ext_error[index_extra_error])
                        dataNIIbnew_err[slit_number]=dataNIIbnew_err[slit_number]+dataNIIbnew[slit_number]*errextraNIIb
                    else:
                        errextraNIIb=float(line_ext_error[index_extra_error])
                        dataNIIbnew_err[slit_number]=dataNIIbnew[slit_number]*errextraNIIb 

            index_line=line_names.index("N2_6583s")
            if lines_available[index_line]=="yes":
                dataNIIcnew[slit_number],dataNIIcnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.NIIc_6584,flux2D_error.NIIc_6584,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("N2_6583e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIc=float(line_ext_error[index_extra_error])
                        dataNIIcnew_err[slit_number]=dataNIIcnew_err[slit_number]+dataNIIcnew[slit_number]*errextraNIIc
                    else:
                        errextraNIIc=float(line_ext_error[index_extra_error])
                        dataNIIcnew_err[slit_number]=dataNIIcnew[slit_number]*errextraNIIc 
            
            index_line=line_names.index("N1_5199s")
            if lines_available[index_line]=="yes":
                dataNInew[slit_number],dataNInew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.NI_5199,flux2D_error.NI_5199,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("N1_5199e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNI=float(line_ext_error[index_extra_error])
                        dataNInew_err[slit_number]=dataNInew_err[slit_number]+dataNInew[slit_number]*errextraNI
                    else:
                        errextraNI=float(line_ext_error[index_extra_error])
                        dataNInew_err[slit_number]=dataNInew[slit_number]*errextraNI
           
            index_line=line_names.index("O3_4363s")
            if lines_available[index_line]=="yes":
                dataOIIIanew[slit_number],dataOIIIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIIIa_4363,flux2D_error.OIIIa_4363,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("O3_4363e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIa=float(line_ext_error[index_extra_error])
                        dataOIIIanew_err[slit_number]=dataOIIIanew_err[slit_number]+dataOIIIanew[slit_number]*errextraOIIIa
                    else:
                        errextraOIIIa=float(line_ext_error[index_extra_error])
                        dataOIIIanew_err[slit_number]=dataOIIIanew[slit_number]*errextraOIIIa
    
            index_line=line_names.index("O3_4959s")
            if lines_available[index_line]=="yes":
                dataOIIIbnew[slit_number],dataOIIIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIIIb_4959,flux2D_error.OIIIb_4959,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("O3_4959e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIb=float(line_ext_error[index_extra_error])
                        dataOIIIbnew_err[slit_number]=dataOIIIbnew_err[slit_number]+dataOIIIbnew[slit_number]*errextraOIIIb
                    else:
                        errextraOIIIb=float(line_ext_error[index_extra_error])
                        dataOIIIbnew_err[slit_number]=dataOIIIbnew[slit_number]*errextraOIIIb

            index_line=line_names.index("O3_5007s")
            if lines_available[index_line]=="yes":
                dataOIIIcnew[slit_number],dataOIIIcnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIIIc_5007,flux2D_error.OIIIc_5007,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 

                index_extra_error=line_names.index("O3_5007e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIc=float(line_ext_error[index_extra_error])
                        dataOIIIcnew_err[slit_number]=dataOIIIcnew_err[slit_number]+dataOIIIcnew[slit_number]*errextraOIIIc
                    else:
                        errextraOIIIc=float(line_ext_error[index_extra_error])
                        dataOIIIcnew_err[slit_number]=dataOIIIcnew[slit_number]*errextraOIIIc
            
            index_line=line_names.index("O2_3727s")
            if lines_available[index_line]=="yes":
                dataOIIanew[slit_number],dataOIIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIIa_3727,flux2D_error.OIIa_3727,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("O2_3727e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIa=float(line_ext_error[index_extra_error])
                        dataOIIanew_err[slit_number]=dataOIIanew_err[slit_number]+dataOIIanew[slit_number]*errextraOIIa
                    else:
                        errextraOIIa=float(line_ext_error[index_extra_error])
                        dataOIIanew_err[slit_number]=dataOIIanew[slit_number]*errextraOIIa
            
            index_line=line_names.index("O2_3729s")
            if lines_available[index_line]=="yes":
                dataOIIbnew[slit_number],dataOIIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIIb_3729,flux2D_error.OIIb_3729,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("O2_3729e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIb=float(line_ext_error[index_extra_error])
                        dataOIIbnew_err[slit_number]=dataOIIbnew_err[slit_number]+dataOIIbnew[slit_number]*errextraOIIb
                    else:
                        errextraOIIb=float(line_ext_error[index_extra_error])
                        dataOIIbnew_err[slit_number]=dataOIIbnew[slit_number]*errextraOIIb
           
            index_line=line_names.index("O2_7320s")
            if lines_available[index_line]=="yes":
                dataOIIcnew[slit_number],dataOIIcnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIIc_7320,flux2D_error.OIIc_7320,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
            
                index_extra_error=line_names.index("O2_7320e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIc=float(line_ext_error[index_extra_error])
                        dataOIIcnew_err[slit_number]=dataOIIcnew_err[slit_number]+dataOIIcnew[slit_number]*errextraOIIc
                    else:
                        errextraOIIc=float(line_ext_error[index_extra_error])
                        dataOIIcnew_err[slit_number]=dataOIIcnew[slit_number]*errextraOIIc
                        
            index_line=line_names.index("O2_7330s")
            if lines_available[index_line]=="yes":
                dataOIIdnew[slit_number],dataOIIdnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIId_7330,flux2D_error.OIId_7330,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("O2_7330e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIId=float(line_ext_error[index_extra_error])
                        dataOIIdnew_err[slit_number]=dataOIIdnew_err[slit_number]+dataOIIdnew[slit_number]*errextraOIId
                    else:
                        errextraOIId=float(line_ext_error[index_extra_error])
                        dataOIIdnew_err[slit_number]=dataOIIdnew[slit_number]*errextraOIId

            index_line=line_names.index("O1_5577s")
            if lines_available[index_line]=="yes":
                dataOIanew[slit_number],dataOIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIa_5577,flux2D_error.OIa_5577,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("O1_5577e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIa=float(line_ext_error[index_extra_error])
                        dataOIanew_err[slit_number]=dataOIanew_err[slit_number]+dataOIanew[slit_number]*errextraOIa
                    else:
                        errextraOIa=float(line_ext_error[index_extra_error])
                        dataOIanew_err[slit_number]=dataOIanew[slit_number]*errextraOIa
                
            index_line=line_names.index("O1_6300s")
            if lines_available[index_line]=="yes":
                dataOIbnew[slit_number],dataOIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIb_6300,flux2D_error.OIb_6300,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 

                index_extra_error=line_names.index("O1_6300e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIb=float(line_ext_error[index_extra_error])
                        dataOIbnew_err[slit_number]=dataOIbnew_err[slit_number]+dataOIbnew[slit_number]*errextraOIb
                    else:
                        errextraOIb=float(line_ext_error[index_extra_error])
                        dataOIbnew_err[slit_number]=dataOIbnew[slit_number]*errextraOIb
                    
            index_line=line_names.index("O1_6363s")
            if lines_available[index_line]=="yes":
                dataOIcnew[slit_number],dataOIcnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.OIc_6363,flux2D_error.OIc_6363,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("O1_6363e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIc=float(line_ext_error[index_extra_error])
                        dataOIcnew_err[slit_number]=dataOIcnew_err[slit_number]+dataOIcnew[slit_number]*errextraOIc
                    else:
                        errextraOIc=float(line_ext_error[index_extra_error])
                        dataOIcnew_err[slit_number]=dataOIcnew[slit_number]*errextraOIc
            
            index_line=line_names.index("S2_6716s")
            if lines_available[index_line]=="yes":
                dataSIIanew[slit_number],dataSIIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.SIIa_6716,flux2D_error.SIIa_6716,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i) 
                
                index_extra_error=line_names.index("S2_6716e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIa=float(line_ext_error[index_extra_error])
                        dataSIIanew_err[slit_number]=dataSIIanew_err[slit_number]+dataSIIanew[slit_number]*errextraSIIa
                    else:
                        errextraSIIa=float(line_ext_error[index_extra_error])
                        dataSIIanew_err[slit_number]=dataSIIanew[slit_number]*errextraSIIa
        
            index_line=line_names.index("S2_6731s")
            if lines_available[index_line]=="yes":
                dataSIIbnew[slit_number],dataSIIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.SIIb_6731,flux2D_error.SIIb_6731,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)

                index_extra_error=line_names.index("S2_6731e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIb=float(line_ext_error[index_extra_error])
                        dataSIIbnew_err[slit_number]=dataSIIbnew_err[slit_number]+dataSIIbnew[slit_number]*errextraSIIb
                    else:
                        errextraSIIb=float(line_ext_error[index_extra_error])
                        dataSIIbnew_err[slit_number]=dataSIIbnew[slit_number]*errextraSIIb
            
            index_line=line_names.index("S3_6312s")
            if lines_available[index_line]=="yes":
                dataSIIIanew[slit_number],dataSIIIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.SIIIa_6312,flux2D_error.SIIIa_6312,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
            
                index_extra_error=line_names.index("S3_6312e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIIa=float(line_ext_error[index_extra_error])
                        dataSIIIanew_err[slit_number]=dataSIIIanew_err[slit_number]+dataSIIIanew[slit_number]*errextraSIIIa
                    else:
                        errextraSIIIa=float(line_ext_error[index_extra_error])
                        dataSIIIanew_err[slit_number]=dataSIIIanew[slit_number]*errextraSIIIa
                        
            index_line=line_names.index("S3_9069s")
            if lines_available[index_line]=="yes":
                dataSIIIbnew[slit_number],dataSIIIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.SIIIb_9069,flux2D_error.SIIIb_9069,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
                
                index_extra_error=line_names.index("S3_9069e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIIb=float(line_ext_error[index_extra_error])
                        dataSIIIbnew_err[slit_number]=dataSIIIbnew_err[slit_number]+dataSIIIbnew[slit_number]*errextraSIIIb
                    else:
                        errextraSIIIb=float(line_ext_error[index_extra_error])
                        dataSIIIbnew_err[slit_number]=dataSIIIbnew[slit_number]*errextraSIIIb
                        
            index_line=line_names.index("Cl3_5517s")
            if lines_available[index_line]=="yes":
                dataClIIIanew[slit_number],dataClIIIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.ClIIIa_5517,flux2D_error.ClIIIa_5517,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
                
                index_extra_error=line_names.index("Cl3_5517e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraClIIIa=float(line_ext_error[index_extra_error])
                        dataClIIIanew_err[slit_number]=dataClIIIanew_err[slit_number]+dataClIIIanew[slit_number]*errextraClIIIa
                    else:
                        errextraClIIIa=float(line_ext_error[index_extra_error])
                        dataClIIIanew_err[slit_number]=dataClIIIanew[slit_number]*errextraClIIIa
                    
            index_line=line_names.index("Cl3_5538s")
            if lines_available[index_line]=="yes":
                dataClIIIbnew[slit_number],dataClIIIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.ClIIIb_5538,flux2D_error.ClIIIb_5538,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)

                index_extra_error=line_names.index("Cl3_5538e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraClIIIb=float(line_ext_error[index_extra_error])
                        dataClIIIbnew_err[slit_number]=dataClIIIbnew_err[slit_number]+dataClIIIbnew[slit_number]*errextraClIIIb
                    else:
                        errextraClIIIb=float(line_ext_error[index_extra_error])
                        dataClIIIbnew_err[slit_number]=dataClIIIbnew[slit_number]*errextraClIIIb

            index_line=line_names.index("Ar3_7136s")
            if lines_available[index_line]=="yes":
                dataArIIInew[slit_number],dataArIIInew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.ArIII_7136,flux2D_error.ArIII_7136,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
            
                index_extra_error=line_names.index("Ar3_7136e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIII=float(line_ext_error[index_extra_error])
                        dataArIIInew_err[slit_number]=dataArIIInew_err[slit_number]+dataArIIInew[slit_number]*errextraArIII
                    else:
                        errextraArIII=float(line_ext_error[index_extra_error])
                        dataArIIInew_err[slit_number]=dataArIIInew[slit_number]*errextraArIII
            
            index_line=line_names.index("Ar4_4712s")
            if lines_available[index_line]=="yes":
                dataArIVanew[slit_number],dataArIVanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.ArIVa_4712,flux2D_error.ArIVa_4712,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
            
                index_extra_error=line_names.index("Ar4_4712e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIVa=float(line_ext_error[index_extra_error])
                        dataArIVanew_err[slit_number]=dataArIVanew_err[slit_number]+dataArIVanew[slit_number]*errextraArIVa
                    else:
                        errextraArIVa=float(line_ext_error[index_extra_error])
                        dataArIVanew_err[slit_number]=dataArIVanew[slit_number]*errextraArIVa
                        
            index_line=line_names.index("Ar4_4740s")
            if lines_available[index_line]=="yes":
                dataArIVbnew[slit_number],dataArIVbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.ArIVb_4740,flux2D_error.ArIVb_4740,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
                
                index_extra_error=line_names.index("Ar4_4740e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIVb=float(line_ext_error[index_extra_error])
                        dataArIVbnew_err[slit_number]=dataArIVbnew_err[slit_number]+dataArIVbnew[slit_number]*errextraArIVb
                    else:
                        errextraArIVb=float(line_ext_error[index_extra_error])
                        dataArIVbnew_err[slit_number]=dataArIVbnew[slit_number]*errextraArIVb
                
            index_line=line_names.index("C1_8727s")
            if lines_available[index_line]=="yes":
                dataCInew[slit_number],dataCInew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.CI_8727,flux2D_error.CI_8727,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
                
                index_extra_error=line_names.index("C1_8727e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraCI=float(line_ext_error[index_extra_error])
                        dataCInew_err[slit_number]=dataCInew_err[slit_number]+dataCInew[slit_number]*errextraCI
                    else:
                        errextraCI=float(line_ext_error[index_extra_error])
                        dataCInew_err[slit_number]=dataCInew[slit_number]*errextraCI
           
            index_line=line_names.index("C2_6461s")
            if lines_available[index_line]=="yes":
                dataCIInew[slit_number],dataCIInew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.CII_6461,flux2D_error.CII_6461,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
                
                index_extra_error=line_names.index("C2_6461e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraCII=float(line_ext_error[index_extra_error])
                        dataCIInew_err[slit_number]=dataCIInew_err[slit_number]+dataCIInew[slit_number]*errextraCII
                    else:
                        errextraCII=float(line_ext_error[index_extra_error])
                        dataCIInew_err[slit_number]=dataCIInew[slit_number]*errextraCII
                
            index_line=line_names.index("Ne3_3868s")
            if lines_available[index_line]=="yes":
                dataNeIIIanew[slit_number],dataNeIIIanew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.NeIIIa_3868,flux2D_error.NeIIIa_3868,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
    
                index_extra_error=line_names.index("Ne3_3868e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNeIIIa=float(line_ext_error[index_extra_error])
                        dataNeIIIanew_err[slit_number]=dataNeIIIanew_err[slit_number]+dataNeIIIanew[slit_number]*errextraNeIIIa
                    else:
                        errextraNeIIIa=float(line_ext_error[index_extra_error])
                        dataNeIIIanew_err[slit_number]=dataNeIIIanew[slit_number]*errextraNeIIIa
    
            index_line=line_names.index("Ne3_3967s")
            if lines_available[index_line]=="yes":
                dataNeIIIbnew[slit_number],dataNeIIIbnew_err[slit_number],num_pix_Ha=slas.analysis_specific_line(flux2D.NeIIIb_3967,flux2D_error.NeIIIb_3967,PA_slit,flux2D.Ha_6563,flux2D.Hb_4861,width_slit,length_slit,xcoor_CS,ycoor_CS,xcoor_slit,ycoor_slit,param_estimated,param_requered,i)
            
                index_extra_error=line_names.index("Ne3_3967e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNeIIIb=float(line_ext_error[index_extra_error])
                        dataNeIIIbnew_err[slit_number]=dataNeIIIbnew_err[slit_number]+dataNeIIIbnew[slit_number]*errextraNeIIIb
                    else:
                        errextraNeIIIb=float(line_ext_error[index_extra_error])
                        dataNeIIIbnew_err[slit_number]=dataNeIIIbnew[slit_number]*errextraNeIIIb
 
 
 
            file3 = open('test.dat','w')
            print >>file3, "LINE test err"

            
            error=nfes.flux_Error(dataHanew[slit_number],dataHanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "H1r_6563A", '{:03.3f}'.format(dataHanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
            
            error=nfes.flux_Error(dataHbnew[slit_number],dataHbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "H1r_4861A", '{:03.3f}'.format(dataHbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
                
            error=nfes.flux_Error(dataHgnew[slit_number],dataHgnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "H1r_4341A", '{:03.3f}'.format(dataHgnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHdnew[slit_number],dataHdnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "H1r_4102A", '{:03.3f}'.format(dataHdnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHeIIanew[slit_number],dataHeIIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "He2r_4686A", '{:03.3f}'.format(dataHeIIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHeIIbnew[slit_number],dataHeIIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "He2r_5411A", '{:03.3f}'.format(dataHeIIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHeIanew[slit_number],dataHeIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "He1r_5876A", '{:03.3f}'.format(dataHeIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataHeIbnew[slit_number],dataHeIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "He1r_6678A", '{:03.3f}'.format(dataHeIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataSIIanew[slit_number],dataSIIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "S2_6716A", '{:03.3f}'.format(dataSIIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataSIIbnew[slit_number],dataSIIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "S2_6731A", '{:03.3f}'.format(dataSIIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataSIIIanew[slit_number],dataSIIIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "S3_6312A", '{:03.3f}'.format(dataSIIIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataSIIIbnew[slit_number],dataSIIIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "S3_9069A", '{:03.3f}'.format(dataSIIIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNInew[slit_number],dataNInew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "N1_5200A", '{:03.3f}'.format(dataNInew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNIIanew[slit_number],dataNIIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "N2_5755A", '{:03.3f}'.format(dataNIIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNIIbnew[slit_number],dataNIIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "N2_6548A", '{:03.3f}'.format(dataNIIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNIIcnew[slit_number],dataNIIcnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "N2_6584A", '{:03.3f}'.format(dataNIIcnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIanew[slit_number],dataOIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O1_5577A", '{:03.3f}'.format(dataOIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIbnew[slit_number],dataOIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O1_6300A", '{:03.3f}'.format(dataOIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIcnew[slit_number],dataOIcnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O1_6364A", '{:03.3f}'.format(dataOIcnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIanew[slit_number],dataOIIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O2_3726A", '{:03.3f}'.format(dataOIIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIbnew[slit_number],dataOIIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O2_3729A", '{:03.3f}'.format(dataOIIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIcnew[slit_number],dataOIIcnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O2_7319A+", '{:03.3f}'.format(dataOIIcnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIdnew[slit_number],dataOIIdnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O2_7330A+", '{:03.3f}'.format(dataOIIdnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIIanew[slit_number],dataOIIIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O3_4363A", '{:03.3f}'.format(dataOIIIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIIbnew[slit_number],dataOIIIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O3_4959A", '{:03.3f}'.format(dataOIIIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataOIIIcnew[slit_number],dataOIIIcnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "O3_5007A", '{:03.3f}'.format(dataOIIIcnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataClIIIanew[slit_number],dataClIIIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])        
            print >>file3, "Cl3_5518A", '{:03.3f}'.format(dataClIIIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
            
            error=nfes.flux_Error(dataClIIIbnew[slit_number],dataClIIIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "Cl3_5538A", '{:03.3f}'.format(dataClIIIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)        
        
            error=nfes.flux_Error(dataNeIIIanew[slit_number],dataNeIIIanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "Ne3_3869A", '{:03.3f}'.format(dataNeIIIanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataNeIIIbnew[slit_number],dataNeIIIbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "Ne3_3968A", '{:03.3f}'.format(dataNeIIIbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)                
        
            error=nfes.flux_Error(dataArIIInew[slit_number],dataArIIInew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "Ar3_7136A", '{:03.3f}'.format(dataArIIInew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataArIVanew[slit_number],dataArIVanew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "Ar4_4711A", '{:03.3f}'.format(dataArIVanew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataArIVbnew[slit_number],dataArIVbnew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "Ar4_4740A", '{:03.3f}'.format(dataArIVbnew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
        
            error=nfes.flux_Error(dataCInew[slit_number],dataCInew_err[slit_number],dataHbnew[slit_number],dataHbnew_err[slit_number])
            print >>file3, "C1_8728A", '{:03.3f}'.format(dataCInew[slit_number]*100/dataHbnew[slit_number]),'{:03.3f}'.format(error)
 
 
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
            
            cHb[slit_number]=float(obs.extinction.cHbeta[0])
            cHb_err[slit_number]=float(np.std(obs_fake.extinction.cHbeta))
            

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
 
        
        
            Ha[slit_number]=float(obs.getIntens()['H1r_6563A'])
            Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
            par1=float(obs.getError()['H1r_6563A'])
            par2=float(obs.getError()['H1r_4861A'])
            par3=float(corHa_er/corHa)
            Ha_err[slit_number]=Ha[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
            Hb_err[slit_number]=0.0
 
            
            index_extra_error=line_names.index("HI_4340e")          
            if lines_available[index_extra_error]=="yes":
                Hg[slit_number]=float(obs.getIntens()['H1r_4341A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['H1r_4341A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHg_er/corHg)
                Hg_err[slit_number]=Hg[slit_number]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("HI_4101e")          
            if lines_available[index_extra_error]=="yes":
                Hd[slit_number]=float(obs.getIntens()['H1r_4102A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['H1r_4102A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHd_er/corHd)
                Hd_err[slit_number]=Hd[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("HeII_4686e")          
            if lines_available[index_extra_error]=="yes":
                HeIIa[slit_number]=float(obs.getIntens()['He2r_4686A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['He2r_4686A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIIa_er/corHeIIa)
                HeIIa_err[slit_number]=HeIIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("HeII_5412e")          
            if lines_available[index_extra_error]=="yes":
                HeIIb[slit_number]=float(obs.getIntens()['He2r_5411A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['He2r_5411A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIIb_er/corHeIIb)
                HeIIb_err[slit_number]=HeIIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("HeI_5876e")          
            print index_extra_error,lines_available[index_extra_error]
            if lines_available[index_extra_error]=="yes":
                HeIa[slit_number]=float(obs.getIntens()['He1r_5876A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['He1r_5876A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIa_er/corHeIa)
                HeIa_err[slit_number]=HeIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
            
            index_extra_error=line_names.index("HeI_6678e")          
            if lines_available[index_extra_error]=="yes":
                HeIb[slit_number]=float(obs.getIntens()['He1r_6678A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['He1r_6678A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHeIb_er/corHeIb)
                HeIb_err[slit_number]=HeIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("N2_5755e")          
            if lines_available[index_extra_error]=="yes":
                NIIa[slit_number]=float(obs.getIntens()['N2_5755A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['N2_5755A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIa_er/corNIIa)
                NIIa_err[slit_number]=NIIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("N2_6548e")          
            if lines_available[index_extra_error]=="yes":
                NIIb[slit_number]=float(obs.getIntens()['N2_6548A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['N2_6548A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIb_er/corNIIb)
                NIIb_err[slit_number]=NIIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("N2_6583e")          
            if lines_available[index_extra_error]=="yes":
                NIIc[slit_number]=float(obs.getIntens()['N2_6584A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['N2_6584A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNIIc_er/corNIIc)
                NIIc_err[slit_number]=NIIc[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                
            index_extra_error=line_names.index("N1_5199e")          
            if lines_available[index_extra_error]=="yes":
                NI[slit_number]=float(obs.getIntens()['N1_5200A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['N1_5200A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNI_er/corNI)
                NI_err[slit_number]=NI[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("O3_4363e")          
            if lines_available[index_extra_error]=="yes":
                OIIIa[slit_number]=float(obs.getIntens()['O3_4363A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O3_4363A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIa_er/corOIIIa)
                OIIIa_err[slit_number]=OIIIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("O3_4959e")          
            if lines_available[index_extra_error]=="yes":
                OIIIb[slit_number]=float(obs.getIntens()['O3_4959A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O3_4959A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIb_er/corOIIIb)
                OIIIb_err[slit_number]=OIIIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("O3_5007e")          
            if lines_available[index_extra_error]=="yes":
                OIIIc[slit_number]=float(obs.getIntens()['O3_5007A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O3_5007A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIIc_er/corOIIIc)
                OIIIc_err[slit_number]=OIIIc[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
     
            index_extra_error=line_names.index("O2_3727e")          
            if lines_available[index_extra_error]=="yes":
                OIIa[slit_number]=float(obs.getIntens()['O2_3726A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O2_3726A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIa_er/corOIIa)
                OIIa_err[slit_number]=OIIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("O2_3729e")          
            if lines_available[index_extra_error]=="yes":
                OIIb[slit_number]=float(obs.getIntens()['O2_3729A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O2_3729A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIb_er/corOIIb)
                OIIb_err[slit_number]=OIIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
     
            index_extra_error=line_names.index("O2_7320e")          
            if lines_available[index_extra_error]=="yes":
                OIIc[slit_number]=float(obs.getIntens()['O2_7319A+'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O2_7319A+'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIIc_er/corOIIc)
                OIIc_err[slit_number]=OIIc[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
     
            index_extra_error=line_names.index("O2_7330e")          
            if lines_available[index_extra_error]=="yes":
                OIId[slit_number]=float(obs.getIntens()['O2_7330A+'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O2_7330A+'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIId_er/corOIId)
                OIId_err[slit_number]=OIId[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
     
            index_extra_error=line_names.index("O1_5577e")          
            if lines_available[index_extra_error]=="yes":
                OIa[slit_number]=float(obs.getIntens()['O1_5577A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O1_5577A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIa_er/corOIa)
                OIa_err[slit_number]=OIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("O1_6300e")          
            if lines_available[index_extra_error]=="yes":
                OIb[slit_number]=float(obs.getIntens()['O1_6300A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O1_6300A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIb_er/corOIb)
                OIb_err[slit_number]=OIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)   
                    
            index_extra_error=line_names.index("O1_6363e")          
            if lines_available[index_extra_error]=="yes":
                OIc[slit_number]=float(obs.getIntens()['O1_6364A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['O1_6364A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corOIc_er/corOIc)
                OIc_err[slit_number]=OIc[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                
            index_extra_error=line_names.index("S2_6716e")          
            if lines_available[index_extra_error]=="yes":
                SIIa[slit_number]=float(obs.getIntens()['S2_6716A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['S2_6716A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIa_er/corSIIa)
                SIIa_err[slit_number]=SIIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
            
            index_extra_error=line_names.index("S2_6731e")          
            if lines_available[index_extra_error]=="yes":
                SIIb[slit_number]=float(obs.getIntens()['S2_6731A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['S2_6731A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIb_er/corSIIb)
                SIIb_err[slit_number]=SIIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("S3_6312e")          
            if lines_available[index_extra_error]=="yes":
                SIIIa[slit_number]=float(obs.getIntens()['S3_6312A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['S3_6312A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIIa_er/corSIIIa)
                SIIIa_err[slit_number]=SIIIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("S3_9069e")          
            if lines_available[index_extra_error]=="yes":
                SIIIb[slit_number]=float(obs.getIntens()['S3_9069A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['S3_9069A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corSIIIb_er/corSIIIb)
                SIIIb_err[slit_number]=SIIIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
            
            index_extra_error=line_names.index("Cl3_5517e")          
            if lines_available[index_extra_error]=="yes":
                ClIIIa[slit_number]=float(obs.getIntens()['Cl3_5518A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Cl3_5518A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corClIIIa_er/corClIIIa)
                ClIIIa_err[slit_number]=ClIIIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("Cl3_5538e")          
            if lines_available[index_extra_error]=="yes":
                ClIIIb[slit_number]=float(obs.getIntens()['Cl3_5538A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Cl3_5538A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corClIIIb_er/corClIIIb)
                ClIIIb_err[slit_number]=ClIIIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("Ar3_7136e")          
            if lines_available[index_extra_error]=="yes":
                ArIII[slit_number]=float(obs.getIntens()['Ar3_7136A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ar3_7136A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIII_er/corArIII)
                ArIII_err[slit_number]=ArIII[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
                    
            index_extra_error=line_names.index("Ar4_4712e")          
            if lines_available[index_extra_error]=="yes":
                ArIVa[slit_number]=float(obs.getIntens()['Ar4_4711A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ar4_4711A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIVa_er/corArIVa)
                ArIVa_err[slit_number]=ArIVa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)
        
            index_extra_error=line_names.index("Ar4_4740e")          
            if lines_available[index_extra_error]=="yes":
                ArIVb[slit_number]=float(obs.getIntens()['Ar4_4740A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ar4_4740A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corArIVb_er/corArIVb)
                ArIVb_err[slit_number]=ArIVb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("C1_8727e")          
            if lines_available[index_extra_error]=="yes":
                CI[slit_number]=float(obs.getIntens()['C1_8728A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['C1_8728A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corCI_er/corCI)
                CI_err[slit_number]=CI[slit_number]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("C2_6461e")          
            if lines_available[index_extra_error]=="yes":
               CII[slit_number]=float(obs.getIntens()['H1r_6563A'])
               Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
               par1=float(obs.getError()['H1r_6563A'])
               par2=float(obs.getError()['H1r_4861A'])
               par3=float(corCII_er/corCII)
               CII_err[slit_number]=CII[slit_number]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("Ne3_3868e")          
            if lines_available[index_extra_error]=="yes":
                NeIIIa[slit_number]=float(obs.getIntens()['Ne3_3869A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ne3_3869A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNeIIIa_er/corNeIIIa)
                NeIIIa_err[slit_number]=NeIIIa[slit_number]*np.sqrt(par1**2+par2**2+par3**2)

            index_extra_error=line_names.index("Ne3_3967e")          
            if lines_available[index_extra_error]=="yes":
                NeIIIb[slit_number]=float(obs.getIntens()['Ne3_3968A'])
                Hb[slit_number]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['Ne3_3968A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corNeIIIb_er/corNeIIIb)
                NeIIIb_err[slit_number]=NeIIIb[slit_number]*np.sqrt(par1**2+par2**2+par3**2)


        
       

    


        flux_angles.Ha_6563.append(dataHanew[slit_number])
        flux_angles.Hb_4861.append(dataHbnew[slit_number])
        flux_angles.Hg_4340.append(dataHgnew[slit_number])
        flux_angles.Hd_4101.append(dataHdnew[slit_number])
        flux_angles.HeIa_5876.append(dataHeIanew[slit_number])
        flux_angles.HeIb_6678.append(dataHeIbnew[slit_number])
        flux_angles.HeIIa_4686.append(dataHeIIanew[slit_number])
        flux_angles.HeIIb_5412.append(dataHeIIbnew[slit_number])
        flux_angles.NIIa_5755.append(dataNIIanew[slit_number])
        flux_angles.NIIb_6548.append(dataNIIbnew[slit_number])
        flux_angles.NIIc_6584.append(dataNIIcnew[slit_number])
        flux_angles.NI_5199.append(dataNInew[slit_number])
        flux_angles.OIIIa_4363.append(dataOIIIanew[slit_number])
        flux_angles.OIIIb_4959.append(dataOIIIbnew[slit_number])
        flux_angles.OIIIc_5007.append(dataOIIIcnew[slit_number])
        flux_angles.OIIa_3727.append(dataOIIanew[slit_number])
        flux_angles.OIIb_3729.append(dataOIIbnew[slit_number])
        flux_angles.OIIc_7320.append(dataOIIcnew[slit_number])
        flux_angles.OIId_7330.append(dataOIIdnew[slit_number])
        flux_angles.OIa_5577.append(dataOIanew[slit_number])
        flux_angles.OIb_6300.append(dataOIbnew[slit_number])
        flux_angles.OIc_6363.append(dataOIcnew[slit_number])
        flux_angles.SIIa_6716.append(dataSIIanew[slit_number])
        flux_angles.SIIb_6731.append(dataSIIbnew[slit_number])
        flux_angles.SIIIa_6312.append(dataSIIIanew[slit_number])
        flux_angles.SIIIb_9069.append(dataSIIIbnew[slit_number])
        flux_angles.ClIIIa_5517.append(dataClIIIanew[slit_number])
        flux_angles.ClIIIb_5538.append(dataClIIIbnew[slit_number])
        flux_angles.ArIII_7136.append(dataArIIInew[slit_number])
        flux_angles.ArIVa_4712.append(dataArIVanew[slit_number])
        flux_angles.ArIVb_4740.append(dataArIVbnew[slit_number])
        flux_angles.CI_8727.append(dataCInew[slit_number])
        flux_angles.CII_6461.append(dataCIInew[slit_number])
        flux_angles.NeIIIa_3868.append(dataNeIIIanew[slit_number])
        flux_angles.NeIIIb_3967.append(dataNeIIIbnew[slit_number])

        flux_angles_error.Ha_6563.append(dataHanew_err[slit_number])
        flux_angles_error.Hb_4861.append(dataHbnew_err[slit_number])
        flux_angles_error.Hg_4340.append(dataHgnew_err[slit_number])
        flux_angles_error.Hd_4101.append(dataHdnew_err[slit_number])
        flux_angles_error.HeIa_5876.append(dataHeIanew_err[slit_number])
        flux_angles_error.HeIb_6678.append(dataHeIbnew_err[slit_number])
        flux_angles_error.HeIIa_4686.append(dataHeIIanew_err[slit_number])
        flux_angles_error.HeIIb_5412.append(dataHeIIbnew_err[slit_number])
        flux_angles_error.NIIa_5755.append(dataNIIanew_err[slit_number])
        flux_angles_error.NIIb_6548.append(dataNIIbnew_err[slit_number])
        flux_angles_error.NIIc_6584.append(dataNIIcnew_err[slit_number])
        flux_angles_error.NI_5199.append(dataNInew_err[slit_number])
        flux_angles_error.OIIIa_4363.append(dataOIIIanew_err[slit_number])
        flux_angles_error.OIIIb_4959.append(dataOIIIbnew_err[slit_number])
        flux_angles_error.OIIIc_5007.append(dataOIIIcnew_err[slit_number])
        flux_angles_error.OIIa_3727.append(dataOIIanew_err[slit_number])
        flux_angles_error.OIIb_3729.append(dataOIIbnew_err[slit_number])
        flux_angles_error.OIIc_7320.append(dataOIIcnew_err[slit_number])
        flux_angles_error.OIId_7330.append(dataOIIdnew_err[slit_number])
        flux_angles_error.OIa_5577.append(dataOIanew_err[slit_number])
        flux_angles_error.OIb_6300.append(dataOIbnew_err[slit_number])
        flux_angles_error.OIc_6363.append(dataOIcnew_err[slit_number])
        flux_angles_error.SIIa_6716.append(dataSIIanew_err[slit_number])
        flux_angles_error.SIIb_6731.append(dataSIIbnew_err[slit_number])
        flux_angles_error.SIIIa_6312.append(dataSIIIanew_err[slit_number])
        flux_angles_error.SIIIb_9069.append(dataSIIIbnew_err[slit_number])
        flux_angles_error.ClIIIa_5517.append(dataClIIIanew_err[slit_number])
        flux_angles_error.ClIIIb_5538.append(dataClIIIbnew_err[slit_number])
        flux_angles_error.ArIII_7136.append(dataArIIInew_err[slit_number])
        flux_angles_error.ArIVa_4712.append(dataArIVanew_err[slit_number])
        flux_angles_error.ArIVb_4740.append(dataArIVbnew_err[slit_number])
        flux_angles_error.CI_8727.append(dataCInew_err[slit_number])
        flux_angles_error.CII_6461.append(dataCIInew_err[slit_number])
        flux_angles_error.NeIIIa_3868.append(dataNeIIIanew_err[slit_number])
        flux_angles_error.NeIIIb_3967.append(dataNeIIIbnew_err[slit_number])

        flux_angles_norm.Ha_6563.append(Ha[slit_number])
        flux_angles_norm.Hb_4861.append(Hb[slit_number])
        flux_angles_norm.Hg_4340.append(Hg[slit_number])
        flux_angles_norm.Hd_4101.append(Hd[slit_number])
        flux_angles_norm.HeIa_5876.append(HeIa[slit_number])
        flux_angles_norm.HeIb_6678.append(HeIb[slit_number])
        flux_angles_norm.HeIIa_4686.append(HeIIa[slit_number])
        flux_angles_norm.HeIIb_5412.append(HeIIb[slit_number])
        flux_angles_norm.NIIa_5755.append(NIIa[slit_number])
        flux_angles_norm.NIIb_6548.append(NIIb[slit_number])
        flux_angles_norm.NIIc_6584.append(NIIc[slit_number])
        flux_angles_norm.NI_5199.append(NI[slit_number])
        flux_angles_norm.OIIIa_4363.append(OIIIa[slit_number])
        flux_angles_norm.OIIIb_4959.append(OIIIb[slit_number])
        flux_angles_norm.OIIIc_5007.append(OIIIc[slit_number])
        flux_angles_norm.OIIa_3727.append(OIIa[slit_number])
        flux_angles_norm.OIIb_3729.append(OIIb[slit_number])
        flux_angles_norm.OIIc_7320.append(OIIc[slit_number])
        flux_angles_norm.OIId_7330.append(OIId[slit_number])
        flux_angles_norm.OIa_5577.append(OIa[slit_number])
        flux_angles_norm.OIb_6300.append(OIb[slit_number])
        flux_angles_norm.OIc_6363.append(OIc[slit_number])
        flux_angles_norm.SIIa_6716.append(SIIa[slit_number])
        flux_angles_norm.SIIb_6731.append(SIIb[slit_number])
        flux_angles_norm.SIIIa_6312.append(SIIIa[slit_number])
        flux_angles_norm.SIIIb_9069.append(SIIIb[slit_number])
        flux_angles_norm.ClIIIa_5517.append(ClIIIa[slit_number])
        flux_angles_norm.ClIIIb_5538.append(ClIIIb[slit_number])
        flux_angles_norm.ArIII_7136.append(ArIII[slit_number])
        flux_angles_norm.ArIVa_4712.append(ArIVa[slit_number])
        flux_angles_norm.ArIVb_4740.append(ArIVb[slit_number])
        flux_angles_norm.CI_8727.append(CI[slit_number])
        flux_angles_norm.CII_6461.append(CII[slit_number])
        flux_angles_norm.NeIIIa_3868.append(NeIIIa[slit_number])
        flux_angles_norm.NeIIIb_3967.append(NeIIIb[slit_number])
        
        flux_angles_norm_error.Ha_6563.append(Ha_err[slit_number])
        flux_angles_norm_error.Hb_4861.append(Hb_err[slit_number])
        flux_angles_norm_error.Hg_4340.append(Hg_err[slit_number])
        flux_angles_norm_error.Hd_4101.append(Hd_err[slit_number])
        flux_angles_norm_error.HeIa_5876.append(HeIa_err[slit_number])
        flux_angles_norm_error.HeIb_6678.append(HeIb_err[slit_number])
        flux_angles_norm_error.HeIIa_4686.append(HeIIa_err[slit_number])
        flux_angles_norm_error.HeIIb_5412.append(HeIIb_err[slit_number])
        flux_angles_norm_error.NIIa_5755.append(NIIa_err[slit_number])
        flux_angles_norm_error.NIIb_6548.append(NIIb_err[slit_number])
        flux_angles_norm_error.NIIc_6584.append(NIIc_err[slit_number])
        flux_angles_norm_error.NI_5199.append(NI_err[slit_number])
        flux_angles_norm_error.OIIIa_4363.append(OIIIa_err[slit_number])
        flux_angles_norm_error.OIIIb_4959.append(OIIIb_err[slit_number])
        flux_angles_norm_error.OIIIc_5007.append(OIIIc_err[slit_number])
        flux_angles_norm_error.OIIa_3727.append(OIIa_err[slit_number])
        flux_angles_norm_error.OIIb_3729.append(OIIb_err[slit_number])
        flux_angles_norm_error.OIIc_7320.append(OIIc_err[slit_number])
        flux_angles_norm_error.OIId_7330.append(OIId_err[slit_number])
        flux_angles_norm_error.OIa_5577.append(OIa_err[slit_number])
        flux_angles_norm_error.OIb_6300.append(OIb_err[slit_number])
        flux_angles_norm_error.OIc_6363.append(OIc_err[slit_number])
        flux_angles_norm_error.SIIa_6716.append(SIIa_err[slit_number])
        flux_angles_norm_error.SIIb_6731.append(SIIb_err[slit_number])
        flux_angles_norm_error.SIIIa_6312.append(SIIIa_err[slit_number])
        flux_angles_norm_error.SIIIb_9069.append(SIIIb_err[slit_number])
        flux_angles_norm_error.ClIIIa_5517.append(ClIIIa_err[slit_number])
        flux_angles_norm_error.ClIIIb_5538.append(ClIIIb_err[slit_number])
        flux_angles_norm_error.ArIII_7136.append(ArIII_err[slit_number])
        flux_angles_norm_error.ArIVa_4712.append(ArIVa_err[slit_number])
        flux_angles_norm_error.ArIVb_4740.append(ArIVb_err[slit_number])
        flux_angles_norm_error.CI_8727.append(CI_err[slit_number])
        flux_angles_norm_error.CII_6461.append(CII_err[slit_number])
        flux_angles_norm_error.NeIIIa_3868.append(NeIIIa_err[slit_number])
        flux_angles_norm_error.NeIIIb_3967.append(NeIIIb_err[slit_number])
        

        index_par=param_estimated.index("log(HeI_5876s_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_HeIa_Ha[slit_number]=np.log10(HeIa[slit_number]/Ha[slit_number])
            par1= (HeIa_err[slit_number]/(HeIa[slit_number]*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_HeIa_Ha_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(HeI_6678s_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_HeIb_Ha[slit_number]=np.log10(HeIb[slit_number]/Ha[slit_number])
            par1= (HeIb_err[slit_number]/(HeIb[slit_number]*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_HeIb_Ha_err[slit_number]=np.sqrt(par1**2+par2**2)
 
        index_par=param_estimated.index("log(HeII_4686s_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_HeIIa_Hb[slit_number]=np.log10(HeIIa[slit_number]/Hb[slit_number])
            par1= (HeIIa_err[slit_number]/(HeIIa[slit_number]*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_HeIIa_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(HeII_5412s_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_HeIIb_Hb[slit_number]=np.log10(HeIIb[slit_number]/Hb[slit_number])    
            par1= (HeIIb_err[slit_number]/(HeIIb[slit_number]*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_HeIIb_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(HeI_5876s_HeII_4686s)")
        if param_requered[index_par]=="yes":
            ratio_HeIa_HeIIa[slit_number]=np.log10(HeIa[slit_number]/HeIIa[slit_number])    
            par1= (HeIa_err[slit_number]/(HeIa[slit_number]*np.log(10)))
            par2= (HeIIa_err[slit_number]/(HeIIa[slit_number]*np.log(10)))    
            ratio_HeIa_HeIIa_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(HeI_5876s_HeII_5412s)")
        if param_requered[index_par]=="yes":
            ratio_HeIa_HeIIb[slit_number]=np.log10(HeIa[slit_number]/HeIIb[slit_number])   
            par1= (HeIa_err[slit_number]/(HeIa[slit_number]*np.log(10)))
            par2= (HeIIb_err[slit_number]/(HeIIb[slit_number]*np.log(10)))    
            ratio_HeIa_HeIIb_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(N2_6583s_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_NIIc_Ha[slit_number]=np.log10(NIIc[slit_number]/Ha[slit_number])   
            par1= (NIIc_err[slit_number]/(NIIc[slit_number]*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_NIIc_Ha_err[slit_number]=np.sqrt(par1**2+par2**2) 
            
        index_par=param_estimated.index("log(N2_6548s+N2_6583s)_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_NIIbc_Ha[slit_number]=np.log10((NIIb[slit_number]+NIIc[slit_number])/Ha[slit_number])
            par1= ((NIIb_err[slit_number]+NIIc_err[slit_number])/((NIIb[slit_number]+NIIc[slit_number])*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_NIIbc_Ha_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((N2_6548s+N2_6583s)_N2_5755s)")
        if param_requered[index_par]=="yes":
            ratio_NIIbc_NIIa[slit_number]=np.log10((NIIb[slit_number]+NIIc[slit_number])/NIIa[slit_number])    
            par1= ((NIIb_err[slit_number]+NIIc_err[slit_number])/((NIIb[slit_number]+NIIc[slit_number])*np.log(10)))
            par2= (NIIa_err[slit_number]/(NIIa[slit_number]*np.log(10)))    
            ratio_NIIbc_NIIa_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((N2_6548s+N2_6583s)_(O3_4959s+O3_5007s))")
        if param_requered[index_par]=="yes":
            ratio_NIIbc_OIIIbc[slit_number]=np.log10((NIIb[slit_number]+NIIc[slit_number])/(OIIIb[slit_number]+OIIIc[slit_number]))     
            par1= ((NIIb_err[slit_number]+NIIc_err[slit_number])/((NIIb[slit_number]+NIIc[slit_number])*np.log(10)))
            par2= ((OIIIb_err[slit_number]+OIIIc_err[slit_number])/((OIIIb[slit_number]+OIIIc[slit_number])*np.log(10)))    
            ratio_NIIbc_OIIIbc_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(N1_5199s_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_NI_Hb[slit_number]=np.log10(NI[slit_number]/Hb[slit_number])  
            par1= (NI_err[slit_number]/(NI[slit_number]*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_NI_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)  
            
        index_par=param_estimated.index("log((S2_6716s+S2_6731s)_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_SIIab_Ha[slit_number]=np.log10((SIIa[slit_number]+SIIb[slit_number])/Ha[slit_number]) 
            par1= ((SIIa_err[slit_number]+SIIb_err[slit_number])/((SIIa[slit_number]+SIIb[slit_number])*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_SIIab_Ha_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(S2_6716s_S2_6731s)")
        if param_requered[index_par]=="yes":
            ratio_SIIa_SIIb[slit_number]=np.log10(SIIa[slit_number]/SIIb[slit_number]) 
            par1= (SIIa_err[slit_number]/(SIIa[slit_number]*np.log(10)))
            par2= (SIIb_err[slit_number]/(SIIb[slit_number]*np.log(10)))    
            ratio_SIIa_SIIb_err[slit_number]=np.sqrt(par1**2+par2**2)  
            
        index_par=param_estimated.index("log((S2_6716s+S2_6731s)_(S3_6312s+S3_9069s))")
        if param_requered[index_par]=="yes":
            ratio_SIIab_SIIIab[slit_number]=np.log10((SIIa[slit_number]+SIIb[slit_number])/(SIIIa[slit_number]+SIIIb[slit_number]))   
            par1= ((SIIa_err[slit_number]+SIIb_err[slit_number])/((SIIa[slit_number]+SIIb[slit_number])*np.log(10)))
            par2= ((SIIIa_err[slit_number]+SIIIb_err[slit_number])/((SIIIa[slit_number]+SIIIb[slit_number])*np.log(10)))    
            ratio_SIIab_SIIIab_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((O1_6300s+O1_6363s)_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_OIbc_Ha[slit_number]=np.log10((OIb[slit_number]+OIc[slit_number])/Ha[slit_number])     
            par1= ((OIb_err[slit_number]+OIc_err[slit_number])/((OIb[slit_number]+OIc[slit_number])*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_OIbc_Ha_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((O1_6300s+O1_6363s)_O1_5577s)")
        if param_requered[index_par]=="yes":
            ratio_OIbc_OIa[slit_number]=np.log10((OIb[slit_number]+OIc[slit_number])/OIa[slit_number])     
            par1= ((OIb_err[slit_number]+OIc_err[slit_number])/((OIb[slit_number]+OIc[slit_number])*np.log(10)))
            par2= (OIa_err[slit_number]/(OIa[slit_number]*np.log(10)))    
            ratio_OIbc_OIa_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((O1_6300s+O1_6363s)_(O3_4959s+O3_5007s))")
        if param_requered[index_par]=="yes":
            ratio_OIbc_OIIIbc[slit_number]=np.log10((OIb[slit_number]+OIc[slit_number])/(OIIIb[slit_number]+OIIIc[slit_number]))  
            par1= ((OIb_err[slit_number]+OIc_err[slit_number])/((OIb[slit_number]+OIc[slit_number])*np.log(10)))
            par2= ((OIIIb_err[slit_number]+OIIIc_err[slit_number])/((OIIIb[slit_number]+OIIIc[slit_number])*np.log(10)))    
            ratio_OIbc_OIIIbc_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((O1_6300s+O1_6363s)_(O2_7320s+O2_7330s))")
        if param_requered[index_par]=="yes":
            ratio_OIbc_OIIcd[slit_number]=np.log10((OIb[slit_number]+OIc[slit_number])/(OIIc[slit_number]+OIId[slit_number]))
            par1= ((OIb_err[slit_number]+OIc_err[slit_number])/((OIb[slit_number]+OIc[slit_number])*np.log(10)))
            par2= ((OIIc_err[slit_number]+OIId_err[slit_number])/((OIIc[slit_number]+OIId[slit_number])*np.log(10)))    
            ratio_OIbc_OIIcd_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((O1_6300s+O1_6363s)_(O2_3727s+O2_3729s))")
        if param_requered[index_par]=="yes":
            ratio_OIbc_OIIab[slit_number]=np.log10((OIb[slit_number]+OIc[slit_number])/(OIIa[slit_number]+OIIb[slit_number]))
            par1= ((OIb_err[slit_number]+OIc_err[slit_number])/((OIb[slit_number]+OIc[slit_number])*np.log(10)))
            par2= ((OIIa_err[slit_number]+OIIb_err[slit_number])/((OIIa[slit_number]+OIIb[slit_number])*np.log(10)))    
            ratio_OIbc_OIIab_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(O1_6300s_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_OIb_Ha[slit_number]=np.log10(OIb[slit_number]/Ha[slit_number])
            par1= (OIb_err[slit_number]/(OIb[slit_number]*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_OIb_Ha_err[slit_number]=np.sqrt(par1**2+par2**2)      
  
        index_par=param_estimated.index("log(O3_5007s_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_OIIIc_Hb[slit_number]=np.log10(OIIIc[slit_number]/Hb[slit_number])
            par1= (OIIIc_err[slit_number]/(OIIIc[slit_number]*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_OIIIc_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)      
  
        index_par=param_estimated.index("log((O3_4959s+O3_5007s)_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_OIIIbc_Hb[slit_number]=np.log10((OIIIb[slit_number]+OIIIc[slit_number])/Hb[slit_number])    
            par1= ((OIIIb_err[slit_number]+OIIIc_err[slit_number])/((OIIIb[slit_number]+OIIIc[slit_number])*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_OIIIbc_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((O3_4959s+O3_5007s)_O3_4363s)")
        if param_requered[index_par]=="yes":
            ratio_OIIIbc_OIIIa[slit_number]=np.log10((OIIIb[slit_number]+OIIIc[slit_number])/OIIIa[slit_number])     
            par1= ((OIIIb_err[slit_number]+OIIIc_err[slit_number])/((OIIIb[slit_number]+OIIIc[slit_number])*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_OIIIbc_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)
        
        index_par=param_estimated.index("log((O2_3727s+O2_3729s)_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_OIIab_Hb[slit_number]=np.log10((OIIa[slit_number]+OIIb[slit_number])/Hb[slit_number])       
            par1= ((OIIa_err[slit_number]+OIIb_err[slit_number])/((OIIa[slit_number]+OIIb[slit_number])*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_OIIab_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((O2_3727s+O2_3729s)_(O3_4959s+O3_5007s))")
        if param_requered[index_par]=="yes":
            ratio_OIIab_OIIIbc[slit_number]=np.log10((OIIa[slit_number]+OIIb[slit_number])/(OIIIb[slit_number]+OIIIc[slit_number]))  
            par1= ((OIIa_err[slit_number]+OIIb_err[slit_number])/((OIIa[slit_number]+OIIb[slit_number])*np.log(10)))
            par2= ((OIIIb_err[slit_number]+OIIIc_err[slit_number])/((OIIIb[slit_number]+OIIIc[slit_number])*np.log(10)))    
            ratio_OIIab_OIIIbc_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((O2_7320s+O2_7330s)_(O3_4959s+O3_5007s))")
        if param_requered[index_par]=="yes":
            ratio_OIIcd_OIIIbc[slit_number]=np.log10((OIIc[slit_number]+OIId[slit_number])/(OIIIb[slit_number]+OIIIc[slit_number]))  
            par1= ((OIIc_err[slit_number]+OIId_err[slit_number])/((OIIc[slit_number]+OIId[slit_number])*np.log(10)))
            par2= ((OIIIb_err[slit_number]+OIIIc_err[slit_number])/((OIIIb[slit_number]+OIIIc[slit_number])*np.log(10)))    
            ratio_OIIcd_OIIIbc_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log(Ar4_4712s_Ar4_4740s)")
        if param_requered[index_par]=="yes":
            ratio_ArIVa_ArIVb[slit_number]=np.log10(ArIVa[slit_number]/ArIVb[slit_number])  
            par1= (ArIVa_err[slit_number]/(ArIVa[slit_number]*np.log(10)))
            par2= (ArIVb_err[slit_number]/(ArIVb[slit_number]*np.log(10)))    
            ratio_ArIVa_ArIVb_err[slit_number]=np.sqrt(par1**2+par2**2)      
        
        index_par=param_estimated.index("log((Ar4_4712s+Ar4_4740s)_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_ArIVab_Hb[slit_number]=np.log10((ArIVa[slit_number]+ArIVb[slit_number])/Hb[slit_number])     
            par1= ((ArIVa_err[slit_number]+ArIVb_err[slit_number])/((ArIVa[slit_number]+ArIVb[slit_number])*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_ArIVab_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)
            
        index_par=param_estimated.index("log((Ne3_3868s+Ne3_3967s)_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_NeIIIab_Hb[slit_number]=np.log10((NeIIIa[slit_number]+NeIIIb[slit_number])/Hb[slit_number])     
            par1= ((NeIIIa_err[slit_number]+NeIIIb_err[slit_number])/((NeIIIa[slit_number]+NeIIIb[slit_number])*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_NeIIIab_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)
    
        index_par=param_estimated.index("log(Cl3_5517s_Cl3_5538s)")
        if param_requered[index_par]=="yes":
            ratio_ClIIIa_ClIIIb[slit_number]=np.log10(ClIIIa[slit_number]/ClIIIb[slit_number]) 
            par1= (ClIIIa_err[slit_number]/(ClIIIa[slit_number]*np.log(10)))
            par2= (ClIIIb_err[slit_number]/(ClIIIb[slit_number]*np.log(10)))    
            ratio_ClIIIa_ClIIIb_err[slit_number]=np.sqrt(par1**2+par2**2)      
  
        index_par=param_estimated.index("log((Cl3_5517s+Cl3_5538s)_HI_4861s)")
        if param_requered[index_par]=="yes":
            ratio_ClIIIab_Hb[slit_number]=np.log10((ClIIIa[slit_number]+ClIIIb[slit_number])/Hb[slit_number])  
            par1= ((ClIIIa_err[slit_number]+ClIIIb_err[slit_number])/((ClIIIa[slit_number]+ClIIIb[slit_number])*np.log(10)))
            par2= (Hb_err[slit_number]/(Hb[slit_number]*np.log(10)))    
            ratio_ClIIIab_Hb_err[slit_number]=np.sqrt(par1**2+par2**2)    
            
        index_par=param_estimated.index("log(C1_8727s_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_CI_Ha[slit_number]=np.log10(CI[slit_number]/Ha[slit_number]) 
            par1= (CI_err[slit_number]/(CI[slit_number]*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_CI_Ha_err[slit_number]=np.sqrt(par1**2+par2**2)      

        index_par=param_estimated.index("log(C2_6461s_HI_6563s)")
        if param_requered[index_par]=="yes":
            ratio_CII_Ha[slit_number]=np.log10(CII[slit_number]/Ha[slit_number]) 
            par1= (CII_err[slit_number]/(CII[slit_number]*np.log(10)))
            par2= (Ha_err[slit_number]/(Ha[slit_number]*np.log(10)))    
            ratio_CII_Ha_err[slit_number]=np.sqrt(par1**2+par2**2)
        
        ratio_angles.HeIa5876_Ha.append(ratio_HeIa_Ha[slit_number])
        ratio_angles.HeIb6678_Ha.append(ratio_HeIb_Ha[slit_number])
        ratio_angles.HeIIa4686_Hb.append(ratio_HeIIa_Hb[slit_number])
        ratio_angles.HeIIb5412_Hb.append(ratio_HeIIb_Hb[slit_number])
        ratio_angles.HeIa5876_HeIIa4686.append(ratio_HeIa_HeIIa[slit_number])
        ratio_angles.HeIa5876_HeIIb5412.append(ratio_HeIa_HeIIb[slit_number])
        ratio_angles.NIIc6583_Ha.append(ratio_NIIc_Ha[slit_number])
        ratio_angles.NIIbc654884_Ha.append(ratio_NIIbc_Ha[slit_number])
        ratio_angles.NIIbc654884_NIIa5755.append(ratio_NIIbc_NIIa[slit_number])
        ratio_angles.NIIbc654884_OIIIbc_495907.append(ratio_NIIbc_OIIIbc[slit_number])
        ratio_angles.NI5199_Hb.append(ratio_NI_Hb[slit_number])  
        ratio_angles.SIIab671631_Ha.append(ratio_SIIab_Ha[slit_number])
        ratio_angles.SIIa6716_SIIb6731.append(ratio_SIIa_SIIb[slit_number])
        ratio_angles.SIIab671631_SIIIab_631269.append(ratio_SIIab_SIIIab[slit_number]) 
        ratio_angles.OIbc630063_Ha.append(ratio_OIbc_Ha[slit_number])     
        ratio_angles.OIbc630063_OIa5577.append(ratio_OIbc_OIa[slit_number])   
        ratio_angles.OIbc630063_OIIIbc495907.append(ratio_OIbc_OIIIbc[slit_number])
        ratio_angles.OIbc630063_OIIcd732030.append(ratio_OIbc_OIIcd[slit_number])
        ratio_angles.OIbc630063_OIIab372729.append(ratio_OIbc_OIIab[slit_number])
        ratio_angles.OIb6300_Ha.append(ratio_OIb_Ha[slit_number])
        ratio_angles.OIIIc5007_Hb.append(ratio_OIIIc_Hb[slit_number])
        ratio_angles.OIIIbc495907_Hb.append(ratio_OIIIbc_Hb[slit_number])       
        ratio_angles.OIIIbc45907_OIIIa4363.append(ratio_OIIIbc_OIIIa[slit_number])   
        ratio_angles.OIIab372729_Hb.append(ratio_OIIab_Hb[slit_number])
        ratio_angles.OIIab372729_OIIIbc495907.append(ratio_OIIab_OIIIbc[slit_number]) 
        ratio_angles.OIIcd732030_OIIIbc495907.append(ratio_OIIcd_OIIIbc[slit_number])
        ratio_angles.ArIVa4712_ArIVb4740.append(ratio_ArIVa_ArIVb[slit_number])
        ratio_angles.ArIVab471240_Hb.append(ratio_ArIVab_Hb[slit_number])
        ratio_angles.NeIIIab386867_Hb.append(ratio_NeIIIab_Hb[slit_number])   
        ratio_angles.ClIIIa5517_ClIIIb5538.append(ratio_ClIIIa_ClIIIb[slit_number])
        ratio_angles.ClIIIab551738_Hb.append(ratio_ClIIIab_Hb[slit_number])
        ratio_angles.CI8727_Ha.append(ratio_CI_Ha[slit_number])
        ratio_angles.CII6461_Ha.append(ratio_CII_Ha[slit_number])
 
        ratio_angles_error.HeIa5876_Ha.append(ratio_HeIa_Ha_err[slit_number])
        ratio_angles_error.HeIb6678_Ha.append(ratio_HeIb_Ha_err[slit_number])
        ratio_angles_error.HeIIa4686_Hb.append(ratio_HeIIa_Hb_err[slit_number])
        ratio_angles_error.HeIIb5412_Hb.append(ratio_HeIIb_Hb_err[slit_number])
        ratio_angles_error.HeIa5876_HeIIa4686.append(ratio_HeIa_HeIIa_err[slit_number])
        ratio_angles_error.HeIa5876_HeIIb5412.append(ratio_HeIa_HeIIb_err[slit_number])
        ratio_angles_error.NIIc6583_Ha.append(ratio_NIIc_Ha_err[slit_number])
        ratio_angles_error.NIIbc654884_Ha.append(ratio_NIIbc_Ha_err[slit_number])
        ratio_angles_error.NIIbc654884_NIIa5755.append(ratio_NIIbc_NIIa_err[slit_number])
        ratio_angles_error.NIIbc654884_OIIIbc_495907.append(ratio_NIIbc_OIIIbc_err[slit_number])
        ratio_angles_error.NI5199_Hb.append(ratio_NI_Hb_err[slit_number])  
        ratio_angles_error.SIIab671631_Ha.append(ratio_SIIab_Ha_err[slit_number])
        ratio_angles_error.SIIa6716_SIIb6731.append(ratio_SIIa_SIIb_err[slit_number])
        ratio_angles_error.SIIab671631_SIIIab_631269.append(ratio_SIIab_SIIIab_err[slit_number]) 
        ratio_angles_error.OIbc630063_Ha.append(ratio_OIbc_Ha_err[slit_number])     
        ratio_angles_error.OIbc630063_OIa5577.append(ratio_OIbc_OIa_err[slit_number])   
        ratio_angles_error.OIbc630063_OIIIbc495907.append(ratio_OIbc_OIIIbc_err[slit_number])
        ratio_angles_error.OIbc630063_OIIcd732030.append(ratio_OIbc_OIIcd_err[slit_number])
        ratio_angles_error.OIbc630063_OIIab372729.append(ratio_OIbc_OIIab_err[slit_number])
        ratio_angles_error.OIb6300_Ha.append(ratio_OIb_Ha_err[slit_number])
        ratio_angles_error.OIIIc5007_Hb.append(ratio_OIIIc_Hb_err[slit_number])
        ratio_angles_error.OIIIbc495907_Hb.append(ratio_OIIIbc_Hb_err[slit_number])       
        ratio_angles_error.OIIIbc45907_OIIIa4363.append(ratio_OIIIbc_OIIIa_err[slit_number])   
        ratio_angles_error.OIIab372729_Hb.append(ratio_OIIab_Hb_err[slit_number])
        ratio_angles_error.OIIab372729_OIIIbc495907.append(ratio_OIIab_OIIIbc_err[slit_number]) 
        ratio_angles_error.OIIcd732030_OIIIbc495907.append(ratio_OIIcd_OIIIbc_err[slit_number])
        ratio_angles_error.ArIVa4712_ArIVb4740.append(ratio_ArIVa_ArIVb_err[slit_number])
        ratio_angles_error.ArIVab471240_Hb.append(ratio_ArIVab_Hb_err[slit_number])
        ratio_angles_error.NeIIIab386867_Hb.append(ratio_NeIIIab_Hb_err[slit_number])   
        ratio_angles_error.ClIIIa5517_ClIIIb5538.append(ratio_ClIIIa_ClIIIb_err[slit_number])
        ratio_angles_error.ClIIIab551738_Hb.append(ratio_ClIIIab_Hb_err[slit_number])
        ratio_angles_error.CI8727_Ha.append(ratio_CI_Ha_err[slit_number])
        ratio_angles_error.CII6461_Ha.append(ratio_CII_Ha_err[slit_number])
        
        print >>file1, "######################################"," specific slit number: ", slit_number+1 ," ##############################"

    
        print >>file1, '{:03.0f}'.format(ang[slit_number]),'{:04.3f}'.format(cHb[slit_number]),'{:07.3f}'.format(OIIa[slit_number]), '{:07.3f}'.format(OIIb[slit_number]),'{:09.3f}'.format(NeIIIa[slit_number]),'{:09.3f}'.format(NeIIIb[slit_number]),'{:06.2f}'.format(Hd[slit_number]),'{:06.2f}'.format(Hg[slit_number]),'{:08.3f}'.format(OIIIa[slit_number]),'{:08.3f}'.format(HeIIa[slit_number]),'{:08.3f}'.format(ArIVa[slit_number]),'{:08.3f}'.format(ArIVb[slit_number]),'{:06.2f}'.format(Hb[slit_number]),'{:08.3f}'.format(OIIIb[slit_number]),'{:08.3f}'.format(OIIIc[slit_number]),'{:06.3f}'.format(NI[slit_number]),'{:07.3f}'.format(HeIIb[slit_number]),'{:09.3f}'.format(ClIIIa[slit_number]),'{:09.3f}'.format(ClIIIb[slit_number]),'{:06.3f}'.format(OIa[slit_number]),'{:07.3f}'.format(NIIa[slit_number]),'{:07.3f}'.format(HeIa[slit_number]),'{:06.3f}'.format(OIb[slit_number]),'{:08.3f}'.format(SIIIa[slit_number]), '{:06.3f}'.format(OIc[slit_number]),'{:07.3f}'.format(CII[slit_number]),'{:07.3f}'.format(NIIb[slit_number]),'{:06.3f}'.format(Ha[slit_number]),'{:07.3f}'.format(NIIc[slit_number]),'{:07.3f}'.format(HeIb[slit_number]),'{:07.3f}'.format(SIIa[slit_number]),'{:07.3f}'.format(SIIb[slit_number]),'{:09.3f}'.format(ArIII[slit_number]),'{:07.3f}'.format(OIIc[slit_number]),'{:07.3f}'.format(OIId[slit_number]),'{:06.3f}'.format(CI[slit_number]),'{:08.3f}'.format(SIIIb[slit_number]),'{:05.4e}'.format(flux_angles.Hb_4861[i]*energy_conversion)

        print >>file1, '{:03.0f}'.format(ang[slit_number]),'{:04.3f}'.format(cHb_err[slit_number]),'{:07.3f}'.format(OIIa_err[slit_number]), '{:07.3f}'.format(OIIb_err[slit_number]),'{:09.3f}'.format(NeIIIa_err[slit_number]),'{:09.3f}'.format(NeIIIb_err[slit_number]),'{:06.2f}'.format(Hd_err[slit_number]),'{:06.2f}'.format(Hg_err[slit_number]),'{:08.3f}'.format(OIIIa_err[slit_number]),'{:08.3f}'.format(HeIIa_err[slit_number]),'{:08.3f}'.format(ArIVa_err[slit_number]),'{:08.3f}'.format(ArIVb_err[slit_number]),'{:06.2f}'.format(Hb_err[slit_number]),'{:08.3f}'.format(OIIIb_err[slit_number]),'{:08.3f}'.format(OIIIc_err[slit_number]),'{:06.3f}'.format(NI_err[slit_number]),'{:07.3f}'.format(HeIIb_err[slit_number]),'{:09.3f}'.format(ClIIIa_err[slit_number]),'{:09.3f}'.format(ClIIIb_err[slit_number]),'{:06.3f}'.format(OIa_err[slit_number]),'{:07.3f}'.format(NIIa_err[slit_number]),'{:07.3f}'.format(HeIa_err[slit_number]),'{:06.3f}'.format(OIb_err[slit_number]),'{:08.3f}'.format(SIIIa_err[slit_number]), '{:06.3f}'.format(OIc_err[slit_number]),'{:07.3f}'.format(CII_err[slit_number]),'{:07.3f}'.format(NIIb_err[slit_number]),'{:06.3f}'.format(Ha_err[slit_number]),'{:07.3f}'.format(NIIc_err[slit_number]),'{:07.3f}'.format(HeIb_err[slit_number]),'{:07.3f}'.format(SIIa_err[slit_number]),'{:07.3f}'.format(SIIb_err[slit_number]),'{:09.3f}'.format(ArIII_err[slit_number]),'{:07.3f}'.format(OIIc_err[slit_number]),'{:07.3f}'.format(OIId_err[slit_number]),'{:06.3f}'.format(CI_err[slit_number]),'{:08.3f}'.format(SIIIb_err[slit_number]),'{:05.4e}'.format(flux_angles_error.Hb_4861[i]*energy_conversion)
        
#fluxes in ergs/s/cm^2 are the sum from the images.  
        print >>file2, "######################################"," specific slit number: ", slit_number+1 ," ##############################"

        print >>file2,'{:03.0f}'.format(ang[slit_number]),'{:04.3f}'.format(cHb[slit_number]),'{:04.3f}'.format(cHg[slit_number]),'{:04.3f}'.format(cHd[slit_number]),'{:05.3f}'.format(ratio_HeIa_Ha[slit_number]),'{:05.3f}'.format(ratio_HeIb_Ha[slit_number]),'{:05.3f}'.format(ratio_HeIIa_Hb[slit_number]),'{:05.3f}'.format(ratio_HeIIb_Hb[slit_number]),'{:05.3f}'.format(ratio_HeIa_HeIIa[slit_number]),'{:05.3f}'.format(ratio_HeIa_HeIIb[slit_number]),'{:05.3f}'.format(ratio_NIIc_Ha[slit_number]),'{:05.3f}'.format(ratio_NIIbc_Ha[slit_number]),'{:05.3f}'.format(ratio_NIIbc_NIIa[slit_number]),'{:05.3f}'.format(ratio_NIIbc_OIIIbc[slit_number]),'{:05.3f}'.format(ratio_NI_Hb[slit_number]),'{:05.3f}'.format(ratio_SIIab_Ha[slit_number]),'{:05.3f}'.format(ratio_SIIa_SIIb[slit_number]),'{:05.3f}'.format(ratio_SIIab_SIIIab[slit_number]),'{:05.3f}'.format(ratio_OIbc_Ha[slit_number]),'{:05.3f}'.format(ratio_OIbc_OIa[slit_number]),'{:05.3f}'.format(ratio_OIbc_OIIIbc[slit_number]),'{:05.3f}'.format(ratio_OIbc_OIIcd[slit_number]),'{:05.3f}'.format(ratio_OIbc_OIIab[slit_number]),'{:05.3f}'.format(ratio_OIb_Ha[slit_number]),'{:05.3f}'.format(ratio_OIIIc_Hb[slit_number]),'{:05.3f}'.format(ratio_OIIIbc_Hb[slit_number]),'{:05.3f}'.format(ratio_OIIIbc_OIIIa[slit_number]),'{:05.3f}'.format(ratio_OIIab_Hb[slit_number]),'{:05.3f}'.format(ratio_OIIab_OIIIbc[slit_number]),'{:05.3f}'.format(ratio_OIIcd_OIIIbc[slit_number]),'{:05.3f}'.format(ratio_ArIVa_ArIVb[slit_number]),'{:05.3f}'.format(ratio_ArIVab_Hb[slit_number]),'{:05.3f}'.format(ratio_NeIIIab_Hb[slit_number]),'{:05.3f}'.format(ratio_ClIIIa_ClIIIb[slit_number]),'{:05.3f}'.format(ratio_ClIIIab_Hb[slit_number]),'{:05.3f}'.format(ratio_CI_Ha[slit_number]),'{:05.3f}'.format(ratio_CII_Ha[slit_number])
    
        print >>file2,'{:03.0f}'.format(ang[slit_number]),'{:04.3f}'.format(cHb_err[slit_number]),'{:04.3f}'.format(cHg_err[slit_number]),'{:04.3f}'.format(cHd_err[slit_number]),'{:05.3f}'.format(ratio_HeIa_Ha_err[slit_number]),'{:05.3f}'.format(ratio_HeIb_Ha_err[slit_number]),'{:05.3f}'.format(ratio_HeIIa_Hb_err[slit_number]),'{:05.3f}'.format(ratio_HeIIb_Hb_err[slit_number]),'{:05.3f}'.format(ratio_HeIa_HeIIa_err[slit_number]),'{:05.3f}'.format(ratio_HeIa_HeIIb_err[slit_number]),'{:05.3f}'.format(ratio_NIIc_Ha_err[slit_number]),'{:05.3f}'.format(ratio_NIIbc_Ha_err[slit_number]),'{:05.3f}'.format(ratio_NIIbc_NIIa_err[slit_number]),'{:05.3f}'.format(ratio_NIIbc_OIIIbc_err[slit_number]),'{:05.3f}'.format(ratio_NI_Hb_err[slit_number]),'{:05.3f}'.format(ratio_SIIab_Ha_err[slit_number]),'{:05.3f}'.format(ratio_SIIa_SIIb_err[slit_number]),'{:05.3f}'.format(ratio_SIIab_SIIIab_err[slit_number]),'{:05.3f}'.format(ratio_OIbc_Ha_err[slit_number]),'{:05.3f}'.format(ratio_OIbc_OIa_err[slit_number]),'{:05.3f}'.format(ratio_OIbc_OIIIbc_err[slit_number]),'{:05.3f}'.format(ratio_OIbc_OIIcd_err[slit_number]),'{:05.3f}'.format(ratio_OIbc_OIIab_err[slit_number]),'{:05.3f}'.format(ratio_OIb_Ha_err[slit_number]),'{:05.3f}'.format(ratio_OIIIc_Hb_err[slit_number]),'{:05.3f}'.format(ratio_OIIIbc_Hb_err[slit_number]),'{:05.3f}'.format(ratio_OIIIbc_OIIIa_err[slit_number]),'{:05.3f}'.format(ratio_OIIab_Hb_err[slit_number]),'{:05.3f}'.format(ratio_OIIab_OIIIbc_err[slit_number]),'{:05.3f}'.format(ratio_OIIcd_OIIIbc_err[slit_number]),'{:05.3f}'.format(ratio_ArIVa_ArIVb_err[slit_number]),'{:05.3f}'.format(ratio_ArIVab_Hb_err[slit_number]),'{:05.3f}'.format(ratio_NeIIIab_Hb_err[slit_number]),'{:05.3f}'.format(ratio_ClIIIa_ClIIIb_err[slit_number]),'{:05.3f}'.format(ratio_ClIIIab_Hb_err[slit_number]),'{:05.3f}'.format(ratio_CI_Ha_err[slit_number]),'{:05.3f}'.format(ratio_CII_Ha_err[slit_number])
        
    
        slit_number=slit_number+1
        index_PA_slit=index_PA_slit+5
        index_width_slit=index_width_slit+5
        index_length_slit=index_length_slit+5
        index_xcoor_slit=index_xcoor_slit+5
        index_ycoor_slit= index_ycoor_slit+5


    file1.close() 
    file2.close() 
    prs.print_output()   

    return flux_angles,flux_angles_error,flux_angles_norm,flux_angles_norm_error,ratio_angles,ratio_angles_error
