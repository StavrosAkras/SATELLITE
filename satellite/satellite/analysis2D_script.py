# analysis2D_script.py:
# Construct all the emission line diagnostic diagrams
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
import seaborn as sns 
import pyneb as pn

from satellite import diagnotic_diagrams_script as dds
from satellite import generate_2D_lineratio_maps_script as g2Dlrs
from satellite import histogram_script as hs
from satellite import calculations_excluding_outliers_script as ceos 
from satellite import norm_flux_error_script as nfes


        
def analysis2D(flux2D,flux2D_error,flux_angles_norm,ang,line_names,line_ext_error,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,DD_name,DD_avail,DDxmin,DDxmax,DDymin,DDymax,hdr,flux_spec_slit_norm,flux_spec_slit_norm_error):
    
    index_size=param_mod_name.index("total_num_pixels_horiz")
    sizex=param_model_values[index_size]
    sizey=param_model_values[index_size]            
    
    
    class flux_norm_2D:
        Ha_6563=np.zeros((sizex, sizey))
        Hb_4861=np.zeros((sizex, sizey))
        Hg_4340=np.zeros((sizex, sizey))
        Hd_4101=np.zeros((sizex, sizey))
        HeIa_5876=np.zeros((sizex, sizey))
        HeIb_6678=np.zeros((sizex, sizey))
        HeIIa_4686=np.zeros((sizex, sizey))
        HeIIb_5412=np.zeros((sizex, sizey))
        NIIa_5755=np.zeros((sizex, sizey))
        NIIb_6548=np.zeros((sizex, sizey))
        NIIc_6584=np.zeros((sizex, sizey))
        NI_5199=np.zeros((sizex, sizey))
        OIIIa_4363=np.zeros((sizex, sizey))
        OIIIb_4959=np.zeros((sizex, sizey))
        OIIIc_5007=np.zeros((sizex, sizey))
        OIIa_3727=np.zeros((sizex, sizey))
        OIIb_3729=np.zeros((sizex, sizey))
        OIIc_7320=np.zeros((sizex, sizey))
        OIId_7330=np.zeros((sizex, sizey))
        OIa_5577=np.zeros((sizex, sizey))
        OIb_6300=np.zeros((sizex, sizey))
        OIc_6363=np.zeros((sizex, sizey))
        SIIa_6716=np.zeros((sizex, sizey))
        SIIb_6731=np.zeros((sizex, sizey))
        SIIIa_6312=np.zeros((sizex, sizey))
        SIIIb_9069=np.zeros((sizex, sizey))
        ClIIIa_5517=np.zeros((sizex, sizey))
        ClIIIb_5538=np.zeros((sizex, sizey))
        ArIII_7136=np.zeros((sizex, sizey))
        ArIVa_4712=np.zeros((sizex, sizey))
        ArIVb_4740=np.zeros((sizex, sizey))
        CI_8727=np.zeros((sizex, sizey))
        CII_6461=np.zeros((sizex, sizey))
        NeIIIa_3868=np.zeros((sizex, sizey))
        NeIIIb_3967=np.zeros((sizex, sizey))

    class flux_norm_2D_error:
        Ha_6563=np.zeros((sizex, sizey))
        Hb_4861=np.zeros((sizex, sizey))
        Hg_4340=np.zeros((sizex, sizey))
        Hd_4101=np.zeros((sizex, sizey))
        HeIa_5876=np.zeros((sizex, sizey))
        HeIb_6678=np.zeros((sizex, sizey))
        HeIIa_4686=np.zeros((sizex, sizey))
        HeIIb_5412=np.zeros((sizex, sizey))
        NIIa_5755=np.zeros((sizex, sizey))
        NIIb_6548=np.zeros((sizex, sizey))
        NIIc_6584=np.zeros((sizex, sizey))
        NI_5199=np.zeros((sizex, sizey))
        OIIIa_4363=np.zeros((sizex, sizey))
        OIIIb_4959=np.zeros((sizex, sizey))
        OIIIc_5007=np.zeros((sizex, sizey))
        OIIa_3727=np.zeros((sizex, sizey))
        OIIb_3729=np.zeros((sizex, sizey))
        OIIc_7320=np.zeros((sizex, sizey))
        OIId_7330=np.zeros((sizex, sizey))
        OIa_5577=np.zeros((sizex, sizey))
        OIb_6300=np.zeros((sizex, sizey))
        OIc_6363=np.zeros((sizex, sizey))
        SIIa_6716=np.zeros((sizex, sizey))
        SIIb_6731=np.zeros((sizex, sizey))
        SIIIa_6312=np.zeros((sizex, sizey))
        SIIIb_9069=np.zeros((sizex, sizey))
        ClIIIa_5517=np.zeros((sizex, sizey))
        ClIIIb_5538=np.zeros((sizex, sizey))
        ArIII_7136=np.zeros((sizex, sizey))
        ArIVa_4712=np.zeros((sizex, sizey))
        ArIVb_4740=np.zeros((sizex, sizey))
        CI_8727=np.zeros((sizex, sizey))
        CII_6461=np.zeros((sizex, sizey))
        NeIIIa_3868=np.zeros((sizex, sizey))
        NeIIIb_3967=np.zeros((sizex, sizey))

        
        
    pn.log_.open_file('my_log.txt')
    
    Chbeta_pyneb=np.zeros((sizex, sizey))
    Chbeta_error_pyneb=np.zeros((sizex, sizey))
#########################################################################
### compute the line intensities (corrected for intestellar extinction)
### normalized to Hbeta=100 for all the pixels in the 2D maps
#########################################################################
    
    for j in range(0,sizey):
        for i in range(0,sizex):
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.Ha_6563[i,j]>2.85*flux2D.Hb_4861[i,j]):
                print("spaxles:",i,j)
                file3 = open('test.dat','w')
                
                print >>file3, "LINE test err"
                    
                index_extra_error=line_names.index("HI_6563e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHa=float(line_ext_error[index_extra_error])
                        flux2D_error.Ha_6563[i,j]=flux2D_error.Ha_6563[i,j]+flux2D.Ha_6563[i,j]*errextraHa
                    else:
                        errextraHa=float(line_ext_error[index_extra_error])
                        flux2D_error.Ha_6563[i,j]=flux2D.Ha_6563[i,j]*errextraHa

                index_extra_error=line_names.index("HI_4861e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHb=float(line_ext_error[index_extra_error])
                        flux2D_error.Hb_4861[i,j]=flux2D_error.Hb_4861[i,j]+flux2D.Hb_4861[i,j]*errextraHb
                    else:
                        errextraHb=float(line_ext_error[index_extra_error])
                        flux2D_error.Hb_4861[i,j]=flux2D.Hb_4861[i,j]*errextraHb
                error=nfes.flux_Error(flux2D.Ha_6563[i,j],flux2D_error.Ha_6563[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "H1r_6563A", '{:03.3f}'.format(flux2D.Ha_6563[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
                error=nfes.flux_Error(flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "H1r_4861A", '{:03.3f}'.format(flux2D.Hb_4861[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
                

                index_extra_error=line_names.index("HI_4340e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHg=float(line_ext_error[index_extra_error])
                        flux2D_error.Hg_4340[i,j]=flux2D_error.Hg_4340[i,j]+flux2D.Hg_4340[i,j]*errextraHg
                    else:
                        errextraHg=float(line_ext_error[index_extra_error])
                        flux2D_error.Hg_4340[i,j]=flux2D.Hg_4340[i,j]*errextraHg  

                error=nfes.flux_Error(flux2D.Hg_4340[i,j],flux2D_error.Hg_4340[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "H1r_4341A", '{:03.3f}'.format(flux2D.Hg_4340[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "H1r_4341Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("HI_4101e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHd=float(line_ext_error[index_extra_error])
                        flux2D_error.Hd_4101[i,j]=flux2D_error.Hd_4101[i,j]+flux2D.Hd_4101[i,j]*errextraHd
                    else:
                        errextraHd=float(line_ext_error[index_extra_error])
                        flux2D_error.Hd_4101[i,j]=flux2D.Hd_4101[i,j]*errextraHd
                        
                error=nfes.flux_Error(flux2D.Hd_4101[i,j],flux2D_error.Hd_4101[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "H1r_4102A", '{:03.3f}'.format(flux2D.Hd_4101[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "H1r_4102Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("HeII_4686e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.HeIIa_4686[i,j]=flux2D_error.HeIIa_4686[i,j]+flux2D.HeIIa_4686[i,j]*errextraHeIIa
                    else:
                        errextraHeIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.HeIIa_4686[i,j]=flux2D.HeIIa_4686[i,j]*errextraHeIIa 
                error=nfes.flux_Error(flux2D.HeIIa_4686[i,j],flux2D_error.HeIIa_4686[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "He2r_4686A", '{:03.3f}'.format(flux2D.HeIIa_4686[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "He2r_4686Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("HeII_5412e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.HeIIb_5412[i,j]=flux2D_error.HeIIb_5412[i,j]+flux2D.HeIIb_5412[i,j]*errextraHeIIb
                    else:
                        errextraHeIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.HeIIb_5412[i,j]=flux2D.HeIIb_5412[i,j]*errextraHeIIb
                        
                error=nfes.flux_Error(flux2D.HeIIb_5412[i,j],flux2D_error.HeIIb_5412[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "He2r_5411A", '{:03.3f}'.format(flux2D.HeIIb_5412[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "He2r_5411Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("HeI_5876e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIa=float(line_ext_error[index_extra_error])
                        flux2D_error.HeIa_5876[i,j]=flux2D_error.HeIa_5876[i,j]+flux2D.HeIa_5876[i,j]*errextraHeIa
                    else:
                        errextraHeIa=float(line_ext_error[index_extra_error])
                        flux2D_error.HeIa_5876[i,j]=flux2D.HeIa_5876[i,j]*errextraHeIa
                        
                error=nfes.flux_Error(flux2D.HeIa_5876[i,j],flux2D_error.HeIa_5876[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "He1r_5876A", '{:03.3f}'.format(flux2D.HeIa_5876[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "He1r_5876Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("HeI_6678e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraHeIb=float(line_ext_error[index_extra_error])
                        flux2D_error.HeIb_6678[i,j]=flux2D_error.HeIb_6678[i,j]+flux2D.HeIb_6678[i,j]*errextraHeIb
                    else:
                        errextraHeIb=float(line_ext_error[index_extra_error])
                        flux2D_error.HeIb_6678[i,j]=flux2D.HeIb_6678[i,j]*errextraHeIb 
                        
                error=nfes.flux_Error(flux2D.HeIb_6678[i,j],flux2D_error.HeIb_6678[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "He1r_6678A", '{:03.3f}'.format(flux2D.HeIb_6678[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "He1r_6678Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("S2_6716e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.SIIa_6716[i,j]=flux2D_error.SIIa_6716[i,j]+flux2D.SIIa_6716[i,j]*errextraSIIa
                    else:
                        errextraSIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.SIIa_6716[i,j]=flux2D.SIIa_6716[i,j]*errextraSIIa
                error=nfes.flux_Error(flux2D.SIIa_6716[i,j],flux2D_error.SIIa_6716[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "S2_6716A", '{:03.3f}'.format(flux2D.SIIa_6716[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "S2_6716Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("S2_6731e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.SIIb_6731[i,j]=flux2D_error.SIIb_6731[i,j]+flux2D.SIIb_6731[i,j]*errextraSIIb
                    else:
                        errextraSIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.SIIb_6731[i,j]=flux2D.SIIb_6731[i,j]*errextraSIIb
                error=nfes.flux_Error(flux2D.SIIb_6731[i,j],flux2D_error.SIIb_6731[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "S2_6731A", '{:03.3f}'.format(flux2D.SIIb_6731[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "S2_6731Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("S3_6312e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.SIIIa_6312[i,j]=flux2D_error.SIIIa_6312[i,j]+flux2D.SIIIa_6312[i,j]*errextraSIIIa
                    else:
                        errextraSIIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.SIIIa_6312[i,j]=flux2D.SIIIa_6312[i,j]*errextraSIIIa
                error=nfes.flux_Error(flux2D.SIIIa_6312[i,j],flux2D_error.SIIIa_6312[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "S3_6312A", '{:03.3f}'.format(flux2D.SIIIa_6312[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "S3_6312Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("S3_9069e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraSIIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.SIIIb_9069[i,j]=flux2D_error.SIIIb_9069[i,j]+flux2D.SIIIb_9069[i,j]*errextraSIIIb
                    else:
                        errextraSIIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.SIIIb_9069[i,j]=flux2D.SIIIb_9069[i,j]*errextraSIIIb
                error=nfes.flux_Error(flux2D.SIIIb_9069[i,j],flux2D_error.SIIIb_9069[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "S3_9069A", '{:03.3f}'.format(flux2D.SIIIb_9069[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "S3_9069Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("N1_5199e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNI=float(line_ext_error[index_extra_error])
                        flux2D_error.NI_5199[i,j]=flux2D_error.NI_5199[i,j]+flux2D.NI_5199[i,j]*errextraNI
                    else:
                        errextraNI=float(line_ext_error[index_extra_error])
                        flux2D_error.NI_5199[i,j]=flux2D.NI_5199[i,j]*errextraNI
                error=nfes.flux_Error(flux2D.NI_5199[i,j],flux2D_error.NI_5199[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "N1_5200A", '{:03.3f}'.format(flux2D.NI_5199[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "N1_5200Ae", '{:03.3f}'.format(error)
        
                index_extra_error=line_names.index("N2_5755e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.NIIa_5755[i,j]=flux2D_error.NIIa_5755[i,j]+flux2D.NIIa_5755[i,j]*errextraNIIa
                    else:
                        errextraNIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.NIIa_5755[i,j]=flux2D.NIIa_5755[i,j]*errextraNIIa
                error=nfes.flux_Error(flux2D.NIIa_5755[i,j],flux2D_error.NIIa_5755[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "N2_5755A", '{:03.3f}'.format(flux2D.NIIa_5755[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
        
                index_extra_error=line_names.index("N2_6548e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.NIIb_6548[i,j]=flux2D_error.NIIb_6548[i,j]+flux2D.NIIb_6548[i,j]*errextraNIIb
                    else:
                        errextraNIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.NIIb_6548[i,j]=flux2D.NIIb_6548[i,j]*errextraNIIb 
                error=nfes.flux_Error(flux2D.NIIb_6548[i,j],flux2D_error.NIIb_6548[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "N2_6548A", '{:03.3f}'.format(flux2D.NIIb_6548[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("N2_6583e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNIIc=float(line_ext_error[index_extra_error])
                        flux2D_error.NIIc_6584[i,j]=flux2D_error.NIIc_6584[i,j]+flux2D.NIIc_6584[i,j]*errextraNIIc
                    else:
                        errextraNIIc=float(line_ext_error[index_extra_error])
                        flux2D_error.NIIc_6584[i,j]=flux2D.NIIc_6584[i,j]*errextraNIIc
                error=nfes.flux_Error(flux2D.NIIc_6584[i,j],flux2D_error.NIIc_6584[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "N2_6584A", '{:03.3f}'.format(flux2D.NIIc_6584[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "N2_6584Ae", '{:03.3f}'.format(error)
            
            
                index_extra_error=line_names.index("O1_5577e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIa=float(line_ext_error[index_extra_error])
                        flux2D_error.OIa_5577[i,j]=flux2D_error.OIa_5577[i,j]+flux2D.OIa_5577[i,j]*errextraOIa
                    else:
                        errextraOIa=float(line_ext_error[index_extra_error])
                        flux2D_error.OIa_5577[i,j]=flux2D.OIa_5577[i,j]*errextraOIa
                error=nfes.flux_Error(flux2D.OIa_5577[i,j],flux2D_error.OIa_5577[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O1_5577A", '{:03.3f}'.format(flux2D.OIa_5577[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O1_5577Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("O1_6300e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIb=float(line_ext_error[index_extra_error])
                        flux2D_error.OIb_6300[i,j]=flux2D_error.OIb_6300[i,j]+flux2D.OIb_6300[i,j]*errextraOIb
                    else:
                        errextraOIb=float(line_ext_error[index_extra_error])
                        flux2D_error.OIb_6300[i,j]=flux2D.OIb_6300[i,j]*errextraOIb
                error=nfes.flux_Error(flux2D.OIb_6300[i,j],flux2D_error.OIb_6300[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O1_6300A", '{:03.3f}'.format(flux2D.OIb_6300[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O1_6300Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("O1_6363e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIc=float(line_ext_error[index_extra_error])
                        flux2D_error.OIc_6363[i,j]=flux2D_error.OIc_6363[i,j]+flux2D.OIc_6363[i,j]*errextraOIc
                    else:
                        errextraOIc=float(line_ext_error[index_extra_error])
                        flux2D_error.OIc_6363[i,j]=flux2D.OIc_6363[i,j]*errextraOIc
                error=nfes.flux_Error(flux2D.OIc_6363[i,j],flux2D_error.OIc_6363[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O1_6364A", '{:03.3f}'.format(flux2D.OIc_6363[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O1_6364Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("O2_3727e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIa_3727[i,j]=flux2D_error.OIIa_3727[i,j]+flux2D.OIIa_3727[i,j]*errextraOIIa
                    else:
                        errextraOIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIa_3727[i,j]=flux2D.OIIa_3727[i,j]*errextraOIIa
                error=nfes.flux_Error(flux2D.OIIa_3727[i,j],flux2D_error.OIIa_3727[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O2_3726A", '{:03.3f}'.format(flux2D.OIIa_3727[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O2_3726Ae", '{:03.3f}'.format(error)
            
                index_extra_error=line_names.index("O2_3729e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIb_3729[i,j]=flux2D_error.OIIb_3729[i,j]+flux2D.OIIb_3729[i,j]*errextraOIIb
                    else:
                        errextraOIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIb_3729[i,j]=flux2D.OIIb_3729[i,j]*errextraOIIb
                error=nfes.flux_Error(flux2D.OIIb_3729[i,j],flux2D_error.OIIb_3729[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O2_3729A", '{:03.3f}'.format(flux2D.OIIb_3729[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O2_3729Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("O2_7320e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIc=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIc_7320[i,j]=flux2D_error.OIIc_7320[i,j]+flux2D.OIIc_7320[i,j]*errextraOIIc
                    else:
                        errextraOIIc=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIc_7320[i,j]=flux2D.OIIc_7320[i,j]*errextraOIIc
                        
                error=nfes.flux_Error(flux2D.OIIc_7320[i,j],flux2D_error.OIIc_7320[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O2_7319A+", '{:03.3f}'.format(flux2D.OIIc_7320[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O2_7319A+e", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("O2_7330e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIId=float(line_ext_error[index_extra_error])
                        flux2D_error.OIId_7330[i,j]=flux2D_error.OIId_7330[i,j]+flux2D.OIId_7330[i,j]*errextraOIId
                    else:
                        errextraOIId=float(line_ext_error[index_extra_error])
                        flux2D_error.OIId_7330[i,j]=flux2D.OIId_7330[i,j]*errextraOIId
                        
                error=nfes.flux_Error(flux2D.OIId_7330[i,j],flux2D_error.OIId_7330[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O2_7330A+", '{:03.3f}'.format(flux2D.OIId_7330[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O2_7330A+e", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("O3_4363e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIIa_4363[i,j]=flux2D_error.OIIIa_4363[i,j]+flux2D.OIIIa_4363[i,j]*errextraOIIIa
                    else:
                        errextraOIIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIIa_4363[i,j]=flux2D.OIIIa_4363[i,j]*errextraOIIIa
                error=nfes.flux_Error(flux2D.OIIIa_4363[i,j],flux2D_error.OIIIa_4363[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O3_4363A", '{:03.3f}'.format(flux2D.OIIIa_4363[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O3_4363Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("O3_4959e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIIb_4959[i,j]=flux2D_error.OIIIb_4959[i,j]+flux2D.OIIIb_4959[i,j]*errextraOIIIb
                    else:
                        errextraOIIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIIb_4959[i,j]=flux2D.OIIIb_4959[i,j]*errextraOIIIb
                error=nfes.flux_Error(flux2D.OIIIb_4959[i,j],flux2D_error.OIIIb_4959[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O3_4959A", '{:03.3f}'.format(flux2D.OIIIb_4959[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O3_4959Ae", '{:03.3f}'.format(error)
        
                index_extra_error=line_names.index("O3_5007e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraOIIIc=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIIc_5007[i,j]=flux2D_error.OIIIc_5007[i,j]+flux2D.OIIIc_5007[i,j]*errextraOIIIc
                    else:
                        errextraOIIIc=float(line_ext_error[index_extra_error])
                        flux2D_error.OIIIc_5007[i,j]=flux2D.OIIIc_5007[i,j]*errextraOIIIc
                error=nfes.flux_Error(flux2D.OIIIc_5007[i,j],flux2D_error.OIIIc_5007[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "O3_5007A", '{:03.3f}'.format(flux2D.OIIIc_5007[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "O3_5007Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("Cl3_5517e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraClIIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.ClIIIa_5517[i,j]=flux2D_error.ClIIIa_5517[i,j]+flux2D.ClIIIa_5517[i,j]*errextraClIIIa
                    else:
                        errextraClIIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.ClIIIa_5517[i,j]=flux2D.ClIIIa_5517[i,j]*errextraClIIIa
                error=nfes.flux_Error(flux2D.ClIIIa_5517[i,j],flux2D_error.ClIIIa_5517[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])        
                print >>file3, "Cl3_5518A", '{:03.3f}'.format(flux2D.ClIIIa_5517[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "Cl3_5518Ae", '{:03.3f}'.format(error)
        
                index_extra_error=line_names.index("Cl3_5538e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraClIIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.ClIIIb_5538[i,j]=flux2D_error.ClIIIb_5538[i,j]+flux2D.ClIIIb_5538[i,j]*errextraClIIIb
                    else:
                        errextraClIIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.ClIIIb_5538[i,j]=flux2D.ClIIIb_5538[i,j]*errextraClIIIb
                error=nfes.flux_Error(flux2D.ClIIIb_5538[i,j],flux2D_error.ClIIIb_5538[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "Cl3_5538A", '{:03.3f}'.format(flux2D.ClIIIb_5538[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)        
#               print >>file3, "Cl3_5538Ae", '{:03.3f}'.format(error)        
                
                index_extra_error=line_names.index("Ne3_3868e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNeIIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.NeIIIa_3868[i,j]=flux2D_error.NeIIIa_3868[i,j]+flux2D.NeIIIa_3868[i,j]*errextraNeIIIa
                    else:
                        errextraNeIIIa=float(line_ext_error[index_extra_error])
                        flux2D_error.NeIIIa_3868[i,j]=flux2D.NeIIIa_3868[i,j]*errextraNeIIIa
                error=nfes.flux_Error(flux2D.NeIIIa_3868[i,j],flux2D_error.NeIIIa_3868[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "Ne3_3869A", '{:03.3f}'.format(flux2D.NeIIIa_3868[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "Ne3_3869Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("Ne3_3967e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraNeIIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.NeIIIb_3967[i,j]=flux2D_error.NeIIIb_3967[i,j]+flux2D.NeIIIb_3967[i,j]*errextraNeIIIb
                    else:
                        errextraNeIIIb=float(line_ext_error[index_extra_error])
                        flux2D_error.NeIIIb_3967[i,j]=flux2D.NeIIIb_3967[i,j]*errextraNeIIIb
                error=nfes.flux_Error(flux2D.NeIIIb_3967[i,j],flux2D_error.NeIIIb_3967[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "Ne3_3968A", '{:03.3f}'.format(flux2D.NeIIIb_3967[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)                
#               print >>file3, "Ne3_3968Ae", '{:03.3f}'.format(error)                
                
                index_extra_error=line_names.index("Ar3_7136e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIII=float(line_ext_error[index_extra_error])
                        flux2D_error.ArIII_7136[i,j]=flux2D_error.ArIII_7136[i,j]+flux2D.ArIII_7136[i,j]*errextraArIII
                    else:
                        errextraArIII=float(line_ext_error[index_extra_error])
                        flux2D_error.ArIII_7136[i,j]=flux2D.ArIII_7136[i,j]*errextraArIII
                error=nfes.flux_Error(flux2D.ArIII_7136[i,j],flux2D_error.ArIII_7136[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "Ar3_7136A", '{:03.3f}'.format(flux2D.ArIII_7136[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "Ar3_7136Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("Ar4_4712e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIVa=float(line_ext_error[index_extra_error])
                        flux2D_error.ArIVa_4712[i,j]=flux2D_error.ArIVa_4712[i,j]+flux2D.ArIVa_4712[i,j]*errextraArIVa
                    else:
                        errextraArIVa=float(line_ext_error[index_extra_error])
                        flux2D_error.ArIVa_4712[i,j]=flux2D.ArIVa_4712[i,j]*errextraArIVa
                error=nfes.flux_Error(flux2D.ArIVa_4712[i,j],flux2D_error.ArIVa_4712[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "Ar4_4711A", '{:03.3f}'.format(flux2D.ArIVa_4712[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "Ar4_4711Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("Ar4_4740e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraArIVb=float(line_ext_error[index_extra_error])
                        flux2D_error.ArIVb_4740[i,j]=flux2D_error.ArIVb_4740[i,j]+flux2D.ArIVb_4740[i,j]*errextraArIVb
                    else:
                        errextraArIVb=float(line_ext_error[index_extra_error])
                        flux2D_error.ArIVb_4740[i,j]=flux2D.ArIVb_4740[i,j]*errextraArIVb
                error=nfes.flux_Error(flux2D.ArIVb_4740[i,j],flux2D_error.ArIVb_4740[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "Ar4_4740A", '{:03.3f}'.format(flux2D.ArIVb_4740[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "Ar4_4740Ae", '{:03.3f}'.format(error)
                
                index_extra_error=line_names.index("C1_8727e")          
                if lines_available[index_extra_error]=="yes":
                    if line_ext_error[index_extra_error]!=0:
                        errextraCI=float(line_ext_error[index_extra_error])
                        flux2D_error.CI_8727[i,j]=flux2D_error.CI_8727[i,j]+flux2D.CI_8727[i,j]*errextraCI
                    else:
                        errextraCI=float(line_ext_error[index_extra_error])
                        flux2D_error.CI_8727[i,j]=flux2D.CI_8727[i,j]*errextraCI
                error=nfes.flux_Error(flux2D.CI_8727[i,j],flux2D_error.CI_8727[i,j],flux2D.Hb_4861[i,j],flux2D_error.Hb_4861[i,j])
                print >>file3, "C1_8728A", '{:03.3f}'.format(flux2D.CI_8727[i,j]*100/flux2D.Hb_4861[i,j]),'{:03.3f}'.format(error)
#               print >>file3, "C1_8728Ae", '{:03.3f}'.format(error)
#               print >>file3, "C2_????A", '{:03.3f}'.format(flux2D.CII_6461[i,j]*100/flux2D.Hb_4861[i,j])
                
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
                obs_fake.readData(obs_data, fileFormat='lines_in_rows_err_cols', errIsRelative=False)#, err_default=0.05)
                obs.readData(obs_data, fileFormat='lines_in_rows_err_cols', errIsRelative=False)#, err_default=0.05)
                
                index_MC=param_mod_name.index("number_fake_obs_MC")
                obs_fake.addMonteCarloObs(N = param_model_values[index_MC])

                
                obs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
                obs_fake.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)

                if obs.extinction.cHbeta[0]<2:
                    Chbeta_pyneb[i,j]=obs.extinction.cHbeta[0]
                    Chbeta_error_pyneb[i,j]=np.std(obs_fake.extinction.cHbeta)
                else:
                    Chbeta_pyneb[i,j]=0
                    Chbeta_error_pyneb[i,j]=0
    

                obs.extinction.law = extinction_law
                obs.extinction.R_V=RV_ext
                obs.correctData(normWave=4861.)
                
                obs_fake.extinction.law = extinction_law
                obs_fake.extinction.R_V=RV_ext
                obs_fake.correctData(normWave=4861.)



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
            

                flux_norm_2D.Ha_6563[i,j]=float(obs.getIntens()['H1r_6563A'])
                flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                par1=float(obs.getError()['H1r_6563A'])
                par2=float(obs.getError()['H1r_4861A'])
                par3=float(corHa_er/corHa)
                flux_norm_2D_error.Ha_6563[i,j]=flux_norm_2D.Ha_6563[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                flux_norm_2D_error.Hb_4861[i,j]=0.0

                index_extra_error=line_names.index("HI_4340e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.Hg_4340[i,j]=float(obs.getIntens()['H1r_4341A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['H1r_4341A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corHg_er/corHg)
                    flux_norm_2D_error.Hg_4340[i,j]=flux_norm_2D.Hg_4340[i,j]*np.sqrt(par1**2+par2**2+par3**2)

                index_extra_error=line_names.index("HI_4101e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.Hd_4101[i,j]=float(obs.getIntens()['H1r_4102A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['H1r_4102A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corHd_er/corHd)
                    flux_norm_2D_error.Hd_4101[i,j]=flux_norm_2D.Hd_4101[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("HeII_4686e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.HeIIa_4686[i,j]=float(obs.getIntens()['He2r_4686A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['He2r_4686A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corHeIIa_er/corHeIIa)
                    flux_norm_2D_error.HeIIa_4686[i,j]=flux_norm_2D.HeIIa_4686[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("HeII_5412e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.HeIIb_5412[i,j]=float(obs.getIntens()['He2r_5411A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['He2r_5411A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corHeIIb_er/corHeIIb)
                    flux_norm_2D_error.HeIIb_5412[i,j]=flux_norm_2D.HeIIb_5412[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("HeI_5876e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.HeIa_5876[i,j]=float(obs.getIntens()['He1r_5876A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['He1r_5876A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corHeIa_er/corHeIa)
                    flux_norm_2D_error.HeIa_5876[i,j]=flux_norm_2D.HeIa_5876[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                
                index_extra_error=line_names.index("HeI_6678e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.HeIb_6678[i,j]=float(obs.getIntens()['He1r_6678A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['He1r_6678A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corHeIb_er/corHeIb)
                    flux_norm_2D_error.HeIb_6678[i,j]=flux_norm_2D.HeIb_6678[i,j]*np.sqrt(par1**2+par2**2+par3**2)
    
                index_extra_error=line_names.index("N2_5755e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.NIIa_5755[i,j]=float(obs.getIntens()['N2_5755A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['N2_5755A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corNIIa_er/corNIIa)
                    flux_norm_2D_error.NIIa_5755[i,j]=flux_norm_2D.NIIa_5755[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("N2_6548e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.NIIb_6548[i,j]=float(obs.getIntens()['N2_6548A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['N2_6548A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corNIIb_er/corNIIb)
                    flux_norm_2D_error.NIIb_6548[i,j]=flux_norm_2D.NIIb_6548[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("N2_6583e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.NIIc_6584[i,j]=float(obs.getIntens()['N2_6584A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['N2_6584A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corNIIc_er/corNIIc)
                    flux_norm_2D_error.NIIc_6584[i,j]=flux_norm_2D.NIIc_6584[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("N1_5199e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.NI_5199[i,j]=float(obs.getIntens()['N1_5200A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['N1_5200A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corNI_er/corNI)
                    flux_norm_2D_error.NI_5199[i,j]=flux_norm_2D.NI_5199[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("O3_4363e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIIIa_4363[i,j]=float(obs.getIntens()['O3_4363A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O3_4363A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIIIa_er/corOIIIa)
                    flux_norm_2D_error.OIIIa_4363[i,j]=flux_norm_2D.OIIIa_4363[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("O3_4959e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIIIb_4959[i,j]=float(obs.getIntens()['O3_4959A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O3_4959A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIIIb_er/corOIIIb)
                    flux_norm_2D_error.OIIIb_4959[i,j]=flux_norm_2D.OIIIb_4959[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("O3_5007e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIIIc_5007[i,j]=float(obs.getIntens()['O3_5007A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O3_5007A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIIIc_er/corOIIIc)
                    flux_norm_2D_error.OIIIc_5007[i,j]=flux_norm_2D.OIIIc_5007[i,j]*np.sqrt(par1**2+par2**2+par3**2)
     
                index_extra_error=line_names.index("O2_3727e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIIa_3727[i,j]=float(obs.getIntens()['O2_3726A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O2_3726A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIIa_er/corOIIa)
                    flux_norm_2D_error.OIIa_3727[i,j]=flux_norm_2D.OIIa_3727[i,j]*np.sqrt(par1**2+par2**2+par3**2)

                index_extra_error=line_names.index("O2_3729e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIIb_3729[i,j]=float(obs.getIntens()['O2_3729A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O2_3729A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIIb_er/corOIIb)
                    flux_norm_2D_error.OIIb_3729[i,j]=flux_norm_2D.OIIb_3729[i,j]*np.sqrt(par1**2+par2**2+par3**2)
     
                index_extra_error=line_names.index("O2_7320e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIIc_7320[i,j]=float(obs.getIntens()['O2_7319A+'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O2_7319A+'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIIc_er/corOIIc)
                    flux_norm_2D_error.OIIc_7320[i,j]=flux_norm_2D.OIIc_7320[i,j]*np.sqrt(par1**2+par2**2+par3**2)
     
                index_extra_error=line_names.index("O2_7330e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIId_7330[i,j]=float(obs.getIntens()['O2_7330A+'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O2_7330A+'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIId_er/corOIId)
                    flux_norm_2D_error.OIId_7330[i,j]=flux_norm_2D.OIId_7330[i,j]*np.sqrt(par1**2+par2**2+par3**2)
        
                index_extra_error=line_names.index("O1_5577e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIa_5577[i,j]=float(obs.getIntens()['O1_5577A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O1_5577A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIa_er/corOIa)
                    flux_norm_2D_error.OIa_5577[i,j]=flux_norm_2D.OIa_5577[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("O1_6300e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIb_6300[i,j]=float(obs.getIntens()['O1_6300A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O1_6300A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIb_er/corOIb)
                    flux_norm_2D_error.OIb_6300[i,j]=flux_norm_2D.OIb_6300[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("O1_6363e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.OIc_6363[i,j]=float(obs.getIntens()['O1_6364A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['O1_6364A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corOIc_er/corOIc)
                    flux_norm_2D_error.OIc_6363[i,j]=flux_norm_2D.OIc_6363[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                
                index_extra_error=line_names.index("S2_6716e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.SIIa_6716[i,j]=float(obs.getIntens()['S2_6716A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['S2_6716A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corSIIa_er/corSIIa)
                    flux_norm_2D_error.SIIa_6716[i,j]=flux_norm_2D.SIIa_6716[i,j]*np.sqrt(par1**2+par2**2+par3**2)
            
                index_extra_error=line_names.index("S2_6731e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.SIIb_6731[i,j]=float(obs.getIntens()['S2_6731A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['S2_6731A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corSIIb_er/corSIIb)
                    flux_norm_2D_error.SIIb_6731[i,j]=flux_norm_2D.SIIb_6731[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("S3_6312e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.SIIIa_6312[i,j]=float(obs.getIntens()['S3_6312A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['S3_6312A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corSIIIa_er/corSIIIa)
                    flux_norm_2D_error.SIIIa_6312[i,j]=flux_norm_2D.SIIIa_6312[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("S3_9069e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.SIIIb_9069[i,j]=float(obs.getIntens()['S3_9069A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['S3_9069A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corSIIIb_er/corSIIIb)
                    flux_norm_2D_error.SIIIb_9069[i,j]=flux_norm_2D.SIIIb_9069[i,j]*np.sqrt(par1**2+par2**2+par3**2)
            
                index_extra_error=line_names.index("Cl3_5517e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.ClIIIa_5517[i,j]=float(obs.getIntens()['Cl3_5518A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['Cl3_5518A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corClIIIa_er/corClIIIa)
                    flux_norm_2D_error.ClIIIa_5517[i,j]=flux_norm_2D.ClIIIa_5517[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("Cl3_5538e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.ClIIIb_5538[i,j]=float(obs.getIntens()['Cl3_5538A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['Cl3_5538A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corClIIIb_er/corClIIIb)
                    flux_norm_2D_error.ClIIIb_5538[i,j]=flux_norm_2D.ClIIIb_5538[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("Ar3_7136e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.ArIII_7136[i,j]=float(obs.getIntens()['Ar3_7136A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['Ar3_7136A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corArIII_er/corArIII)
                    flux_norm_2D_error.ArIII_7136[i,j]=flux_norm_2D.ArIII_7136[i,j]*np.sqrt(par1**2+par2**2+par3**2)
                    
                index_extra_error=line_names.index("Ar4_4712e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.ArIVa_4712[i,j]=float(obs.getIntens()['Ar4_4711A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['Ar4_4711A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corArIVa_er/corArIVa)
                    flux_norm_2D_error.ArIVa_4712[i,j]=flux_norm_2D.ArIVa_4712[i,j]*np.sqrt(par1**2+par2**2+par3**2)
        
                index_extra_error=line_names.index("Ar4_4740e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.ArIVb_4740[i,j]=float(obs.getIntens()['Ar4_4740A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['Ar4_4740A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corArIVb_er/corArIVb)
                    flux_norm_2D_error.ArIVb_4740[i,j]=flux_norm_2D.ArIVb_4740[i,j]*np.sqrt(par1**2+par2**2+par3**2)

                index_extra_error=line_names.index("C1_8727e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.CI_8727[i,j]=float(obs.getIntens()['C1_8728A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['C1_8728A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corCI_er/corCI)
                    flux_norm_2D_error.CI_8727[i,j]=flux_norm_2D.CI_8727[i,j]*np.sqrt(par1**2+par2**2+par3**2)

                index_extra_error=line_names.index("C2_6461e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.CII_6461[i,j]=float(obs.getIntens()['H1r_6563A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['H1r_6563A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corCII_er/corCII)
                    flux_norm_2D_error.CII_6461[i,j]=flux_norm_2D.CII_6461[i,j]*np.sqrt(par1**2+par2**2+par3**2)

                index_extra_error=line_names.index("Ne3_3868e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.NeIIIa_3868[i,j]=float(obs.getIntens()['Ne3_3869A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['Ne3_3869A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corNeIIIa_er/corNeIIIa)
                    flux_norm_2D_error.NeIIIa_3868[i,j]=flux_norm_2D.NeIIIa_3868[i,j]*np.sqrt(par1**2+par2**2+par3**2)

                index_extra_error=line_names.index("Ne3_3967e")          
                if lines_available[index_extra_error]=="yes":
                    flux_norm_2D.NeIIIb_3967[i,j]=float(obs.getIntens()['Ne3_3968A'])
                    flux_norm_2D.Hb_4861[i,j]=float(obs.getIntens()['H1r_4861A'])
                    par1=float(obs.getError()['Ne3_3968A'])
                    par2=float(obs.getError()['H1r_4861A'])
                    par3=float(corNeIIIb_er/corNeIIIb)
                    flux_norm_2D_error.NeIIIb_3967[i,j]=flux_norm_2D.NeIIIb_3967[i,j]*np.sqrt(par1**2+par2**2+par3**2)


#########################################################################
### create 2D line ratio maps for all the line ratios defined by the user
#########################################################################

    g2Dlrs.lineratios(flux_norm_2D,line_names,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,hdr)



    maxHb=np.amax(flux_norm_2D.Hb_4861)
    maxper=95
    minper=5
    ypar,xpar,ang2=[],[],[]

    index_DD_name=DD_name.index("Ha/NII+_vs_Ha/SII+")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes":
        
#######################################################################################################
#SII 6717+6731/Ha VS NII 6548+6584/Ha
#######################################################################################################
        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]

        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if(flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.Ha_6563[i]/(flux_angles_norm.NIIb_6548[i]+flux_angles_norm.NIIc_6584[i]))
                xpar[number]=np.log10(flux_angles_norm.Ha_6563[i]/(flux_angles_norm.SIIa_6716[i]+flux_angles_norm.SIIb_6731[i]))
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1

                
        dataNII6548_6584=flux_norm_2D.NIIb_6548+flux_norm_2D.NIIc_6584
        dataSII6716_6731=flux_norm_2D.SIIa_6716+flux_norm_2D.SIIb_6731
    
#                         dds.NIISIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.NIISIIplot(dataNII6548_6584,flux_norm_2D.Ha_6563,dataSII6716_6731,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax1=plt.subplots(figsize=(13,13))

        ax1.add_artist(Ellipse(xy=(1.25, 0.30), width=1.0, height=3.1, linewidth=3,facecolor="none",edgecolor='sienna',alpha=1.0,angle=-45))
        ax1.text(1.8,0.8,"PNe",color='sienna',fontsize=24)

        ax1.add_artist(Rectangle(xy=(0.5, 0.2), width=0.4, height=0.45,linewidth=3,facecolor="none",edgecolor='g',alpha=1.0))
        ax1.text(0.53,0.5,"H II",color='g',fontsize=24)
        ax1.set_xlim(DD_xmin,DD_xmax)
        ax1.set_ylim(DD_ymin,DD_ymax)

        ax1.plot([-0.1, -0.1], [-0.5, 0.05], 'b-',linewidth=3)
        ax1.plot([0.15, 0.15], [-0.5, -0.05], 'b-',linewidth=3)
        ax1.plot([-0.1, 0.15], [-0.5, -0.5], 'b-',linewidth=3)
        ax1.plot([0.3, 0.3], [0,0.23], 'b-',linewidth=3)
        ax1.plot([-0.1, 0.3], [0.05,0.23], 'b-',linewidth=3)
        ax1.plot([0.15, 0.3], [-0.05,0], 'b-',linewidth=3)
        ax1.text(-0.08,-0.0,"SNR",color='b',fontsize=21, rotation=0)

        ax1.scatter(dataposx,dataposy,s=15,c="teal",marker=".")
        
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax1.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")
       
        ax1.set_xlabel(r'log(H$\alpha$/[S II] (6716+6731))',fontsize=21,labelpad=16)
        ax1.set_ylabel(r'log(H$\alpha$/[N II] (6548+6584))',fontsize=21,labelpad=16)

        ax1.tick_params(axis='x', labelsize=20, pad=12)
        ax1.tick_params(axis='y', labelsize=20, pad=12)
        
########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.Ha_6563[kk]/(flux_spec_slit_norm.NIIb_6548[kk]+flux_spec_slit_norm.NIIc_6584[kk]))
            spec_x[kk]=np.log10(flux_spec_slit_norm.Ha_6563[kk]/(flux_spec_slit_norm.SIIa_6716[kk]+flux_spec_slit_norm.SIIb_6731[kk]))

        ax1.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")
        
########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.59, 0.15, 0.28, 0.28]

            ax2 = fig.add_axes([left, bottom, width, height])
            diag=ax2.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            
            ax2.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax2.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax2.tick_params(axis='x', labelsize=18, pad=8)
            ax2.tick_params(axis='y', labelsize=18, pad=8)
            
        fig.savefig("output_Diagnostic_Diagrams/fig_NIIHa_SIIHa.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_NIIHa_SIIHa.pdf")
    
    
    
    
    index_DD_name=DD_name.index("OIII5007/Hb_vs_NII6584/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes":
        
#######################################################################################################
#NII 6584/Ha VS OIII 5007/Hb
#######################################################################################################
  
        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]

        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.NIIc_6584[i]/flux_angles_norm.Ha_6563[i]) 
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
                
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataNII6584=flux_norm_2D.NIIc_6584
    
#                         dds.OIIINIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIINIIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataNII6584,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
        t1 = np.arange(-3.0, 0.0, 0.01)   
        t2 = np.arange(-3.0, 0.3, 0.01)   

        t3 = np.arange(-1.5, -0.15, 0.01)   
        t4 = np.arange(-2.5, -1.5, 0.01)
        
        index_DD_BPT_Kau=DD_name.index("Kauffmann2003_BPT_NII")
        DD_BPT_Kau=DD_avail[index_DD_BPT_Kau]
        index_DD_BPT_Kew=DD_name.index("Kewley2001_BPT_NII")
        DD_BPT_Kew=DD_avail[index_DD_BPT_Kew]


        if (DD_BPT_Kau=="no" and DD_BPT_Kew=="no"):   
            ax3.plot([-3.0, 1.5], [1.3, 1.3], 'sienna',linewidth=3)
            ax3.text(0.25,0.80,"PNe",color='sienna',fontsize=21, rotation=0) 

            ax3.plot([-0.5, -0.5], [-2.5, 0.3], 'b-',linewidth=3)
            ax3.plot([-0.5, -0.1], [0.3, 2.0], 'b-',linewidth=3)
            ax3.text(-0.25,0.5,"SNR",color='b',fontsize=21, rotation=0) 
            
            ax3.plot(t3,-2.19*t3**6-15.61*t3**5-44.8*t3**4-66.3*t3**3-54.1*t3**2-24.23*t3-4.23,'g',linewidth=3)  #similar to Sabin et al. 2010
            ax3.plot(t4,0.945*t4/t4,'g',linewidth=3)  #similar to Sabin et al. 2010
            ax3.text(-1.5,0.85,"H II",color='g',fontsize=21, rotation=0) 
    
        if (DD_BPT_Kau=="yes"):   
            ax3.plot(t1,(0.61/(t1-0.05))+1.3,label="Kauffmann03")  #Kauf03
            ax3.text(-1.0,0.6,"H II",color='g',fontsize=21, rotation=0) 
        if (DD_BPT_Kew=="yes"):   
            ax3.plot(t2,(0.61/(t2-0.47))+1.19,label="Kewley01") #kew01
            ax3.text(-0.7,0.75,"AGN",color='g',fontsize=21, rotation=0) 
        ax3.legend()
    
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 

        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([N II] 6584/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)
        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/(flux_spec_slit_norm.Hb_4861[kk]))
            spec_x[kk]=np.log10(flux_spec_slit_norm.NIIc_6584[kk]/(flux_spec_slit_norm.Ha_6563[kk]))

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.22, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
            
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NIIHa.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NIIHa.pdf")
        
    
    
    index_DD_name=DD_name.index("OIII5007/Hb_vs_SII+/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes":
        
#######################################################################################################
#SII 6716+6731/Ha VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]

        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]): 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10((flux_angles_norm.SIIa_6716[i]+flux_angles_norm.SIIb_6731[i])/flux_angles_norm.Ha_6563[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataSII6716_6731=flux_norm_2D.SIIa_6716+flux_norm_2D.SIIb_6731
    
#                         dds.OIIISIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIISIIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataSII6716_6731,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)
        
        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
        t1 = np.arange(-3.0, 0.1, 0.01)   
        t2 = np.arange(-0.31, 0.75, 0.01)   
 
 
        index_DD_BPT_AGN=DD_name.index("main_AGN_line_BPT_SII")
        DD_BPT_AGN=DD_avail[index_DD_BPT_AGN]
        index_DD_BPT_LIN=DD_name.index("LINER/Sy2_line_BPT_SII")
        DD_BPT_LIN=DD_avail[index_DD_BPT_LIN]
 
        t3 = np.arange(-1.8, -0.4, 0.01)   
        t4 = np.arange(-3.0, -1.8, 0.01)
        if (DD_BPT_AGN=="no" and DD_BPT_LIN=="no"): 
            ax3.plot([-3.0, -0.4], [1.3, 1.3],'sienna',linewidth=3)
            ax3.plot([-0.4, -0.25], [1.3, 1.15], 'sienna',linewidth=3)
            ax3.plot([-0.25, 0.0], [1.15, 0.9], 'sienna',linewidth=3)
            ax3.plot([0.0, 0.0], [0.9, 0.0], 'sienna',linewidth=3)
            ax3.text(-0.5,0.8,"PNe",color='sienna',fontsize=21, rotation=0) 

            ax3.plot([-0.4, -0.4], [-2.5, 0.35], 'b-',linewidth=3)
            ax3.plot([-0.4, 0.1], [0.35, 2.0], 'b-',linewidth=3)
            ax3.text(0.15,0.25,"SNR",color='b',fontsize=21, rotation=0) 
        
            ax3.plot(t3,-0.31*t3*t3-1.12*t3-0.05,'g',linewidth=3)  #similar to Sabin et al. 2010
            ax3.plot(t4,0.96*t4/t4,'g',linewidth=3)  #similar to Sabin et al. 2010
            ax3.text(-0.8,0.3,"H II",color='g',fontsize=21, rotation=0) 

        if (DD_BPT_AGN=="yes"): 
            ax3.plot(t1,(0.72/(t1-0.32))+1.3,label="Kauffmann03")  #Kauf03
            ax3.text(-1.5,-0.0,"H II",color='g',fontsize=21, rotation=0) 
        if (DD_BPT_LIN=="yes"):
            ax3.plot(t2,(1.89*(t2))+0.76,label="Kewley01") #kew01  
            ax3.text(-0.4,1.25,"Seyfert",color='g',fontsize=21, rotation=0) 
            ax3.text(0.15,0.0,"LINER",color='g',fontsize=21, rotation=0) 
        
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([S II] 6716+6731/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/(flux_spec_slit_norm.Hb_4861[kk]))
            spec_x[kk]=np.log10((flux_spec_slit_norm.SIIa_6716[kk]+flux_spec_slit_norm.SIIb_6731[kk])/flux_spec_slit_norm.Ha_6563[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.22, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
            
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_SIIHa.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_SIIHa.pdf")
    
    
    index_DD_name=DD_name.index("OIII5007/Hb_vs_HeII4686/Hb")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes":
        
#######################################################################################################
#HeII 4686/Hb VS OIII 5007/Hb
#######################################################################################################
  
        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if(flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]): 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.HeIIa_4686[i]/flux_angles_norm.Hb_4861[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
                
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataHeII4686=flux_norm_2D.HeIIa_4686
    
#                         dds.OIIIHeII4686plot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIHeII4686plot(dataOIII5007,flux_norm_2D.Hb_4861,dataHeII4686,flux_norm_2D.Hb_4861,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log(He II 4686/H$\beta$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)
        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.HeIIa_4686[kk]/flux_spec_slit_norm.Hb_4861[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_HeII4686Hb.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_HeII4686Hb.pdf")
    
    
    index_DD_name=DD_name.index("OIII5007/Hb_vs_HeII5412/Hb")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#HeII 5412/Hb VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]

        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.HeIIb_5412[i]/flux_angles_norm.Hb_4861[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataHeII5412=flux_norm_2D.HeIIb_5412

#                         dds.NIISIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIHeII5412plot(dataOIII5007,flux_norm_2D.Hb_4861,dataHeII5412,flux_norm_2D.Hb_4861,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log(He II 5412/H$\beta$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)
        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.HeIIb_5412[kk]/flux_spec_slit_norm.Hb_4861[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
            
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_HeII5412Hb.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_HeII5412Hb.pdf")

       
    index_DD_name=DD_name.index("OIII5007/Hb_vs_NI5199/Hb")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#NI 5199/Hb VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.NI_5199[i]/flux_angles_norm.Hb_4861[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataNI5199=flux_norm_2D.NI_5199

#                         dds.OIIINIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIINIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataNI5199,flux_norm_2D.Hb_4861,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 

        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([N I] 5199/H$\beta$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)
        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.NI_5199[kk]/flux_spec_slit_norm.Hb_4861[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NI5199Hb.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NI5199Hb.pdf")



    index_DD_name=DD_name.index("OIII5007/Hb_vs_OI6300/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#OI 6300/Ha VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.OIb_6300[i]/flux_angles_norm.Ha_6563[i])
                number=number+1
            else: 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataOI6300=flux_norm_2D.OIb_6300

#                         dds.OIIIOIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIOIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataOI6300,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
        t1 = np.arange(-2.5, 0.0, 0.01)   
        t2 = np.arange(-0.53, 0.0, 0.01)  
 
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([O I] 6300/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)


        index_DD_BPT_AGN=DD_name.index("main_AGN_line_BPT_OI")
        DD_BPT_AGN=DD_avail[index_DD_BPT_AGN]
        index_DD_BPT_LIN=DD_name.index("LINER/Sy2_line_BPT_OI")
        DD_BPT_LIN=DD_avail[index_DD_BPT_LIN]
    
        if (DD_BPT_AGN=="yes"):   
            ax3.plot(t1,(0.73/(t1-0.59))+1.33,label="Kauffmann03")  #Kauf03
            ax3.text(-1.5,0.75,"HII",color='g',fontsize=21, rotation=0) 
        if (DD_BPT_LIN=="yes"):   
            ax3.plot(t2,(1.18*(t2))+1.30,label="Kewley01") #kew01
            ax3.text(-0.7,1.15,"Seyfert",color='g',fontsize=21, rotation=0) 
            ax3.text(-0.35,0.75,"LINER",color='g',fontsize=21, rotation=0) 


        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)


########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.OIb_6300[kk]/flux_spec_slit_norm.Ha_6563[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OI6300Ha.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OI6300Ha.pdf")

    
    
    
    index_DD_name=DD_name.index("OIII5007/Hb_vs_NII5755/Hb")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#NII 5755/Hb VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.NIIa_5755[i]/flux_angles_norm.Hb_4861[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataNII5755=flux_norm_2D.NIIa_5755

#                         dds.OIIINII5755plot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIINII5755plot(dataOIII5007,flux_norm_2D.Hb_4861,dataNII5755,flux_norm_2D.Hb_4861,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([N II] 5755/H$\beta$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)


########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.NIIa_5755[kk]/flux_spec_slit_norm.Hb_4861[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NII5755Hb.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NII5755Hb.pdf")

    
    index_DD_name=DD_name.index("SII6716/6731_vs_Ha/SII+")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
#######################################################################################################
#Ha/SII+ VS SII 6716/6731
#######################################################################################################
        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=flux_angles_norm.SIIa_6716[i]/flux_angles_norm.SIIb_6731[i]
                xpar[number]=np.log10(flux_angles_norm.Ha_6563[i]/(flux_angles_norm.SIIa_6716[i]+flux_angles_norm.SIIb_6731[i]))
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataSII6716=flux_norm_2D.SIIa_6716
        dataSII6731=flux_norm_2D.SIIb_6731
        dataSII=flux_norm_2D.SIIa_6716+flux_norm_2D.SIIb_6731

#                         dds.SII67166731HaSIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.SII67166731HaSIIplot(dataSII6716,dataSII6731,flux_norm_2D.Ha_6563,dataSII,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log(H$\alpha$/[S II] 6716+6731)',fontsize=21,labelpad=16)
        plt.ylabel(r'[S II] 6716/6731',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=(flux_spec_slit_norm.SIIa_6716[kk]/flux_spec_slit_norm.SIIb_6731[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.Ha_6563[kk]/(flux_spec_slit_norm.SIIa_6716[kk]+flux_spec_slit_norm.SIIb_6731[kk]))

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")
            
            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_SII_HaSII.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_SII_HaSII.pdf")

     
     
    index_DD_name=DD_name.index("SII6716/6731_vs_Ha/NII+")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#Ha/NII+ VS SII 6716/6731
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=flux_angles_norm.SIIa_6716[i]/flux_angles_norm.SIIb_6731[i]
                xpar[number]=np.log10(flux_angles_norm.Ha_6563[i]/(flux_angles_norm.NIIb_6548[i]+flux_angles_norm.NIIc_6584[i]))
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataSII6716=flux_norm_2D.SIIa_6716
        dataSII6731=flux_norm_2D.SIIb_6731
        dataNII=flux_norm_2D.NIIb_6548+flux_norm_2D.NIIc_6584

#                         dds.SII67166731HaNIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.SII67166731HaNIIplot(dataSII6716,dataSII6731,flux_norm_2D.Ha_6563,dataNII,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log(H$\alpha$/[N II] 6548+6584)',fontsize=21,labelpad=16)
        plt.ylabel(r'[S II] 6716/6731',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)
                
########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=flux_spec_slit_norm.SIIa_6716[kk]/flux_spec_slit_norm.SIIb_6731[kk]
            spec_x[kk]=np.log10(flux_spec_slit_norm.Ha_6563[kk]/(flux_spec_slit_norm.NIIb_6548[kk]+flux_spec_slit_norm.NIIc_6584[kk]))

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_SII_HaNII.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_SII_HaNII.pdf")


    
    index_DD_name=DD_name.index("NII+/SII+_vs_SII+/OI+")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
#######################################################################################################
#SII+/OI+ VS NII+ 6548+6584/SII+ 6716+6731
#######################################################################################################
        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):  
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10((flux_angles_norm.NIIb_6548[i]+flux_angles_norm.NIIc_6584[i])/(flux_angles_norm.SIIa_6716[i]+flux_angles_norm.SIIb_6731[i]))
                xpar[number]=np.log10((flux_angles_norm.SIIa_6716[i]+flux_angles_norm.SIIb_6731[i])/(flux_angles_norm.OIb_6300[i]))
                number=number+1
            else: 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataSII=flux_norm_2D.SIIa_6716+flux_norm_2D.SIIb_6731
        dataNII=flux_norm_2D.NIIb_6548+flux_norm_2D.NIIc_6584
        dataOI=flux_norm_2D.OIb_6300

#                         dds.NIISIISIIOIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.NIISIISIIOIplot(dataNII,dataSII,dataSII,dataOI,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([N II] 6548+6584/[S II] 6716+6731)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([S II] 6716+6731/[O I] 6300)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)
                
########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10((flux_spec_slit_norm.NIIb_6548[kk]+flux_spec_slit_norm.NIIc_6584[kk])/(flux_spec_slit_norm.SIIa_6716[kk]+flux_spec_slit_norm.SIIb_6731[kk]))
            spec_x[kk]=np.log10((flux_spec_slit_norm.SIIa_6716[kk]+flux_spec_slit_norm.SIIb_6731[kk])/(flux_spec_slit_norm.OIb_6300[kk]))

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_SIIOI_NIISII.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_SIIOI_NIISII.pdf")



    
    index_DD_name=DD_name.index("NII+/SII+_vs_NII+/OI+")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#NII+/OI+ VS NII+ 6548+6584/SII+ 6716+6731
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]): 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10((flux_angles_norm.NIIb_6548[i]+flux_angles_norm.NIIc_6584[i])/(flux_angles_norm.SIIa_6716[i]+flux_angles_norm.SIIb_6731[i]))
                xpar[number]=np.log10((flux_angles_norm.NIIb_6548[i]+flux_angles_norm.NIIc_6584[i])/(flux_angles_norm.OIb_6300[i]))
                number=number+1
            else: 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataSII=flux_norm_2D.SIIa_6716+flux_norm_2D.SIIb_6731
        dataNII=flux_norm_2D.NIIb_6548+flux_norm_2D.NIIc_6584
        dataOI=flux_norm_2D.OIb_6300

#                         dds.NIISIINIIOIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.NIISIINIIOIplot(dataNII,dataSII,dataNII,dataOI,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")


        plt.xlabel(r'log([N II] 6548+6584/[S II] 6716+6731)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([N II] 6548+6584/[O I] 6300)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10((flux_spec_slit_norm.NIIb_6548[kk]+flux_spec_slit_norm.NIIc_6584[kk])/(flux_spec_slit_norm.SIIa_6716[kk]+flux_spec_slit_norm.SIIb_6731[kk]))
            spec_x[kk]=np.log10((flux_spec_slit_norm.NIIb_6548[kk]+flux_spec_slit_norm.NIIc_6584[kk])/(flux_spec_slit_norm.OIb_6300[kk]))

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_NIIOI_NIISII.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_NIIOI_NIISII.pdf")


   
     
    index_DD_name=DD_name.index("OIII5007/Hb_vs_NII6584/OI6300")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
#######################################################################################################
#NII 6584/OI 6300 VS OIII 5007/Hb
#######################################################################################################
        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.NIIc_6584[i]/flux_angles_norm.OIb_6300[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataNII6584=flux_norm_2D.NIIc_6584
        dataOI6300=flux_norm_2D.OIb_6300
        
#                         dds.OIIINII6584OIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIINII6584OIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataNII6584,dataOI6300,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([N II] 6584/[O I] 6300)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.NIIc_6584[kk]/flux_spec_slit_norm.OIb_6300[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])

            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NII6584OI.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NII6584OI.pdf")

    
    index_DD_name=DD_name.index("OIII5007/Hb_vs_OII7320+/OIII5007")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#OII 7320+7330/OIII 5007 VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10((flux_angles_norm.OIIc_7320[i]+flux_angles_norm.OIId_7330[i])/flux_angles_norm.OIIIc_5007[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataOII732030=flux_norm_2D.OIIc_7320+flux_norm_2D.OIId_7330
        
#                         dds.OIIIOII732030OIIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIOII732030OIIIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataOII732030,dataOIII5007,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([O II] 7320+7330/[O III] 5007)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10((flux_spec_slit_norm.OIIc_7320[kk]+flux_spec_slit_norm.OIId_7330[kk])/flux_spec_slit_norm.OIIIc_5007[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OII732030OIII.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OII732030OIII.pdf")


    index_DD_name=DD_name.index("OIII5007/Hb_vs_OII3727+/OIII5007")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#OII 3727+3729/OIII 5007 VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]): 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10((flux_angles_norm.OIIa_3727[i]+flux_angles_norm.OIIb_3729[i])/flux_angles_norm.OIIIc_5007[i])
                number=number+1
            else: 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataOII372729=flux_norm_2D.OIIa_3727+flux_norm_2D.OIIb_3729
        
#                         dds.OIIIOII372729OIIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIOII372729OIIIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataOII372729,dataOIII5007,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([O II] 3727+3729/[O III] 5007)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10((flux_spec_slit_norm.OIIa_3727[kk]+flux_spec_slit_norm.OIIb_3729[kk])/flux_spec_slit_norm.OIIIc_5007[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OII372729OIII.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OII372729OIII.pdf")

    
    index_DD_name=DD_name.index("OIII5007/Hb_vs_OII3727+/Hb")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#OII 3727+3729/Hb VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10((flux_angles_norm.OIIa_3727[i]+flux_angles_norm.OIIb_3729[i])/flux_angles_norm.Hb_4861[i])
                number=number+1
            else: 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
                
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataOII372729=flux_norm_2D.OIIa_3727+flux_norm_2D.OIIb_3729
        
#                         dds.OIIIOII372729Hbplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIOII372729Hbplot(dataOIII5007,flux_norm_2D.Hb_4861,dataOII372729,flux_norm_2D.Hb_4861,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([O II] 3727+3729/H$\beta$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10((flux_spec_slit_norm.OIIa_3727[kk]+flux_spec_slit_norm.OIIb_3729[kk])/flux_spec_slit_norm.Hb_4861[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")
    
            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OII372729Hb.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OII372729Hb.pdf")

       
    index_DD_name=DD_name.index("OIII5007/Hb_vs_OII7320+/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#OII 3727+3729/Hb VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10((flux_angles_norm.OIIc_7320[i]+flux_angles_norm.OIId_7330[i])/flux_angles_norm.Ha_6563[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataOII732030=flux_norm_2D.OIIc_7320+flux_norm_2D.OIId_7330
        
#                         dds.OIIIOII732030Haplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIOII732030Haplot(dataOIII5007,flux_norm_2D.Hb_4861,dataOII732030,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([O II] 7320+7330/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10((flux_spec_slit_norm.OIIc_7320[kk]+flux_spec_slit_norm.OIId_7330[kk])/flux_spec_slit_norm.Ha_6563[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OII732030Ha.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OII732030Ha.pdf")
      
       
    index_DD_name=DD_name.index("OIII5007/Hb_vs_HeI5876/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#HeI 5876/Ha VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.HeIa_5876[i]/flux_angles_norm.Ha_6563[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataHeI5876=flux_norm_2D.HeIa_5876
      
#                         dds.OIIIHeIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIHeIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataHeI5876,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log(He I 5876/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.HeIa_5876[kk]/flux_spec_slit_norm.Ha_6563[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_HeIHa.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_HeIHa.pdf")

    
        
    index_DD_name=DD_name.index("OIII5007/Hb_vs_ArIII7136/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#Ar III 7136/Ha VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.ArIII_7136[i]/flux_angles_norm.Ha_6563[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataArIII7136=flux_norm_2D.ArIII_7136
        
#                         dds.OIIIArIIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIArIIIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataArIII7136,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([Ar III] 7136/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.ArIII_7136[kk]/flux_spec_slit_norm.Ha_6563[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_ArIIIHa.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_ArIIIHa.pdf")
    
        
    index_DD_name=DD_name.index("OIII5007/Hb_vs_OIII5007/OI6300")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#OIII 5007/OI 6300 VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.OIb_6300[i])
                number=number+1
            else: 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataOI6300=flux_norm_2D.OIb_6300
        
#                         dds.OIIIOIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIOIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataOIII5007,dataOI6300,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([O III] 5007/[O I] 6300)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.OIb_6300[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OIIIOI6300.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OIIIOI6300.pdf")
    
          
    index_DD_name=DD_name.index("OIII5007/Hb_vs_NI5199/NII6584")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#NI 5199/NII 6584 VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.NI_5199[i]/flux_angles_norm.NIIc_6584[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=100
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataNI5199=flux_norm_2D.NI_5199
        dataNII6584=flux_norm_2D.NIIc_6584

#                         dds.OIIININIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIININIIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataNI5199,dataNII6584,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([N I] 5199/[N II] 6584)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.NI_5199[kk]/flux_spec_slit_norm.NIIc_6584[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NINII.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NINII.pdf")
    
  
          
    index_DD_name=DD_name.index("OIII5007/Hb_vs_NeIII3869/Hb")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#NeIII3869/Hb VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.NeIIIa_3868[i]/flux_angles_norm.Hb_4861[i])
                number=number+1
            else: 
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataNeIII3868=flux_norm_2D.NeIIIa_3868

#                         dds.OIIINeIIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIINeIIIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataNeIII3868,flux_norm_2D.Hb_4861,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([Ne III] 3869/H$\beta$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.NeIIIa_3868[kk]/flux_spec_slit_norm.Hb_4861[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NeIIIHb.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_NeIIIHb.pdf")

 
          
    index_DD_name=DD_name.index("OIII5007/Hb_vs_OIII4363/Hg")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#OIII 4363/Hg VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.OIIIa_4363[i]/flux_angles_norm.Hg_4340[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
            
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataOIII4363=flux_norm_2D.OIIIa_4363

#                         dds.OIIIOIIIHgplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIOIIIHgplot(dataOIII5007,flux_norm_2D.Hb_4861,dataOIII4363,flux_norm_2D.Hg_4340,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([O III] 4363/H$\gamma$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)


########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.OIIIa_4363[kk]/flux_spec_slit_norm.Hg_4340[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OIIIHg.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_OIIIHg.pdf")


    
          
    index_DD_name=DD_name.index("OIII5007/Hb_vs_ArIV4712/4740")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#ArIV 4712/4740 VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=flux_angles_norm.ArIVa_4712[i]/flux_angles_norm.ArIVb_4740[i]
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataArIV4712=flux_norm_2D.ArIVa_4712
        dataArIV4740=flux_norm_2D.ArIVb_4740

#                         dds.OIIIArIVplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIIArIVplot(dataOIII5007,flux_norm_2D.Hb_4861,dataArIV4712,dataArIV4740,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'[Ar IV] 4712/4740',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)
        
########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=flux_spec_slit_norm.ArIVa_4712[kk]/flux_spec_slit_norm.ArIVb_4740[kk]

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")
    
########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################    
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_ArIV.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_ArIV.pdf")

          
    index_DD_name=DD_name.index("OIII5007/Hb_vs_CI8727/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#CI 8727/Ha VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.CI_8727[i]/flux_angles_norm.Ha_6563[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataCI8727=flux_norm_2D.CI_8727

#                         dds.OIIICIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIICIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataCI8727,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([C I] 8727/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.CI_8727[kk]/flux_spec_slit_norm.Ha_6563[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")
    
            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_CIHa.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_CIHa.pdf")

    
             
    index_DD_name=DD_name.index("OIII5007/Hb_vs_CII6461/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#CII 6461/Ha VS OIII 5007/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIIIc_5007[i]/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.CII_6461[i]/flux_angles_norm.Ha_6563[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
            
        dataOIII5007=flux_norm_2D.OIIIc_5007
        dataCII6461=flux_norm_2D.CII_6461

#                         dds.OIIICIIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OIIICIIplot(dataOIII5007,flux_norm_2D.Hb_4861,dataCII6461,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([C II] 6461/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O III] 5007/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIIIc_5007[kk]/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.CII_6461[kk]/flux_spec_slit_norm.Ha_6563[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_CIIHa.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIIIHb_CIIHa.pdf")


    
             
    index_DD_name=DD_name.index("OI6300/Ha_vs_CI8727/Ha")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#OI 8727/Ha VS OI 6300/Ha
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10(flux_angles_norm.OIb_6300[i]/flux_angles_norm.Ha_6563[i])
                xpar[number]=np.log10(flux_angles_norm.CI_8727[i]/flux_angles_norm.Ha_6563[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataOI6300=flux_norm_2D.OIb_6300
        dataCI8727=flux_norm_2D.CI_8727

#                         dds.OICIplot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.OICIplot(dataOI6300,flux_norm_2D.Ha_6563,dataCI8727,flux_norm_2D.Ha_6563,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log([C I] 8727/H$\alpha$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([O I] 6300/H$\alpha$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10(flux_spec_slit_norm.OIb_6300[kk]/flux_spec_slit_norm.Ha_6563[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.CI_8727[kk]/flux_spec_slit_norm.Ha_6563[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_OIHa_CIHa.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_OIHa_CIHa.pdf")

          
    index_DD_name=DD_name.index("ArIV+/Hb_vs_HeII4686/Hb")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#He II 4686/Hb VS ArIV 4712+4740/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10((flux_angles_norm.ArIVa_4712[i]+flux_angles_norm.ArIVb_4740[i])/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.HeIIa_4686[i]/flux_angles_norm.Hb_4861[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataArIV471240=flux_norm_2D.ArIVa_4712+flux_norm_2D.ArIVb_4740
        dataHeII4686=flux_norm_2D.HeIIa_4686

#                         dds.ArIVHeII4686plot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.ArIVHeII4686plot(dataArIV471240,flux_norm_2D.Hb_4861,dataHeII4686,flux_norm_2D.Hb_4861,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log(He II 4686/H$\beta$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([Ar IV] 4712+4740/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10((flux_spec_slit_norm.ArIVa_4712[kk]+flux_spec_slit_norm.ArIVb_4740[kk])/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.HeIIa_4686[kk]/flux_spec_slit_norm.Hb_4861[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_ArIVHb_HeII4686Hb.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_ArIVHb_HeII4686Hb.pdf")


        
    index_DD_name=DD_name.index("ArIV+/Hb_vs_HeII5412/Hb")
    DD_param=DD_avail[index_DD_name]
    if DD_param=="yes": 
        
#######################################################################################################
#He II 5412/Hb VS ArIV 4712+4740/Hb
#######################################################################################################

        DD_xmin=float(DDxmin[index_DD_name])
        DD_xmax=float(DDxmax[index_DD_name])
        DD_ymin=float(DDymin[index_DD_name])
        DD_ymax=float(DDymax[index_DD_name])
        ypar,xpar,ang2=[],[],[]
        
        number=0
        for i in range(0,len(flux_angles_norm.Ha_6563)):
            if (flux_angles_norm.Ha_6563[i]>0 and flux_angles_norm.Ha_6563[i]>2.84*flux_angles_norm.Hb_4861[i]):
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ang2[number]=ang[i]
                ypar[number]=np.log10((flux_angles_norm.ArIVa_4712[i]+flux_angles_norm.ArIVb_4740[i])/flux_angles_norm.Hb_4861[i])
                xpar[number]=np.log10(flux_angles_norm.HeIIb_5412[i]/flux_angles_norm.Hb_4861[i])
                number=number+1
            else:
                ypar.append(-99)
                xpar.append(-99)
                ang2.append(-99)
                ypar[number]=1000
                xpar[number]=1000
                ang2[number]=ang[i]
                number=number+1
        
        dataArIV471240=flux_norm_2D.ArIVa_4712+flux_norm_2D.ArIVb_4740
        dataHeII5412=flux_norm_2D.HeIIb_5412

#                         dds.ArIVHeII4686plot(y_numer,y_denom,x_numer,x_denom,sizex,sizey)
        dataposx,dataposy=dds.ArIVHeII4686plot(dataArIV471240,flux_norm_2D.Hb_4861,dataHeII5412,flux_norm_2D.Hb_4861,sizex,sizey)
        dataposx_mean,dataposx_std=ceos.estimates_without_outliers1D(dataposx,number)
        dataposy_mean,dataposy_std=ceos.estimates_without_outliers1D(dataposy,number)

        cm2 = plt.cm.get_cmap('Greys')
        cm = plt.cm.get_cmap('Accent')
        fig,ax3=plt.subplots(figsize=(13,13))
 
        ax3.scatter(dataposx,dataposy,s=15,c="teal",marker=".") 
        
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            diag1=ax3.scatter(xpar,ypar,s=100,facecolor="m",edgecolor="m",marker="o")

        plt.xlabel(r'log(He II 5412/H$\beta$)',fontsize=21,labelpad=16)
        plt.ylabel(r'log([Ar IV] 4712+4740/H$\beta$)',fontsize=21,labelpad=16)

        ax3.set_xlim(DD_xmin,DD_xmax)
        ax3.set_ylim(DD_ymin,DD_ymax)
        
        ax3.tick_params(axis='x', labelsize=20, pad=12)
        ax3.tick_params(axis='y', labelsize=20, pad=12)

########################################################################################
#       plot the line ratio values from the specific slits and rotation analysis modules
########################################################################################

        spec_x,spec_y=[],[]
        for kk in range(0,len(flux_spec_slit_norm.Ha_6563)):
            spec_x.append(-99)
            spec_y.append(-99)
            spec_y[kk]=np.log10((flux_spec_slit_norm.ArIVa_4712[kk]+flux_spec_slit_norm.ArIVb_4740[kk])/flux_spec_slit_norm.Hb_4861[kk])
            spec_x[kk]=np.log10(flux_spec_slit_norm.HeIIb_5412[kk]/flux_spec_slit_norm.Hb_4861[kk])

        ax3.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

########################################################################################
#      define inner plot in the diagnostic diagram
########################################################################################
        index=param_estimated.index('angular_analysis_task')
        if param_requered[index]=="yes":
            left, bottom, width, height = [0.53, 0.15, 0.28, 0.28]
            ax4 = fig.add_axes([left, bottom, width, height])
            diag=ax4.scatter(xpar,ypar,s=100,c=ang,cmap=cm,marker="o")

            bounds=[0,45,90,135,180,225,270,315,360]
            cbar=plt.colorbar(diag, cmap=cm, boundaries=bounds, ticks=[0,45,90,135,180,225,270,315,360])
            ax4.scatter(np.mean(xpar),np.mean(ypar),s=300,c="OrangeRed",marker="x",linewidths=4) # average value from integrated spectrum
            ax4.scatter(spec_x,spec_y,s=100,c="orange",edgecolors="black",marker="d")

            ax4.tick_params(axis='x', labelsize=18, pad=8)
            ax4.tick_params(axis='y', labelsize=18, pad=8)
        
        fig.savefig("output_Diagnostic_Diagrams/fig_ArIVHb_HeII5412Hb.png")
        fig.savefig("output_Diagnostic_Diagrams/fig_ArIVHb_HeII5412Hb.pdf")



    return      

        
