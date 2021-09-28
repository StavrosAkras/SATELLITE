# TeNe_2D_script.py:
# Conputes c(Hb), Te and Ne parameters for each spaxel in 2D maps using the PyNeb package (Luridiana et al. 2015) as well as 
# the ionic/elemental abundances and ICFs. The routine also returns 2D maps in FITS format and histrograms for all these parameters.
# Finally, the results of 5%, 25%, 50%, 75%, 95% percentile, mean value and standard deviation
# for all the maps are saved in the ascii file (statisticnumbers_1.txt)
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
import os, fnmatch
from numpy import inf


from satellite import diagnotic_diagrams_script as dds
from satellite import save_FITSimages_script as sfs 
from satellite import ionicabundances_2D_script as ia2Ds
from satellite import element_abundances_ICFs_2D_script as eaI2Ds
from satellite import norm_flux_error_script as nfes
from satellite import histogram_script as hs
from satellite import statistics_calculations_script as scs
from satellite import mean_value_script as meanvs

def TeNe(flux2D,flux2D_error,line_names,line_ext_error,lines_available,param_mod_name,param_model_values,param_estimated,param_requered,hdr,Te_PA,Ne_PA):
    index_size=param_mod_name.index("total_num_pixels_horiz")
    sizex=param_model_values[index_size]
    sizey=param_model_values[index_size]
    


    class Te:
        NIISII=np.zeros((sizex, sizey))
        OISII=np.zeros((sizex, sizey))
        OIISII=np.zeros((sizex, sizey))
        OIIISII=np.zeros((sizex, sizey))
        SIIISII=np.zeros((sizex, sizey))
        OIIOII=np.zeros((sizex, sizey))
        NIIOII=np.zeros((sizex, sizey))
        OIOII=np.zeros((sizex, sizey))
        OIIIClIII=np.zeros((sizex, sizey))
        SIIIClIII=np.zeros((sizex, sizey))
        OIIIArVI=np.zeros((sizex, sizey))
        SIIIArVI=np.zeros((sizex, sizey))
        NIIClIII=np.zeros((sizex, sizey))
        
    class Ne:
        NIISII=np.zeros((sizex, sizey))
        OISII=np.zeros((sizex, sizey))
        OIISII=np.zeros((sizex, sizey))
        OIIISII=np.zeros((sizex, sizey))
        SIIISII=np.zeros((sizex, sizey))
        OIIOII=np.zeros((sizex, sizey))
        NIIOII=np.zeros((sizex, sizey))
        OIOII=np.zeros((sizex, sizey))
        OIIIClIII=np.zeros((sizex, sizey))
        SIIIClIII=np.zeros((sizex, sizey))
        OIIIArVI=np.zeros((sizex, sizey))
        SIIIArVI=np.zeros((sizex, sizey))
        NIIClIII=np.zeros((sizex, sizey))
        
    class Te_error:
        NIISII=np.zeros((sizex, sizey))
        OISII=np.zeros((sizex, sizey))
        OIISII=np.zeros((sizex, sizey))
        OIIISII=np.zeros((sizex, sizey))
        SIIISII=np.zeros((sizex, sizey))
        OIIOII=np.zeros((sizex, sizey))
        NIIOII=np.zeros((sizex, sizey))
        OIOII=np.zeros((sizex, sizey))
        OIIIClIII=np.zeros((sizex, sizey))
        SIIIClIII=np.zeros((sizex, sizey))
        OIIIArVI=np.zeros((sizex, sizey))
        SIIIArVI=np.zeros((sizex, sizey))
        NIIClIII=np.zeros((sizex, sizey))
        
    class Ne_error:
        NIISII=np.zeros((sizex, sizey))
        OISII=np.zeros((sizex, sizey))
        OIISII=np.zeros((sizex, sizey))
        OIIISII=np.zeros((sizex, sizey))
        SIIISII=np.zeros((sizex, sizey))
        OIIOII=np.zeros((sizex, sizey))
        NIIOII=np.zeros((sizex, sizey))
        OIOII=np.zeros((sizex, sizey))
        OIIIClIII=np.zeros((sizex, sizey))
        SIIIClIII=np.zeros((sizex, sizey))
        OIIIArVI=np.zeros((sizex, sizey))
        SIIIArVI=np.zeros((sizex, sizey))
        NIIClIII=np.zeros((sizex, sizey))
    
    class Te_fake:
        NIISII=np.zeros((sizex, sizey))
        OISII=np.zeros((sizex, sizey))
        OIISII=np.zeros((sizex, sizey))
        OIIISII=np.zeros((sizex, sizey))
        SIIISII=np.zeros((sizex, sizey))
        OIIOII=np.zeros((sizex, sizey))
        NIIOII=np.zeros((sizex, sizey))
        OIOII=np.zeros((sizex, sizey))
        OIIIClIII=np.zeros((sizex, sizey))
        SIIIClIII=np.zeros((sizex, sizey))
        OIIIArVI=np.zeros((sizex, sizey))
        SIIIArVI=np.zeros((sizex, sizey))
        NIIClIII=np.zeros((sizex, sizey))
        
    class Ne_fake:
        NIISII=np.zeros((sizex, sizey))
        OISII=np.zeros((sizex, sizey))
        OIISII=np.zeros((sizex, sizey))
        OIIISII=np.zeros((sizex, sizey))
        SIIISII=np.zeros((sizex, sizey))
        OIIOII=np.zeros((sizex, sizey))
        NIIOII=np.zeros((sizex, sizey))
        OIOII=np.zeros((sizex, sizey))
        OIIIClIII=np.zeros((sizex, sizey))
        SIIIClIII=np.zeros((sizex, sizey))
        OIIIArVI=np.zeros((sizex, sizey))
        SIIIArVI=np.zeros((sizex, sizey))
        NIIClIII=np.zeros((sizex, sizey))    
    
    class ion_abun:
        HeIa=np.zeros((sizex, sizey))
        HeIb=np.zeros((sizex, sizey))
        HeIIa=np.zeros((sizex, sizey))
        HeIIb=np.zeros((sizex, sizey))
        NI=np.zeros((sizex, sizey))
        NIIa=np.zeros((sizex, sizey))
        NIIb=np.zeros((sizex, sizey))
        NIIc=np.zeros((sizex, sizey))
        OIa=np.zeros((sizex, sizey))
        OIb=np.zeros((sizex, sizey))
        OIc=np.zeros((sizex, sizey))
        OIIa=np.zeros((sizex, sizey))
        OIIb=np.zeros((sizex, sizey))
        OIIc=np.zeros((sizex, sizey))
        OIId=np.zeros((sizex, sizey))
        OIIIa=np.zeros((sizex, sizey))
        OIIIb=np.zeros((sizex, sizey))
        OIIIc=np.zeros((sizex, sizey))
        SIIa=np.zeros((sizex, sizey))
        SIIb=np.zeros((sizex, sizey))
        SIIIa=np.zeros((sizex, sizey))
        SIIIb=np.zeros((sizex, sizey))
        NeIIIa=np.zeros((sizex, sizey))
        NeIIIb=np.zeros((sizex, sizey))
        ArIII=np.zeros((sizex, sizey))
        ArIVa=np.zeros((sizex, sizey))
        ArIVb=np.zeros((sizex, sizey))
        ClIIIa=np.zeros((sizex, sizey))
        ClIIIb=np.zeros((sizex, sizey))
    
    class ion_abun_error:
        HeIa=np.zeros((sizex, sizey))
        HeIb=np.zeros((sizex, sizey))
        HeIIa=np.zeros((sizex, sizey))
        HeIIb=np.zeros((sizex, sizey))
        NI=np.zeros((sizex, sizey))
        NIIa=np.zeros((sizex, sizey))
        NIIb=np.zeros((sizex, sizey))
        NIIc=np.zeros((sizex, sizey))
        OIa=np.zeros((sizex, sizey))
        OIb=np.zeros((sizex, sizey))
        OIc=np.zeros((sizex, sizey))
        OIIa=np.zeros((sizex, sizey))
        OIIb=np.zeros((sizex, sizey))
        OIIc=np.zeros((sizex, sizey))
        OIId=np.zeros((sizex, sizey))
        OIIIa=np.zeros((sizex, sizey))
        OIIIb=np.zeros((sizex, sizey))
        OIIIc=np.zeros((sizex, sizey))
        SIIa=np.zeros((sizex, sizey))
        SIIb=np.zeros((sizex, sizey))
        SIIIa=np.zeros((sizex, sizey))
        SIIIb=np.zeros((sizex, sizey))
        NeIIIa=np.zeros((sizex, sizey))
        NeIIIb=np.zeros((sizex, sizey))
        ArIII=np.zeros((sizex, sizey))
        ArIVa=np.zeros((sizex, sizey))
        ArIVb=np.zeros((sizex, sizey))
        ClIIIa=np.zeros((sizex, sizey))
        ClIIIb=np.zeros((sizex, sizey))
    
    class elem_abun_KB:
        He=np.zeros((sizex, sizey))
        N=np.zeros((sizex, sizey))
        O=np.zeros((sizex, sizey))
        S=np.zeros((sizex, sizey))
        Ne=np.zeros((sizex, sizey))
        Ar=np.zeros((sizex, sizey))
        Cl=np.zeros((sizex, sizey))
    
    class elem_abun_KB_ratio:
        NO=np.zeros((sizex, sizey))
        SO=np.zeros((sizex, sizey))
        NeO=np.zeros((sizex, sizey))
        ArO=np.zeros((sizex, sizey))
        ClO=np.zeros((sizex, sizey))
    
    class elem_ICF_KB:
        He=np.ones((sizex, sizey))
        N=np.ones((sizex, sizey))
        O=np.ones((sizex, sizey))
        S=np.ones((sizex, sizey))
        Ne=np.ones((sizex, sizey))
        Ar=np.ones((sizex, sizey))
        Cl=np.ones((sizex, sizey))
        
    class elem_abun_DIMS:
        He=np.zeros((sizex, sizey))
        N=np.zeros((sizex, sizey))
        O=np.zeros((sizex, sizey))
        S=np.zeros((sizex, sizey))
        Ne=np.zeros((sizex, sizey))
        Ar=np.zeros((sizex, sizey))
        Cl=np.zeros((sizex, sizey))
    
    class elem_abun_DIMS_ratio:
        NO=np.zeros((sizex, sizey))
        SO=np.zeros((sizex, sizey))
        NeO=np.zeros((sizex, sizey))
        ArO=np.zeros((sizex, sizey))
        ClO=np.zeros((sizex, sizey))
        
    class elem_ICF_DIMS:
        He=np.ones((sizex, sizey))
        N=np.ones((sizex, sizey))
        O=np.ones((sizex, sizey))
        S=np.ones((sizex, sizey))
        Ne=np.ones((sizex, sizey))
        Ar=np.ones((sizex, sizey))
        Cl=np.ones((sizex, sizey))
    
    class elem_abun_KB_error:
        He=np.zeros((sizex, sizey))
        N=np.zeros((sizex, sizey))
        O=np.zeros((sizex, sizey))
        S=np.zeros((sizex, sizey))
        Ne=np.zeros((sizex, sizey))
        Ar=np.zeros((sizex, sizey))
        Cl=np.zeros((sizex, sizey))
    
    class elem_abun_KB_ratio_error:
        NO=np.zeros((sizex, sizey))
        SO=np.zeros((sizex, sizey))
        NeO=np.zeros((sizex, sizey))
        ArO=np.zeros((sizex, sizey))
        ClO=np.zeros((sizex, sizey))
    
    class elem_ICF_KB_error:
        He=np.zeros((sizex, sizey))
        N=np.zeros((sizex, sizey))
        O=np.zeros((sizex, sizey))
        S=np.zeros((sizex, sizey))
        Ne=np.zeros((sizex, sizey))
        Ar=np.zeros((sizex, sizey))
        Cl=np.zeros((sizex, sizey))
        
    class elem_abun_DIMS_error:
        He=np.zeros((sizex, sizey))
        N=np.zeros((sizex, sizey))
        O=np.zeros((sizex, sizey))
        S=np.zeros((sizex, sizey))
        Ne=np.zeros((sizex, sizey))
        Ar=np.zeros((sizex, sizey))
        Cl=np.zeros((sizex, sizey))
    
    class elem_abun_DIMS_ratio_error:
        NO=np.zeros((sizex, sizey))
        SO=np.zeros((sizex, sizey))
        NeO=np.zeros((sizex, sizey))
        ArO=np.zeros((sizex, sizey))
        ClO=np.zeros((sizex, sizey))
    
    class elem_ICF_DIMS_error:
        He=np.zeros((sizex, sizey))
        N=np.zeros((sizex, sizey))
        O=np.zeros((sizex, sizey))
        S=np.zeros((sizex, sizey))
        Ne=np.zeros((sizex, sizey))
        Ar=np.zeros((sizex, sizey))
        Cl=np.zeros((sizex, sizey))
    
    ###################################################
    ### for all the angle in the rotation angle task###
    ###################################################

    pn.log_.open_file('my_log.txt')

    index_atomicData=12
    name_atomicdata=str(param_mod_name[index_atomicData])
    
    if (name_atomicdata!="Chianti"):
        if (name_atomicdata!="IRAF_09"):
            pn.atomicData.includeFitsPath()
            pn.atomicData.includeDeprecatedPath()
            pn.atomicData.setDataFileDict(name_atomicdata)
        else:
            pn.atomicData.includeFitsPath()
            pn.atomicData.setDataFileDict(name_atomicdata)
    else:
        pn.atomicData.addAllChianti()

    Chbeta_pyneb=np.zeros((sizex, sizey))
    Chbeta_error_pyneb=np.zeros((sizex, sizey))
    
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
                
                index_extra_error=line_names.index("Cl3_5517s")          
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
                # deredden data with fitspatric'law
                # instantiate the Diagnostics class
                diags = pn.Diagnostics()
                # include in diags the relevant line ratios

                ########################
                ### Density diagnostics
                ########################
                index_line1=line_names.index("S2_6716s")
                index_line2=line_names.index("S2_6731s")
                index_Ne_SII="no"
                if lines_available[index_line1]=="yes" and  lines_available[index_line2]=="yes" and flux2D.SIIa_6716[i,j]>0 and flux2D.SIIb_6731[i,j]>0:
                    diags.addDiag([
                        '[SII] 6731/6716',
                        ])    
                    diags.addClabel('[SII] 6731/6716', '[SII]a')
                    index_Ne_SII="yes"

                index_line1=line_names.index("O2_3727s")
                index_line2=line_names.index("O2_3729s")
                index_Ne_OII="no"
                if lines_available[index_line1]=="yes" and  lines_available[index_line2]=="yes"  and flux2D.OIIa_3727[i,j]>0 and flux2D.OIIb_3729[i,j]>0:
                    diags.addDiag([
                        '[OII] 3726/3729',
                        ])    
                    index_Ne_OII="yes"

                index_line1=line_names.index("Cl3_5517s")
                index_line2=line_names.index("Cl3_5538s")
                index_Ne_ClIII="no"
                if lines_available[index_line1]=="yes" and  lines_available[index_line2]=="yes" and flux2D.ClIIIa_5517[i,j]>0 and flux2D.ClIIIb_5538[i,j]>0:
                    diags.addDiag([
                        '[ClIII] 5538/5518',
                        ])    
                    index_Ne_ClIII="yes"

                index_line1=line_names.index("Ar4_4712s")
                index_line2=line_names.index("Ar4_4740s")
                index_Ne_ArVI="no"
                if lines_available[index_line1]=="yes" and  lines_available[index_line2]=="yes" and flux2D.ArIVa_4712[i,j]>0 and flux2D.ArIVb_4740[i,j]>0:
                    diags.addDiag([
                        '[ArIV] 4740/4711',
                        ])   
                    index_Ne_ArVI="yes"

                ###########################
                ### Temperature diagnostics
                ###########################

                index_line1=line_names.index("N2_5755s")
                index_line2=line_names.index("N2_6548s")
                index_line3=line_names.index("N2_6583s")
                index_Te_NII="no"
                if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes" and flux2D.NIIa_5755[i,j]>0 and flux2D.NIIb_6548[i,j]>0 and flux2D.NIIc_6584[i,j]>0:
                    diags.addDiag([
                        '[NII] 5755/6584+',
                        ])    
                    index_Te_NII="yes"

                index_line1=line_names.index("O3_4363s")
                index_line2=line_names.index("O3_4959s")
                index_line3=line_names.index("O3_5007s")
                index_Te_OIII="no"
                if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes" and flux2D.OIIa_3727[i,j]>0 and flux2D.OIIIb_4959[i,j]>0 and flux2D.OIIIc_5007[i,j]>0:
                    diags.addDiag([
                        '[OIII] 4363/5007+',
                        ])      
                    index_Te_OIII="yes"

                index_line1=line_names.index("O1_6300s")
                index_line2=line_names.index("O1_6363s")
                index_line3=line_names.index("O1_5577s")
                index_Te_OI="no"
                index_Te_OI2="no"
                if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes" and flux2D.OIa_5577[i,j]>0 and flux2D.OIb_6300[i,j]>0 and flux2D.OIc_6363[i,j]>0:
                    diags.addDiag([
                        '[OI] 5577/6300+',
                        ]) 
                    index_Te_OI="yes"
                elif lines_available[index_line1]=="yes" and lines_available[index_line2]=="no" and lines_available[index_line3]=="yes" and flux2D.OIa_5577[i,j]>0 and flux2D.OIb_6300[i,j]>0:
                    diags.addDiag([
                        '[OI] 5577/6300',
                        ])       
                    index_Te_OI2="yes"

                index_line1=line_names.index("O2_3727s")
                index_line2=line_names.index("O2_3729s")
                index_line3=line_names.index("O2_7320s")
                index_line4=line_names.index("O2_7330s")
                index_Te_OII="no"
                if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes" and lines_available[index_line4]=="yes" and flux2D.OIIa_3727[i,j]>0 and flux2D.OIIb_3729[i,j]>0 and flux2D.OIIc_7320[i,j]>0 and flux2D.OIId_7330[i,j]>0:
                    diags.addDiag([
                        '[OII] 3727+/7325+',
                        ])
                    index_Te_OII="yes"
    
                index_line1=line_names.index("S3_6312s")
                index_line2=line_names.index("S3_9069s")
                #index_line3=line_names.index("S3_9592s")
                index_Te_SIII="no"
                if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and flux2D.SIIIa_6312[i,j]>0 and flux2D.SIIIb_9069[i,j]>0:
                    diags.addDiag([
                        '[SIII] 6312/9069',
                        ])
                    index_Te_SIII="yes"
#               elif lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes":
#                   diags.addDiag([
#                        '[SIII] 6312/9200+',
#                       ])  
#                   index_Te_SIII="yes"

                print("######################################################################################")
                #The observed ratio can be automatically extracted from an Observation object named obs:
                if (index_Te_NII=="yes" and index_Ne_SII=="yes"):
                    a, b = diags.getCrossTemDen('[NII] 5755/6584+', '[SII] 6731/6716', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+', '[SII] 6731/6716', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.NIISII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.NIISII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.NIISII=a
                        Ne_fake.NIISII=b
                        
                        Te.NIISII[i,j]=a1
                        Ne.NIISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.NIISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.NIISII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.NIISII=a
                        Ne_fake.NIISII=b
                        
                        Te.NIISII[i,j]=a1
                        Ne.NIISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.NIISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.NIISII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: low zones (NII-SII): Te_NII = {0:5.0f} K, Ne_SII = {1:7.1f} cm-1'.format(Te.NIISII[i,j], Ne.NIISII[i,j]))
    
                if (index_Te_OI=="yes" and index_Te_OI2=="no" and index_Ne_SII=="yes"):
                    a, b = diags.getCrossTemDen('[OI] 5577/6300+', '[SII] 6731/6716', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OI] 5577/6300+', '[SII] 6731/6716', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OISII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OISII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:

                        Te_fake.OISII=a
                        Ne_fake.OISII=b
                        
                        Te.OISII[i,j]=a1
                        Ne.OISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OISII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OISII=a
                        Ne_fake.OISII=b
                        
                        Te.OISII[i,j]=a1
                        Ne.OISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OISII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: low zones (OI-SII): Te_OI = {0:5.0f} K, Ne_SII = {1:7.1f} cm-1'.format(Te.OISII[i,j], Ne.OISII[i,j]))
                
                if (index_Te_OI=="no" and index_Te_OI2=="yes" and index_Ne_SII=="yes"):
                    a, b = diags.getCrossTemDen('[OI] 5577/6300', '[SII] 6731/6716', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OI] 5577/6300', '[SII] 6731/6716', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OISII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OISII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:

                        Te_fake.OISII=a
                        Ne_fake.OISII=b
                        
                        Te.OISII[i,j]=a1
                        Ne.OISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OISII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OISII=a
                        Ne_fake.OISII=b
                        
                        Te.OISII[i,j]=a1
                        Ne.OISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OISII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: low zones (OI-SII): Te_OI = {0:5.0f} K, Ne_SII = {1:7.1f} cm-1'.format(Te.OISII[i,j], Ne.OISII[i,j]))
                
                if (index_Te_OII=="yes" and index_Ne_SII=="yes"):
                    a, b = diags.getCrossTemDen('[OII] 3727+/7325+', '[SII] 6731/6716', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OII] 3727+/7325+', '[SII] 6731/6716', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OIISII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OIISII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.OIISII=a
                        Ne_fake.OIISII=b
                        
                        Te.OIISII[i,j]=a1
                        Ne.OIISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIISII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OIISII=a
                        Ne_fake.OIISII=b
                        
                        Te.OIISII[i,j]=a1
                        Ne.OIISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIISII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: low zones (OII-SII): Te_OII = {0:5.0f} K, Ne_SII = {1:7.1f} cm-1'.format(Te.OIISII[i,j], Ne.OIISII[i,j]))
    
                if (index_Te_OIII=="yes" and index_Ne_SII=="yes"):
                    a, b = diags.getCrossTemDen('[OIII] 4363/5007+', '[SII] 6731/6716', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+', '[SII] 6731/6716', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OIIISII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OIIISII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.OIIISII=a
                        Ne_fake.OIIISII=b
                        
                        Te.OIIISII[i,j]=a1
                        Ne.OIIISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIIISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIIISII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OIIISII=a
                        Ne_fake.OIIISII=b
                        
                        Te.OIIISII[i,j]=a1
                        Ne.OIIISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIIISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIIISII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: medium zones (OIII-SII): Te_OIII = {0:5.0f} K, Ne_SII = {1:7.1f} cm-1'.format(Te.OIIISII[i,j], Ne.OIIISII[i,j]))

                if (index_Te_SIII=="yes" and index_Ne_SII=="yes"):
                    a, b = diags.getCrossTemDen('[SIII] 6312/9069', '[SII] 6731/6716', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069', '[SII] 6731/6716', obs=obs)

                    mean_PAs_Te=meanvs.meanvalue(Te_PA.SIIISII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.SIIISII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.SIIISII=a
                        Ne_fake.SIIISII=b

                        Te.SIIISII[i,j]=a1
                        Ne.SIIISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.SIIISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.SIIISII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.SIIISII=a
                        Ne_fake.SIIISII=b
                                                
                        Te.SIIISII[i,j]=a1
                        Ne.SIIISII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.SIIISII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.SIIISII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: medium zones (SIII-SII): Te_SIII = {0:5.0f} K, Ne_SII = {1:7.1f} cm-1'.format(Te.SIIISII[i,j], Ne.SIIISII[i,j]))    
                
                if (index_Te_OII=="yes" and index_Ne_OII=="yes"):
                    a, b = diags.getCrossTemDen('[OII] 3727+/7325+', '[OII] 3726/3729', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OII] 3727+/7325+', '[OII] 3726/3729', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OIIOII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OIIOII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.OIIOII=a
                        Ne_fake.OIIOII=b
                        
                        Te.OIIOII[i,j]=a1
                        Ne.OIIOII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIIOII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIIOII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OIIOII=a
                        Ne_fake.OIIOII=b
                        
                        Te.OIIOII[i,j]=a1
                        Ne.OIIOII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIIOII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIIOII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: low zones (OII-OII): Te_OII = {0:5.0f} K, Ne_OII = {1:7.1f} cm-1'.format(Te.OIIOII[i,j], Ne.OIIOII[i,j]))
    
                if (index_Te_NII=="yes" and index_Ne_OII=="yes"):
                    a, b = diags.getCrossTemDen('[NII] 5755/6584+', '[OII] 3726/3729', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+', '[OII] 3726/3729', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.NIIOII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.NIIOII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.NIIOII=a
                        Ne_fake.NIIOII=b
                        
                        Te.NIIOII[i,j]=a1
                        Ne.NIIOII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.NIIOII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.NIIOII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.NIIOII=a
                        Ne_fake.NIIOII=b
                        
                        Te.NIIOII[i,j]=a1
                        Ne.NIIOII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.NIIOII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.NIIOII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: low zones (NII-OII): Te_NII = {0:5.0f} K, Ne_OII = {1:7.1f} cm-1'.format(Te.NIIOII[i,j], Ne.NIIOII[i,j]))
                            
                if (index_Te_OI=="yes" and index_Te_OI2=="no" and index_Ne_OII=="yes"):
                    a, b = diags.getCrossTemDen('[OI] 5577/6300+', '[OII] 3726/3729', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OI] 5577/6300+', '[OII] 3726/3729', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OIOII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OIOII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.OIOII=a
                        Ne_fake.OIOII=b
                        
                        Te.OIOII[i,j]=a1
                        Ne.OIOII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIOII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIOII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OIOII=a
                        Ne_fake.OIOII=b
                        
                        Te.OIOII[i,j]=a1
                        Ne.OIOII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIOII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIOII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: low zones (OI-OII): Te_OI = {0:5.0f} K, Ne_OII = {1:7.1f} cm-1'.format(Te.OIOII[i,j], Ne.OIOII[i,j]))
  
                
                if (index_Te_OI=="no" and index_Te_OI2=="yes" and index_Ne_OII=="yes"):
                    a, b = diags.getCrossTemDen('[OI] 5577/6300', '[OII] 3726/3729', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OI] 5577/6300', '[OII] 3726/3729', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OIOII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OIOII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.OIOII=a
                        Ne_fake.OIOII=b
                        
                        Te.OIOII[i,j]=a1
                        Ne.OIOII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIOII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIOII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OIOII=a
                        Ne_fake.OIOII=b
                        
                        Te.OIOII[i,j]=a1
                        Ne.OIOII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIOII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIOII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: low zones (OI-OII): Te_OI = {0:5.0f} K, Ne_OII = {1:7.1f} cm-1'.format(Te.OIOII[i,j], Ne.OIOII[i,j]))

                if (index_Te_OIII=="yes" and index_Ne_ClIII=="yes"):
                    a, b = diags.getCrossTemDen('[OIII] 4363/5007+', '[ClIII] 5538/5518', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+', '[ClIII] 5538/5518', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OIIIClIII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OIIIClIII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.OIIIClIII=a
                        Ne_fake.OIIIClIII=b
                        
                        Te.OIIIClIII[i,j]=a1
                        Ne.OIIIClIII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIIIClIII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIIIClIII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OIIIClIII=a
                        Ne_fake.OIIIClIII=b
                        
                        Te.OIIIClIII[i,j]=a1
                        Ne.OIIIClIII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIIIClIII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIIIClIII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: medium zones (OIII-ClIII): Te_OIII = {0:5.0f} K, Ne_ClIII = {1:7.1f} cm-1'.format(Te.OIIIClIII[i,j], Ne.OIIIClIII[i,j]))          
                  
                if (index_Te_SIII=="yes" and index_Ne_ClIII=="yes"):
                    a, b = diags.getCrossTemDen('[SIII] 6312/9069', '[ClIII] 5538/5518', obs=obs_fake)         
                    a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069', '[ClIII] 5538/5518', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.SIIIClIII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.SIIIClIII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.SIIIClIII=a
                        Ne_fake.SIIIClIII=b
                        
                        Te.SIIIClIII[i,j]=a1
                        Ne.SIIIClIII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.SIIIClIII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.SIIIClIII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.SIIIClIII=a
                        Ne_fake.SIIIClIII=b
                        
                        Te.SIIIClIII[i,j]=a1
                        Ne.SIIIClIII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.SIIIClIII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.SIIIClIII[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: medium zones (SIII-ClIII): Te_SIII = {0:5.0f} K, Ne_ClIII = {1:7.1f} cm-1'.format(Te.SIIIClIII[i,j], Ne.SIIIClIII[i,j]))
    
                if (index_Te_OIII=="yes" and index_Ne_ArVI=="yes"):
                    a, b = diags.getCrossTemDen('[OIII] 4363/5007+', '[ArIV] 4740/4711', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+', '[ArIV] 4740/4711', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.OIIIArVI,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.OIIIArVI,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.OIIIArVI=a
                        Ne_fake.OIIIArVI=b
                        
                        Te.OIIIArVI[i,j]=a1
                        Ne.OIIIArVI[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIIIArVI[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIIIArVI[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.OIIIArVI=a
                        Ne_fake.OIIIArVI=b
                        
                        Te.OIIIArVI[i,j]=a1
                        Ne.OIIIArVI[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.OIIIArVI[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.OIIIArVI[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: medium zones (OIII-ArVI): Te_OIII = {0:5.0f} K, Ne_ArVI = {1:7.1f} cm-1'.format(Te.OIIIArVI[i,j], Ne.OIIIArVI[i,j]))          
        
                if (index_Te_SIII=="yes" and index_Ne_ArVI=="yes"):
                    a, b = diags.getCrossTemDen('[SIII] 6312/9069', '[ArIV] 4740/4711', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069', '[ArIV] 4740/4711', obs=obs)
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.SIIIArVI,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.SIIIArVI,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.SIIIArVI=a
                        Ne_fake.SIIIArVI=b
                        
                        Te.SIIIArVI[i,j]=a1
                        Ne.SIIIArVI[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.SIIIArVI[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.SIIIArVI[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.SIIIArVI=a
                        Ne_fake.SIIIArVI=b
                        
                        Te.SIIIArVI[i,j]=a1
                        Ne.SIIIArVI[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.SIIIArVI[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.SIIIArVI[i,j]=np.std(np.nan_to_num(b))
#                    print('PyNeb: medium zones (SIII-ArVI): Te_SIII = {0:5.0f} K, Ne_ArVI = {1:7.1f} cm-1'.format(Te.SIIIArVI[i,j], Ne.SIIIArVI[i,j]))
                
                if (index_Te_NII=="yes" and index_Ne_ClIII=="yes"):
                    a, b = diags.getCrossTemDen('[NII] 5755/6584+', '[ClIII] 5538/5518', obs=obs_fake)
                    a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+', '[ClIII] 5538/5518', obs=obs)
                    
                    
                    mean_PAs_Te=meanvs.meanvalue(Te_PA.NIIClIII,0.0)
                    mean_PAs_Ne=meanvs.meanvalue(Ne_PA.NIIClIII,0.0)
                    if (np.isnan([a1]) == True): a1=0
                    if (np.isnan([b1]) == True): b1=0
                    
                    if a1>=0.33*mean_PAs_Te and a1<=3*mean_PAs_Te and b1>=0.1*mean_PAs_Ne and b1<=10*mean_PAs_Ne:
                        Te_fake.NIIClIII=a
                        Ne_fake.NIIClIII=b
                        
                        Te.NIIClIII[i,j]=a1
                        Ne.NIIClIII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.NIIClIII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.NIIClIII[i,j]=np.std(np.nan_to_num(b))
                    else:
                        Te_fake.NIIClIII=a
                        Ne_fake.NIIClIII=b
                        
                        Te.NIIClIII[i,j]=a1
                        Ne.NIIClIII[i,j]=b1
                        
                        a[a == inf] =  0
                        a=np.nan_to_num(a)
                        a = [i2 for i2 in a if i2 != 0]
                        a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

                        b[b == inf] =  0
                        b=np.nan_to_num(b)
                        b = [i2 for i2 in b if i2 != 0]
                        b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
                        
                        Te_error.NIIClIII[i,j]=np.std(np.nan_to_num(a))
                        Ne_error.NIIClIII[i,j]=np.std(np.nan_to_num(b))
#                    print >>file10, 'PyNeb: medium zones (NII-ClIII): Te_NII  = {0:5.0f},Ne_ClIII = {1:7.1f}  cm-1'.format(Te.NIIClIII[i],Ne.NIIClIII[i])
  
  
                ##############################################
                ### calculate, print and save ionic abundances
                ##############################################
                
                ion_abun,ion_abun_error,atom_abun_fake=ia2Ds.ionicabundances(obs,obs_fake,Te,Te_error,Te_fake,Ne,Ne_error,Ne_fake,i,j,line_names,lines_available,ion_abun,ion_abun_error,param_estimated,param_requered,Te_PA,Ne_PA,param_mod_name)
                
                ion_abun.HeIa=np.nan_to_num(ion_abun.HeIa)
                ion_abun.HeIb=np.nan_to_num(ion_abun.HeIb)
                ion_abun.HeIIa=np.nan_to_num(ion_abun.HeIIa)
                ion_abun.HeIIb=np.nan_to_num(ion_abun.HeIIb)
                ion_abun.NIIa=np.nan_to_num(ion_abun.NIIa)
                ion_abun.NIIb=np.nan_to_num(ion_abun.NIIb)
                ion_abun.NIIc=np.nan_to_num(ion_abun.NIIc)
                ion_abun.OIIIa=np.nan_to_num(ion_abun.OIIIa)
                ion_abun.OIIIb=np.nan_to_num(ion_abun.OIIIb)
                ion_abun.OIIIc=np.nan_to_num(ion_abun.OIIIc)
                ion_abun.OIIa=np.nan_to_num(ion_abun.OIIa)
                ion_abun.OIIb=np.nan_to_num(ion_abun.OIIb)
                ion_abun.OIIc=np.nan_to_num(ion_abun.OIIc)
                ion_abun.OIId=np.nan_to_num(ion_abun.OIId)
                ion_abun.SIIa=np.nan_to_num(ion_abun.SIIa)
                ion_abun.SIIb=np.nan_to_num(ion_abun.SIIb)
                ion_abun.SIIIa=np.nan_to_num(ion_abun.SIIIa)
                ion_abun.SIIIb=np.nan_to_num(ion_abun.SIIIb)
                ion_abun.NeIIIa=np.nan_to_num(ion_abun.NeIIIa)
                ion_abun.NeIIIb=np.nan_to_num(ion_abun.NeIIIb)
                ion_abun.ArIII=np.nan_to_num(ion_abun.ArIII)
                ion_abun.ArIVa=np.nan_to_num(ion_abun.ArIVa)
                ion_abun.ArIVb=np.nan_to_num(ion_abun.ArIVb)
                ion_abun.ClIIIa=np.nan_to_num(ion_abun.ClIIIa)
                ion_abun.ClIIIb=np.nan_to_num(ion_abun.ClIIIb)
                
                ion_abun_error.HeIa=np.nan_to_num(ion_abun_error.HeIa)
                ion_abun_error.HeIb=np.nan_to_num(ion_abun_error.HeIb)
                ion_abun_error.HeIIa=np.nan_to_num(ion_abun_error.HeIIa)
                ion_abun_error.HeIIb=np.nan_to_num(ion_abun_error.HeIIb)
                ion_abun_error.NIIa=np.nan_to_num(ion_abun_error.NIIa)
                ion_abun_error.NIIb=np.nan_to_num(ion_abun_error.NIIb)
                ion_abun_error.NIIc=np.nan_to_num(ion_abun_error.NIIc)
                ion_abun_error.OIIIa=np.nan_to_num(ion_abun_error.OIIIa)
                ion_abun_error.OIIIb=np.nan_to_num(ion_abun_error.OIIIb)
                ion_abun_error.OIIIc=np.nan_to_num(ion_abun_error.OIIIc)
                ion_abun_error.OIIa=np.nan_to_num(ion_abun_error.OIIa)
                ion_abun_error.OIIb=np.nan_to_num(ion_abun_error.OIIb)
                ion_abun_error.OIIc=np.nan_to_num(ion_abun_error.OIIc)
                ion_abun_error.OIId=np.nan_to_num(ion_abun_error.OIId)
                ion_abun_error.SIIa=np.nan_to_num(ion_abun_error.SIIa)
                ion_abun_error.SIIb=np.nan_to_num(ion_abun_error.SIIb)
                ion_abun_error.SIIIa=np.nan_to_num(ion_abun_error.SIIIa)
                ion_abun_error.SIIIb=np.nan_to_num(ion_abun_error.SIIIb)
                ion_abun_error.NeIIIa=np.nan_to_num(ion_abun_error.NeIIIa)
                ion_abun_error.NeIIIb=np.nan_to_num(ion_abun_error.NeIIIb)
                ion_abun_error.ArIII=np.nan_to_num(ion_abun_error.ArIII)
                ion_abun_error.ArIVa=np.nan_to_num(ion_abun_error.ArIVa)
                ion_abun_error.ArIVb=np.nan_to_num(ion_abun_error.ArIVb)
                ion_abun_error.ClIIIa=np.nan_to_num(ion_abun_error.ClIIIa)
                ion_abun_error.ClIIIb=np.nan_to_num(ion_abun_error.ClIIIb)
                
                ###########################
                ### calculate, print and save total elemental abundances and ICFs
                ##########################
                elem_abun_KB,elem_ICF_KB,elem_abun_DIMS,elem_ICF_DIMS,elem_abun_KB_error,elem_ICF_KB_error,elem_abun_DIMS_error,elem_ICF_DIMS_error,elem_abun_KB_ratio,elem_abun_KB_ratio_error,elem_abun_DIMS_ratio,elem_abun_DIMS_ratio_error=eaI2Ds.total_abundances_ICFs(ion_abun,ion_abun_error,atom_abun_fake,i,j,line_names,lines_available,elem_abun_KB,elem_ICF_KB,elem_abun_DIMS,elem_ICF_DIMS,elem_abun_KB_error,elem_ICF_KB_error,elem_abun_DIMS_error,elem_ICF_DIMS_error,elem_abun_KB_ratio,elem_abun_KB_ratio_error,elem_abun_DIMS_ratio,elem_abun_DIMS_ratio_error,param_mod_name)
    
                if elem_abun_KB.He[i,j]>1:
                    elem_abun_KB.He[i,j]=0
                if elem_abun_DIMS.He[i,j]>1:
                    elem_abun_DIMS.He[i,j]=0
                if elem_ICF_KB.He[i,j]>100:
                    elem_ICF_KB.He[i,j]=0
                if elem_ICF_DIMS.He[i,j]>100:
                    elem_ICF_DIMS.He[i,j]=0
                if elem_abun_KB.O[i,j]>1:
                    elem_abun_KB.O[i,j]=0
                if elem_abun_DIMS.O[i,j]>1:
                    elem_abun_DIMS.O[i,j]=0
                if elem_ICF_KB.O[i,j]>100:
                    elem_ICF_KB.O[i,j]=0
                if elem_ICF_DIMS.O[i,j]>100:
                    elem_ICF_DIMS.O[i,j]=0
                if elem_abun_KB.N[i,j]>1:
                    elem_abun_KB.N[i,j]=0
                if elem_abun_DIMS.N[i,j]>1:
                    elem_abun_DIMS.N[i,j]=0
                if elem_ICF_KB.N[i,j]>100:
                    elem_ICF_KB.N[i,j]=0
                if elem_ICF_DIMS.N[i,j]>100:
                    elem_ICF_DIMS.N[i,j]=0
                if elem_abun_KB.S[i,j]>1:
                    elem_abun_KB.S[i,j]=0
                if elem_abun_DIMS.S[i,j]>1:
                    elem_abun_DIMS.S[i,j]=0
                if elem_ICF_KB.S[i,j]>100:
                    elem_ICF_KB.S[i,j]=0
                if elem_ICF_DIMS.S[i,j]>100:
                    elem_ICF_DIMS.S[i,j]=0
                if elem_abun_KB.Ne[i,j]>1:
                    elem_abun_KB.Ne[i,j]=0
                if elem_abun_DIMS.Ne[i,j]>1:
                    elem_abun_DIMS.Ne[i,j]=0
                if elem_ICF_KB.Ne[i,j]>100:
                    elem_ICF_KB.Ne[i,j]=0
                if elem_ICF_DIMS.Ne[i,j]>100:
                    elem_ICF_DIMS.Ne[i,j]=0
                if elem_abun_KB.Ar[i,j]>1:
                    elem_abun_KB.Ar[i,j]=0
                if elem_abun_DIMS.Ar[i,j]>1:
                    elem_abun_DIMS.Ar[i,j]=0
                if elem_ICF_KB.Ar[i,j]>100:
                    elem_ICF_KB.Ar[i,j]=0
                if elem_ICF_DIMS.Ar[i,j]>100:
                    elem_ICF_DIMS.Ar[i,j]=0
                if elem_abun_KB.Cl[i,j]>1:
                    elem_abun_KB.Cl[i,j]=0
                if elem_abun_DIMS.Cl[i,j]>1:
                    elem_abun_DIMS.Cl[i,j]=0
                if elem_ICF_KB.Cl[i,j]>100:
                    elem_ICF_KB.Cl[i,j]=0
                if elem_ICF_DIMS.Cl[i,j]>100:
                    elem_ICF_DIMS.Cl[i,j]=0
                
                if elem_abun_KB_ratio.NO[i,j]>1:
                    elem_abun_KB_ratio.NO[i,j]=0
                if elem_abun_DIMS_ratio.NO[i,j]>1:
                    elem_abun_DIMS_ratio.NO[i,j]=0
                if elem_abun_KB_ratio.SO[i,j]>1:
                    elem_abun_KB_ratio.SO[i,j]=0
                if elem_abun_DIMS_ratio.SO[i,j]>1:
                    elem_abun_DIMS_ratio.SO[i,j]=0
                if elem_abun_KB_ratio.NeO[i,j]>1:
                    elem_abun_KB_ratio.NeO[i,j]=0
                if elem_abun_DIMS_ratio.NeO[i,j]>1:
                    elem_abun_DIMS_ratio.NeO[i,j]=0
                if elem_abun_KB_ratio.ArO[i,j]>1:
                    elem_abun_KB_ratio.ArO[i,j]=0
                if elem_abun_DIMS_ratio.ArO[i,j]>1:
                    elem_abun_DIMS_ratio.ArO[i,j]=0
                if elem_abun_KB_ratio.ClO[i,j]>1:
                    elem_abun_KB_ratio.ClO[i,j]=0
                if elem_abun_DIMS_ratio.ClO[i,j]>1:
                    elem_abun_DIMS_ratio.ClO[i,j]=0
                    
                os.remove("test.dat")
    pn.log_.close_file()
       
       
    ##################################################################################################
    ### Save fits images of C(Hb), C(Hg), C(Hd)
    ##################################################################################################
    sfs.save_image2(Chbeta_pyneb,hdr,'c_Hb','c_Hb_pyneb',param_estimated,param_requered)    
    sfs.save_image2(Chbeta_error_pyneb,hdr,'c_Hb_error','c_Hb_error_pyneb',param_estimated,param_requered)    
                

    ##################################################################################################
    ### Save fits images of temperature and density electron diagnostics
    ##################################################################################################
    sfs.save_image(Te.NIISII,hdr,'Te(NII6548_84)_Ne(SII6716_31)',param_estimated,param_requered)
    sfs.save_image(Ne.NIISII,hdr,'Ne(SII6716_31)_Te(NII6548_84)',param_estimated,param_requered)
 
    sfs.save_image2(Te_error.NIISII,hdr,'Te(NII6548_84)_Ne(SII6716_31)','Te(NII6548_84)_Ne(SII6716_31)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.NIISII,hdr,'Ne(SII6716_31)_Te(NII6548_84)','Ne(SII6716_31)_Te(NII6548_84)_err',param_estimated,param_requered)
 
    sfs.save_image(Te.OISII,hdr,'Te(OI6300_63)_Ne(SII6716_31)',param_estimated,param_requered)
    sfs.save_image(Ne.OISII,hdr,'Ne(SII6716_31)_Te(OI6300_63)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.OISII,hdr,'Te(OI6300_63)_Ne(SII6716_31)','Te(OI6300_63)_Ne(SII6716_31)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.OISII,hdr,'Ne(SII6716_31)_Te(OI6300_63)','Ne(SII6716_31)_Te(OI6300_63)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.OIISII,hdr,'Te(OII3727_29_7320_30)_Ne(SII6716_31)',param_estimated,param_requered)
    sfs.save_image(Ne.OIISII,hdr,'Ne(SII6716_31)_Te(OII3727_29_7320_30)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.OIISII,hdr,'Te(OII3727_29_7320_30)_Ne(SII6716_31)','Te(OII3727_29_7320_30)_Ne(SII6716_31)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.OIISII,hdr,'Ne(SII6716_31)_Te(OII3727_29_7320_30)','Ne(SII6716_31)_Te(OII3727_29_7320_30)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.OIIISII,hdr,'Te(OIII4959_5007)_Ne(SII6716_31)',param_estimated,param_requered)
    sfs.save_image(Ne.OIIISII,hdr,'Ne(SII6716_31)_Te(OIII4959_5007)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.OIIISII,hdr,'Te(OIII4959_5007)_Ne(SII6716_31)','Te(OIII4959_5007)_Ne(SII6716_31)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.OIIISII,hdr,'Ne(SII6716_31)_Te(OIII4959_5007)','Ne(SII6716_31)_Te(OIII4959_5007)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.SIIISII,hdr,'Te(SIII6312_9069)_Ne(SII6716_31)',param_estimated,param_requered)
    sfs.save_image(Ne.SIIISII,hdr,'Ne(SII6716_31)_Te(SIII6312_9069)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.SIIISII,hdr,'Te(SIII6312_9069)_Ne(SII6716_31)','Te(SIII6312_9069)_Ne(SII6716_31)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.SIIISII,hdr,'Ne(SII6716_31)_Te(SIII6312_9069)','Ne(SII6716_31)_Te(SIII6312_9069)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.OIIOII,hdr,'Te(OII3727_29_7320_30)_Ne(OII3727_29)',param_estimated,param_requered)
    sfs.save_image(Ne.OIIOII,hdr,'Ne(OII3727_29)_Te(OII3727_29_7320_30)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.OIIOII,hdr,'Te(OII3727_29_7320_30)_Ne(OII3727_29)','Te(OII3727_29_7320_30)_Ne(OII3727_29)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.OIIOII,hdr,'Ne(OII3727_29)_Te(OII3727_29_7320_30)','Ne(OII3727_29)_Te(OII3727_29_7320_30)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.NIIOII,hdr,'Te(NII6548_84)_Ne(OII3727_29)',param_estimated,param_requered)
    sfs.save_image(Ne.NIIOII,hdr,'Ne(OII3727_29)_Te(NII6548_84)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.NIIOII,hdr,'Te(NII6548_84)_Ne(OII3727_29)','Te(NII6548_84)_Ne(OII3727_29)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.NIIOII,hdr,'Ne(OII3727_29)_Te(NII6548_84)','Ne(OII3727_29)_Te(NII6548_84)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.OIOII,hdr,'Te(OI6300_63)_Ne(OII3727_29)',param_estimated,param_requered)
    sfs.save_image(Ne.OIOII,hdr,'Ne(OII3727_29)_Te(OI6300_63)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.OIOII,hdr,'Te(OI6300_63)_Ne(OII3727_29)','Te(OI6300_63)_Ne(OII3727_29)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.OIOII,hdr,'Ne(OII3727_29)_Te(OI6300_63)','Ne(OII3727_29)_Te(OI6300_63)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.OIIIClIII,hdr,'Te(OIII4959_5007)_Ne(ClIII5517_38)',param_estimated,param_requered)
    sfs.save_image(Ne.OIIIClIII,hdr,'Ne(ClIII5517_38)_Te(OIII4959_5007)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.OIIIClIII,hdr,'Te(OIII4959_5007)_Ne(ClIII5517_38)','Te(OIII4959_5007)_Ne(ClIII5517_38)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.OIIIClIII,hdr,'Ne(ClIII5517_38)_Te(OIII4959_5007)','Ne(ClIII5517_38)_Te(OIII4959_5007)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.SIIIClIII,hdr,'Te(SIII6312_9069)_Ne(ClIII5517_38)',param_estimated,param_requered)
    sfs.save_image(Ne.SIIIClIII,hdr,'Ne(ClIII5517_38)_Te(SIII6312_9069)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.SIIIClIII,hdr,'Te(SIII6312_9069)_Ne(ClIII5517_38)','Te(SIII6312_9069)_Ne(ClIII5517_38)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.SIIIClIII,hdr,'Ne(ClIII5517_38)_Te(SIII6312_9069)','Ne(ClIII5517_38)_Te(SIII6312_9069)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.OIIIArVI,hdr,'Te(OIII4959_5007)_Ne(ArVI4712_40)',param_estimated,param_requered)
    sfs.save_image(Ne.OIIIArVI,hdr,'Ne(ArVI4712_40)_Te(OIII4959_5007)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.OIIIArVI,hdr,'Te(OIII4959_5007)_Ne(ArVI4712_40)','Te(OIII4959_5007)_Ne(ArVI4712_40)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.OIIIArVI,hdr,'Ne(ArVI4712_40)_Te(OIII4959_5007)','Ne(ArVI4712_40)_Te(OIII4959_5007)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.SIIIArVI,hdr,'Te(SIII6312_9069)_Ne(ArVI4712_40)',param_estimated,param_requered)
    sfs.save_image(Ne.SIIIArVI,hdr,'Ne(ArVI4712_40)_Te(SIII6312_9069)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.SIIIArVI,hdr,'Te(SIII6312_9069)_Ne(ArVI4712_40)','Te(SIII6312_9069)_Ne(ArVI4712_40)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.SIIIArVI,hdr,'Ne(ArVI4712_40)_Te(SIII6312_9069)','Ne(ArVI4712_40)_Te(SIII6312_9069)_err',param_estimated,param_requered)
    
    sfs.save_image(Te.NIIClIII,hdr,'Te(NII6548_84)_Ne(ClIII5517_38)',param_estimated,param_requered)
    sfs.save_image(Ne.NIIClIII,hdr,'Ne(ClIII5517_38)_Te(NII6548_84)',param_estimated,param_requered)
    
    sfs.save_image2(Te_error.NIIClIII,hdr,'Te(NII6548_84)_Ne(ClIII5517_38)','Te(NII6548_84)_Ne(ClIII5517_38)_err',param_estimated,param_requered)
    sfs.save_image2(Ne_error.NIIClIII,hdr,'Ne(ClIII5517_38)_Te(NII6548_84)','Ne(ClIII5517_38)_Te(NII6548_84)_err',param_estimated,param_requered)
    
    ##################################################################################################
    ### Save fits images of ionic abundances and errors
    ##################################################################################################

    sfs.save_image_ionabun(ion_abun.HeIa,hdr,'HeI_5876s','He+a(5876)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.HeIb,hdr,'HeI_6678s','He+b(6678)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.HeIIa,hdr,'HeII_4686s','He++a(4686)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.HeIIb,hdr,'HeII_5412s','He++b(5412)',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun.OIa,hdr,'O1_5577s','Ooa(5577)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.OIb,hdr,'O1_6300s','Oob(6300)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.OIc,hdr,'O1_6363s','Ooc(6363)',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun.OIIa,hdr,'O2_3727s','O+a(3727)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.OIIb,hdr,'O2_3729s','O+b(3729)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.OIIc,hdr,'O2_7320s','O+c(7320)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.OIId,hdr,'O2_7330s','O+d(7330)',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun.OIIIa,hdr,'O3_4363s','O++a(4363)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.OIIIb,hdr,'O3_4959s','O++b(4959)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.OIIIc,hdr,'O3_5007s','O++c(5007)',param_estimated,param_requered,line_names,lines_available)

    sfs.save_image_ionabun(ion_abun.NI,hdr,'N1_5199s','No(5199)',param_estimated,param_requered,line_names,lines_available)    
    sfs.save_image_ionabun(ion_abun.NIIa,hdr,'N2_5755s','N+a(5755)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.NIIb,hdr,'N2_6548s','N+b(6548)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.NIIc,hdr,'N2_6583s','N+c(6584)',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun.SIIa,hdr,'S2_6716s','S+a(6716)',param_estimated,param_requered,line_names,lines_available)    
    sfs.save_image_ionabun(ion_abun.SIIb,hdr,'S2_6731s','S+b(6731)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.SIIIa,hdr,'S3_6312s','S++a(6312)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.SIIIb,hdr,'S3_9069s','S++b(9069)',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun.ClIIIa,hdr,'Cl3_5517s','Cl++a(5517)',param_estimated,param_requered,line_names,lines_available)    
    sfs.save_image_ionabun(ion_abun.ClIIIb,hdr,'Cl3_5538s','Cl++b(5538)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.ArIII,hdr,'Ar3_7136s','Ar++(7136)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.ArIVa,hdr,'Ar4_4712s','Ar+++a(4712)',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun.ArIVb,hdr,'Ar4_4740s','Ar+++b(4740)',param_estimated,param_requered,line_names,lines_available)

    sfs.save_image_ionabun(ion_abun.NeIIIa,hdr,'Ne3_3868s','Ne++a(3868)',param_estimated,param_requered,line_names,lines_available)    
    sfs.save_image_ionabun(ion_abun.NeIIIb,hdr,'Ne3_3967s','Ne++b(3968)',param_estimated,param_requered,line_names,lines_available)

    ### errors#

    sfs.save_image_ionabun(ion_abun_error.HeIa,hdr,'HeI_5876s','He+a(5876)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.HeIb,hdr,'HeI_6678s','He+b(6678)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.HeIIa,hdr,'HeII_4686s','He++a(4686)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.HeIIb,hdr,'HeII_5412s','He++b(5412)_err',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun_error.OIa,hdr,'O1_5577s','Ooa(5577)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.OIb,hdr,'O1_6300s','Oob(6300)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.OIc,hdr,'O1_6363s','Ooc(6363)_err',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun_error.OIIa,hdr,'O2_3727s','O+a(3727)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.OIIb,hdr,'O2_3729s','O+b(3729)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.OIIc,hdr,'O2_7320s','O+c(7320)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.OIId,hdr,'O2_7330s','O+d(7330)_err',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun_error.OIIIa,hdr,'O3_4363s','O++a(4363)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.OIIIb,hdr,'O3_4959s','O++b(4959)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.OIIIc,hdr,'O3_5007s','O++c(5007)_err',param_estimated,param_requered,line_names,lines_available)

    sfs.save_image_ionabun(ion_abun_error.NI,hdr,'N1_5199s','No(5199)_err',param_estimated,param_requered,line_names,lines_available)    
    sfs.save_image_ionabun(ion_abun_error.NIIa,hdr,'N2_5755s','N+a(5755)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.NIIb,hdr,'N2_6548s','N+b(6548)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.NIIc,hdr,'N2_6583s','N+c(6584)_err',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun_error.SIIa,hdr,'S2_6716s','S+a(6716)_err',param_estimated,param_requered,line_names,lines_available)    
    sfs.save_image_ionabun(ion_abun_error.SIIb,hdr,'S2_6731s','S+b(6731)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.SIIIa,hdr,'S3_6312s','S++a(6312)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.SIIIb,hdr,'S3_9069s','S++b(9069)_err',param_estimated,param_requered,line_names,lines_available)
    
    sfs.save_image_ionabun(ion_abun_error.ClIIIa,hdr,'Cl3_5517s','Cl++a(5517)_err',param_estimated,param_requered,line_names,lines_available)    
    sfs.save_image_ionabun(ion_abun_error.ClIIIb,hdr,'Cl3_5538s','Cl++b(5538)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.ArIII,hdr,'Ar3_7136s','Ar++(7136)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.ArIVa,hdr,'Ar4_4712s','Ar+++a(4712)_err',param_estimated,param_requered,line_names,lines_available)
    sfs.save_image_ionabun(ion_abun_error.ArIVb,hdr,'Ar4_4740s','Ar+++b(4740)_err',param_estimated,param_requered,line_names,lines_available)

    sfs.save_image_ionabun(ion_abun_error.NeIIIa,hdr,'Ne3_3868s','Ne++a(3868)_err',param_estimated,param_requered,line_names,lines_available)    
    sfs.save_image_ionabun(ion_abun_error.NeIIIb,hdr,'Ne3_3967s','Ne++b(3968)_err',param_estimated,param_requered,line_names,lines_available)
 
    ##################################################################################################
    ### Save fits images of total elemental abundances and ICFs
    ##################################################################################################
    
    sfs.save_image_ionabun(elem_abun_KB.He,hdr,'HeI_5876s','He_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
    sfs.save_image_ionabun(elem_abun_DIMS.He,hdr,'HeI_5876s','He_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
    sfs.save_image_ionabun(elem_abun_KB_error.He,hdr,'HeI_5876s','He_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
    sfs.save_image_ionabun(elem_abun_DIMS_error.He,hdr,'HeI_5876s','He_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
    
    index_line1=line_names.index("O3_5007s")
    index_line2=line_names.index("O3_4959s")
    if lines_available[index_line1]=="yes":
        sfs.save_image_ionabun(elem_abun_KB.O,hdr,'O3_5007s','O_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.O,hdr,'O3_5007s','O_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.O,hdr,'O3_5007s','O_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.O,hdr,'O3_5007s','O_ICF_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_error.O,hdr,'O3_5007s','O_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.O,hdr,'O3_5007s','O_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.O,hdr,'O3_5007s','O_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.O,hdr,'O3_5007s','O_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)         
        
        
    elif lines_available[index_line2]=="yes" and lines_available[index_line1]=="no":
        sfs.save_image_ionabun(elem_abun_KB.O,hdr,'O3_4959s','O_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.O,hdr,'O3_4959s','O_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.O,hdr,'O3_4959s','O_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.O,hdr,'O3_4959s','O_ICF_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_error.O,hdr,'O3_4959s','O_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.O,hdr,'O3_4959s','O_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.O,hdr,'O3_4959s','O_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.O,hdr,'O3_4959s','O_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
    
        
    index_line1=line_names.index("N2_6583s")
    index_line2=line_names.index("N2_6548s")
    if lines_available[index_line1]=="yes":
        sfs.save_image_ionabun(elem_abun_KB.N,hdr,'N2_6583s','N_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.N,hdr,'N2_6583s','N_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.N,hdr,'N2_6583s','N_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.N,hdr,'N2_6583s','N_ICF_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_error.N,hdr,'N2_6583s','N_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.N,hdr,'N2_6583s','N_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.N,hdr,'N2_6583s','N_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.N,hdr,'N2_6583s','N_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.NO,hdr,'N2_6583s','NO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.NO,hdr,'N2_6583s','NO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.NO,hdr,'N2_6583s','NO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.NO,hdr,'N2_6583s','NO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 

    elif lines_available[index_line2]=="yes" and lines_available[index_line1]=="no":
        sfs.save_image_ionabun(elem_abun_KB.N,hdr,'N2_6548s','N_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.N,hdr,'N2_6548s','N_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.N,hdr,'N2_6548s','N_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.N,hdr,'N2_6548s','N_ICF_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_error.N,hdr,'N2_6548s','N_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.N,hdr,'N2_6548s','N_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.N,hdr,'N2_6548s','N_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.N,hdr,'N2_6548s','N_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
    
        sfs.save_image_ionabun(elem_abun_KB_ratio.NO,hdr,'N2_6548s','NO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.NO,hdr,'N2_6548s','NO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.NO,hdr,'N2_6548s','NO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.NO,hdr,'N2_6548s','NO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
    
    index_line1=line_names.index("S2_6716s")
    index_line2=line_names.index("S2_6731s")
    if lines_available[index_line1]=="yes":
        sfs.save_image_ionabun(elem_abun_KB.S,hdr,'S2_6716s','S_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.S,hdr,'S2_6716s','S_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.S,hdr,'S2_6716s','S_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.S,hdr,'S2_6716s','S_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.S,hdr,'S2_6716s','S_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.S,hdr,'S2_6716s','S_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.S,hdr,'S2_6716s','S_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.S,hdr,'S2_6716s','S_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.SO,hdr,'S2_6716s','SO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.SO,hdr,'S2_6716s','SO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.SO,hdr,'S2_6716s','SO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.SO,hdr,'S2_6716s','SO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        
    elif lines_available[index_line2]=="yes":
        sfs.save_image_ionabun(elem_abun_KB.S,hdr,'S2_6731s','S_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.S,hdr,'S2_6731s','S_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.S,hdr,'S2_6731s','S_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.S,hdr,'S2_6731s','S_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.S,hdr,'S2_6731s','S_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.S,hdr,'S2_6731s','S_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.S,hdr,'S2_6731s','S_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.S,hdr,'S2_6731s','S_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.SO,hdr,'S2_6731s','SO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.SO,hdr,'S2_6731s','SO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.SO,hdr,'S2_6731s','SO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.SO,hdr,'S2_6731s','SO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        
    index_line1=line_names.index("Ne3_3868s")
    index_line2=line_names.index("Ne3_3967s")
    if lines_available[index_line1]=="yes":    
        sfs.save_image_ionabun(elem_abun_KB.Ne,hdr,'Ne3_3868s','Ne_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.Ne,hdr,'Ne3_3868s','Ne_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.Ne,hdr,'Ne3_3868s','Ne_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.Ne,hdr,'Ne3_3868s','Ne_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.Ne,hdr,'Ne3_3868s','Ne_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.Ne,hdr,'Ne3_3868s','Ne_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.Ne,hdr,'Ne3_3868s','Ne_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.Ne,hdr,'Ne3_3868s','Ne_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.NeO,hdr,'Ne3_3868s','NeO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.NeO,hdr,'Ne3_3868s','NeO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.NeO,hdr,'Ne3_3868s','NeO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.NeO,hdr,'Ne3_3868s','NeO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        
        
    elif lines_available[index_line2]=="yes" and lines_available[index_line1]=="no":    
        sfs.save_image_ionabun(elem_abun_KB.Ne,hdr,'Ne3_3967s','Ne_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.Ne,hdr,'Ne3_3967s','Ne_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.Ne,hdr,'Ne3_3967s','Ne_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.Ne,hdr,'Ne3_3967s','Ne_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.Ne,hdr,'Ne3_3967s','Ne_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.Ne,hdr,'Ne3_3967s','Ne_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.Ne,hdr,'Ne3_3967s','Ne_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.Ne,hdr,'Ne3_3967s','Ne_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
       
        sfs.save_image_ionabun(elem_abun_KB_ratio.NeO,hdr,'Ne3_3967s','NeO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.NeO,hdr,'Ne3_3967s','NeO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.NeO,hdr,'Ne3_3967s','NeO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.NeO,hdr,'Ne3_3967s','NeO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        
        
    index_line1=line_names.index("Ar3_7136s")
    index_line2=line_names.index("Ar4_4712s")
    index_line3=line_names.index("Ar4_4740s")
    if lines_available[index_line1]=="yes" and lines_available[index_line2]=="no" and lines_available[index_line3]=="no":     
        sfs.save_image_ionabun(elem_abun_KB.Ar,hdr,'Ar3_7136s','Ar_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.Ar,hdr,'Ar3_7136s','Ar_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.Ar,hdr,'Ar3_7136s','Ar_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.Ar,hdr,'Ar3_7136s','Ar_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.Ar,hdr,'Ar3_7136s','Ar_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.Ar,hdr,'Ar3_7136s','Ar_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.Ar,hdr,'Ar3_7136s','Ar_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.Ar,hdr,'Ar3_7136s','Ar_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.ArO,hdr,'Ar3_7136s','ArO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.ArO,hdr,'Ar3_7136s','ArO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.ArO,hdr,'Ar3_7136s','ArO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.ArO,hdr,'Ar3_7136s','ArO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        
    elif lines_available[index_line2]=="yes" and lines_available[index_line1]=="no":     
        sfs.save_image_ionabun(elem_abun_KB.Ar,hdr,'Ar4_4712s','Ar_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.Ar,hdr,'Ar4_4712s','Ar_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.Ar,hdr,'Ar4_4712s','Ar_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.Ar,hdr,'Ar4_4712s','Ar_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.Ar,hdr,'Ar4_4712s','Ar_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.Ar,hdr,'Ar4_4712s','Ar_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.Ar,hdr,'Ar4_4712s','Ar_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.Ar,hdr,'Ar4_4712s','Ar_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.ArO,hdr,'Ar4_4712s','ArO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.ArO,hdr,'Ar4_4712s','ArO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.ArO,hdr,'Ar4_4712s','ArO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.ArO,hdr,'Ar4_4712s','ArO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        
    elif lines_available[index_line3]=="yes" and lines_available[index_line1]=="no":     
        sfs.save_image_ionabun(elem_abun_KB.Ar,hdr,'Ar4_4740s','Ar_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.Ar,hdr,'Ar4_4740s','Ar_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.Ar,hdr,'Ar4_4740s','Ar_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.Ar,hdr,'Ar4_4740s','Ar_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.Ar,hdr,'Ar4_4740s','Ar_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.Ar,hdr,'Ar4_4740s','Ar_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.Ar,hdr,'Ar4_4740s','Ar_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.Ar,hdr,'Ar4_4740s','Ar_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.ArO,hdr,'Ar4_4740s','ArO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.ArO,hdr,'Ar4_4740s','ArO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.ArO,hdr,'Ar4_4740s','ArO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.ArO,hdr,'Ar4_4740s','ArO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
    
    
    
    index_line1=line_names.index("Cl3_5517s")
    index_line2=line_names.index("Cl3_5538s")
    if lines_available[index_line1]=="yes":  
        sfs.save_image_ionabun(elem_abun_KB.Cl,hdr,'Cl3_5517s','Cl_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.Cl,hdr,'Cl3_5517s','Cl_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.Cl,hdr,'Cl3_5517s','Cl_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.Cl,hdr,'Cl3_5517s','Cl_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.Cl,hdr,'Cl3_5517s','Cl_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.Cl,hdr,'Cl3_5517s','Cl_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.Cl,hdr,'Cl3_5517s','Cl_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.Cl,hdr,'Cl3_5517s','Cl_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.ClO,hdr,'Cl3_5517s','ClO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.ClO,hdr,'Cl3_5517s','ClO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.ClO,hdr,'Cl3_5517s','ClO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.ClO,hdr,'Cl3_5517s','ClO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
    
    elif lines_available[index_line2]=="yes":  
        sfs.save_image_ionabun(elem_abun_KB.Cl,hdr,'Cl3_5538s','Cl_abundance_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS.Cl,hdr,'Cl3_5538s','Cl_abundance_DIMS',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB.Cl,hdr,'Cl3_5538s','Cl_ICF_KB',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS.Cl,hdr,'Cl3_5538s','Cl_ICF_DIMS',param_estimated,param_requered,line_names,lines_available)
        sfs.save_image_ionabun(elem_abun_KB_error.Cl,hdr,'Cl3_5538s','Cl_abundance_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_error.Cl,hdr,'Cl3_5538s','Cl_abundance_DIMS_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_KB_error.Cl,hdr,'Cl3_5538s','Cl_ICF_KB_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_ICF_DIMS_error.Cl,hdr,'Cl3_5538s','Cl_ICF_DIMS_error',param_estimated,param_requered,line_names,lines_available)
        
        sfs.save_image_ionabun(elem_abun_KB_ratio.ClO,hdr,'Cl3_5538s','ClO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio.ClO,hdr,'Cl3_5538s','ClO_abundance_ratio',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_KB_ratio_error.ClO,hdr,'Cl3_5538s','ClO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
        sfs.save_image_ionabun(elem_abun_DIMS_ratio_error.ClO,hdr,'Cl3_5538s','ClO_abundance_ratio_error',param_estimated,param_requered,line_names,lines_available) 
    

    
    
    hs.histogram_plots(Chbeta_pyneb,Chbeta_error_pyneb,'c_Hb',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.NIISII,Te_error.NIISII,'Te(NII6548_84)_Ne(SII6716_31)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.NIISII,Ne_error.NIISII,'Ne(SII6716_31)_Te(NII6548_84)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.OISII,Te_error.OISII,'Te(OI6300_63)_Ne(SII6716_31)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.OISII,Ne_error.OISII,'Ne(SII6716_31)_Te(OI6300_63)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.OIISII,Te_error.OIISII,'Te(OII3727_29_7320_30)_Ne(SII6716_31)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.OIISII,Ne_error.OIISII,'Ne(SII6716_31)_Te(OII3727_29_7320_30)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.OIIISII,Te_error.OIIISII,'Te(OIII4959_5007)_Ne(SII6716_31)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.OIIISII,Ne_error.OIIISII,'Ne(SII6716_31)_Te(OIII4959_5007)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.SIIISII,Te_error.SIIISII,'Te(SIII6312_9069)_Ne(SII6716_31)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.SIIISII,Ne_error.SIIISII,'Ne(SII6716_31)_Te(SIII6312_9069)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.OIIOII,Te_error.OIIOII,'Te(OII3727_29_7320_30)_Ne(OII3727_29)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.OIIOII,Ne_error.OIIOII,'Ne(OII3727_29)_Te(OII3727_29_7320_30)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.NIIOII,Te_error.NIIOII,'Te(NII6548_84)_Ne(OII3727_29)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.NIIOII,Ne_error.NIIOII,'Ne(OII3727_29)_Te(NII6548_84)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.OIOII,Te_error.OIOII,'Te(OI6300_63)_Ne(OII3727_29)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.OIOII,Ne_error.OIOII,'Ne(OII3727_29)_Te(OI6300_63)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.OIIIClIII,Te_error.OIIIClIII,'Te(OIII4959_5007)_Ne(ClIII5517_38)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.OIIIClIII,Ne_error.OIIIClIII,'Ne(ClIII5517_38)_Te(OIII4959_5007)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.SIIIClIII,Te_error.SIIIClIII,'Te(SIII6312_9069)_Ne(ClIII5517_38)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.SIIIClIII,Ne_error.SIIIClIII,'Ne(ClIII5517_38)_Te(SIII6312_9069)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.OIIIArVI,Te_error.OIIIArVI,'Te(OIII4959_5007)_Ne(ArVI4712_40)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.OIIIArVI,Ne_error.OIIIArVI,'Ne(ArVI4712_40)_Te(OIII4959_5007)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(Te.SIIIArVI,Te_error.SIIIArVI,'Te(SIII6312_9069)_Ne(ArVI4712_40)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.SIIIArVI,Ne_error.SIIIArVI,'Ne(ArVI4712_40)_Te(SIII6312_9069)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)

    hs.histogram_plots(Te.NIIClIII,Te_error.NIIClIII,'Te(NII6548_84)_Ne(ClIII5517_38)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(Ne.NIIClIII,Ne_error.NIIClIII,'Ne(ClIII5517_38)_Te(NII6548_84)',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)

    
    hs.histogram_plots(elem_abun_KB.He,elem_abun_KB_error.He,'He_adund_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS.He,elem_abun_DIMS_error.He,'He_abund_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_KB.He,elem_ICF_KB_error.He,'He_ICF_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered) 
    hs.histogram_plots(elem_ICF_DIMS.He,elem_ICF_DIMS_error.He,'He_ICF_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  

    hs.histogram_plots(elem_abun_KB.O,elem_abun_KB_error.O,'O_adund_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS.O,elem_abun_DIMS_error.O,'O_abund_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_KB.O,elem_ICF_KB_error.O,'O_ICF_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered) 
    hs.histogram_plots(elem_ICF_DIMS.O,elem_ICF_DIMS_error.O,'O_ICF_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)
    
    hs.histogram_plots(elem_abun_KB.N,elem_abun_KB_error.N,'N_adund_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS.N,elem_abun_DIMS_error.N,'N_abund_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_KB.N,elem_ICF_KB_error.N,'N_ICF_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered) 
    hs.histogram_plots(elem_ICF_DIMS.N,elem_ICF_DIMS_error.N,'N_ICF_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)
    
    hs.histogram_plots(elem_abun_KB_ratio.NO,elem_abun_KB_ratio_error.NO,'NO_ratio_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS_ratio.NO,elem_abun_DIMS_ratio_error.NO,'NO_ratio_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    
    hs.histogram_plots(elem_abun_KB.S,elem_abun_KB_error.S,'S_adund_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS.S,elem_abun_DIMS_error.S,'S_abund_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_KB.S,elem_ICF_KB_error.S,'S_ICF_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_DIMS.S,elem_ICF_DIMS_error.S,'S_ICF_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)
    
    hs.histogram_plots(elem_abun_KB_ratio.SO,elem_abun_KB_ratio_error.SO,'SO_ratio_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS_ratio.SO,elem_abun_DIMS_ratio_error.SO,'SO_ratio_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    
    hs.histogram_plots(elem_abun_KB.Ne,elem_abun_KB_error.Ne,'Ne_adund_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS.Ne,elem_abun_DIMS_error.Ne,'Ne_abund_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_KB.Ne,elem_ICF_KB_error.Ne,'Ne_ICF_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_DIMS.Ne,elem_ICF_DIMS_error.Ne,'Ne_ICF_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)
    
    hs.histogram_plots(elem_abun_KB_ratio.NeO,elem_abun_KB_ratio_error.NeO,'NeO_ratio_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS_ratio.NeO,elem_abun_DIMS_ratio_error.NeO,'NeO_ratio_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    
    hs.histogram_plots(elem_abun_KB.Ar,elem_abun_KB_error.Ar,'Ar_adund_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS.Ar,elem_abun_DIMS_error.Ar,'Ar_abund_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_KB.Ar,elem_ICF_KB_error.Ar,'Ar_ICF_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_DIMS.Ar,elem_ICF_DIMS_error.Ar,'Ar_ICF_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)
    
    hs.histogram_plots(elem_abun_KB_ratio.ArO,elem_abun_KB_ratio_error.ArO,'ArO_ratio_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS_ratio.ArO,elem_abun_DIMS_ratio_error.ArO,'ArO_ratio_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    
    hs.histogram_plots(elem_abun_KB.Cl,elem_abun_KB_error.Cl,'Cl_adund_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS.Cl,elem_abun_DIMS_error.Cl,'Cl_abund_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_KB.Cl,elem_ICF_KB_error.Cl,'Cl_ICF_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_ICF_DIMS.Cl,elem_ICF_DIMS_error.Cl,'Cl_ICF_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)
    
    hs.histogram_plots(elem_abun_KB_ratio.ClO,elem_abun_KB_ratio_error.ClO,'ClO_ratio_KB',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    hs.histogram_plots(elem_abun_DIMS_ratio.ClO,elem_abun_DIMS_ratio_error.ClO,'ClO_ratio_DIMS',line_names,lines_available,param_mod_name,param_model_values,param_estimated,param_requered)  
    
    file11 = open('statisticnumbers_1.txt','w')
    print >>file11, '{:32}'.format('Name Parameter'), '{:10}'.format('Nspaxels'),'{:10}'.format('5% value'),'{:10}'.format('Q1 value'), '{:12}'.format('median value'), '{:10}'.format('Q3 value'), '{:10}'.format('95% value'), '{:10}'.format('mean value'), '{:10}'.format('sigma')
    
    index=param_estimated.index('c_Hb')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Chbeta_pyneb,sizex,sizey)
        print >>file11, '{:32}'.format('c(Hbeta)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(NII6548_84)_Ne(SII6716_31)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.NIISII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(NII6548_84)_Ne(SII6716_31)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(SII6716_31)_Te(NII6548_84)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.NIISII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(SII6716_31)_Te(NII6548_84)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(OI6300_63)_Ne(SII6716_31)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.OISII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(OI6300_63)_Ne(SII6716_31)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(SII6716_31)_Te(OI6300_63)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.OISII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(SII6716_31)_Te(OI6300_63)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(OII3727_29_7320_30)_Ne(SII6716_31)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.OIISII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(OII3727_29_7320_30)_Ne(SII6716_31)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(SII6716_31)_Te(OII3727_29_7320_30)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.OIISII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(SII6716_31)_Te(OII3727_29_7320_30)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
        
    index=param_estimated.index('Te(OIII4959_5007)_Ne(SII6716_31)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.OIIISII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(OIII4959_5007)_Ne(SII6716_31)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(SII6716_31)_Te(OIII4959_5007)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.OIIISII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(SII6716_31)_Te(OIII4959_5007)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(SIII6312_9069)_Ne(SII6716_31)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.SIIISII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(SIII6312_9069)_Ne(SII6716_31)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(SII6716_31)_Te(SIII6312_9069)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.SIIISII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(SII6716_31)_Te(SIII6312_9069)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(OII3727_29_7320_30)_Ne(OII3727_29)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.OIIOII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(OII3727_29_7320_30)_Ne(OII3727_29)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(OII3727_29)_Te(OII3727_29_7320_30)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.OIIOII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(OII3727_29)_Te(OII3727_29_7320_30)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(NII6548_84)_Ne(OII3727_29)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.NIIOII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(NII6548_84)_Ne(OII3727_29)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(OII3727_29)_Te(NII6548_84)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.NIIOII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(OII3727_29)_Te(NII6548_84)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(OI6300_63)_Ne(OII3727_29)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.OIOII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(OI6300_63)_Ne(OII3727_29)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(OII3727_29)_Te(OI6300_63)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.OIOII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(OII3727_29)_Te(OI6300_63)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(OIII4959_5007)_Ne(ClIII5517_38)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.OIIIClIII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(OIII4959_5007)_Ne(ClIII5517_38)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(ClIII5517_38)_Te(OIII4959_5007)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.OIIIClIII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(ClIII5517_38)_Te(OIII4959_5007)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(SIII6312_9069)_Ne(ClIII5517_38)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.SIIIClIII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(SIII6312_9069)_Ne(ClIII5517_38)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(ClIII5517_38)_Te(SIII6312_9069)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.SIIIClIII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(ClIII5517_38)_Te(SIII6312_9069)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(OIII4959_5007)_Ne(ArVI4712_40)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.OIIIArVI,sizex,sizey)
        print >>file11, '{:32}'.format('Te(OIII4959_5007)_Ne(ArVI4712_40)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(ArVI4712_40)_Te(OIII4959_5007)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.OIIIArVI,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(ArVI4712_40)_Te(OIII4959_5007)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(SIII6312_9069)_Ne(ArVI4712_40)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.SIIIArVI,sizex,sizey)
        print >>file11, '{:32}'.format('Te(SIII6312_9069)_Ne(ArVI4712_40)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(ArVI4712_40)_Te(SIII6312_9069)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.SIIIArVI,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(ArVI4712_40)_Te(SIII6312_9069)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Te(NII6548_84)_Ne(ClIII5517_38)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Te.NIIClIII,sizex,sizey)
        print >>file11, '{:32}'.format('Te(NII6548_84)_Ne(ClIII5517_38)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index=param_estimated.index('Ne(ClIII5517_38)_Te(NII6548_84)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(Ne.NIIClIII,sizex,sizey)
        print >>file11, '{:32}'.format('Ne(ClIII5517_38)_Te(NII6548_84)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)

    print >>file11, "###################################################################################"
    

    index_line2=line_names.index("HeI_5876s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.HeIa,sizex,sizey)
        print >>file11, '{:32}'.format('HeI_5876s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("HeI_6678s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.HeIb,sizex,sizey)
        print >>file11, '{:32}'.format('HeI_6678s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("HeII_4686s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.HeIIa,sizex,sizey)
        print >>file11, '{:32}'.format('HeII_4686s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("HeII_5412s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.HeIIb,sizex,sizey)
        print >>file11, '{:32}'.format('HeII_5412s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
  
    index_line2=line_names.index("O1_5577s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIa,sizex,sizey)
        print >>file11, '{:32}'.format('O1_5577s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("O1_6300s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIb,sizex,sizey)
        print >>file11, '{:32}'.format('O1_6300s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("O1_6363s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIc,sizex,sizey)
        print >>file11, '{:32}'.format('O1_6363s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
  
    index_line2=line_names.index("O2_3727s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIIa,sizex,sizey)
        print >>file11, '{:32}'.format('O2_3727s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
  
    index_line2=line_names.index("O2_3729s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIIb,sizex,sizey)
        print >>file11, '{:32}'.format('O2_3729s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("O2_7320s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIIc,sizex,sizey)
        print >>file11, '{:32}'.format('O2_7320s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("O2_7330s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIId,sizex,sizey)
        print >>file11, '{:32}'.format('O2_7330s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("O3_4363s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIIIa,sizex,sizey)
        print >>file11, '{:32}'.format('O3_4363s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("O3_4959s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIIIb,sizex,sizey)
        print >>file11, '{:32}'.format('O3_4959s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("O3_5007s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.OIIIc,sizex,sizey)
        print >>file11, '{:32}'.format('O3_5007s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
 
    index_line2=line_names.index("N1_5199s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.NI,sizex,sizey)
        print >>file11, '{:32}'.format('N1_5199s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)

    index_line2=line_names.index("N2_5755s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.NIIa,sizex,sizey)
        print >>file11, '{:32}'.format('N2_5755s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("N2_6548s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.NIIb,sizex,sizey)
        print >>file11, '{:32}'.format('N2_6548s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("N2_6583s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.NIIc,sizex,sizey)
        print >>file11, '{:32}'.format('N2_6583s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
 
    index_line2=line_names.index("S2_6716s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.SIIa,sizex,sizey)
        print >>file11, '{:32}'.format('S2_6716s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)

    index_line2=line_names.index("S2_6731s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.SIIb,sizex,sizey)
        print >>file11, '{:32}'.format('S2_6731s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("S3_6312s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.SIIIa,sizex,sizey)
        print >>file11, '{:32}'.format('S3_6312s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("S3_9069s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.SIIIb,sizex,sizey)
        print >>file11, '{:32}'.format('S3_9069s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
 
    index_line2=line_names.index("Cl3_5517s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.ClIIIa,sizex,sizey)
        print >>file11, '{:32}'.format('Cl3_5517s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("Cl3_5538s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.ClIIIb,sizex,sizey)
        print >>file11, '{:32}'.format('Cl3_5538s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("Ar3_7136s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.ArIII,sizex,sizey)
        print >>file11, '{:32}'.format('Ar3_7136s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
 
    index_line2=line_names.index("Ar4_4712s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.ArIVa,sizex,sizey)
        print >>file11, '{:32}'.format('Ar4_4712s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)

    index_line2=line_names.index("Ar4_4740s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.ArIVb,sizex,sizey)
        print >>file11, '{:32}'.format('Ar4_4740s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("Ne3_3868s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.NeIIIa,sizex,sizey)
        print >>file11, '{:32}'.format('Ne3_3868s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
    
    index_line2=line_names.index("Ne3_3967s")
    if lines_available[index_line2]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(ion_abun.NeIIIb,sizex,sizey)
        print >>file11, '{:32}'.format('Ne3_3967s'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue)
 
    
    file11.close() 
    
    return Te,Ne       
