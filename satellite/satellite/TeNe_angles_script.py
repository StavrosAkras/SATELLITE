# TeNe_angles_script.py:
# Computes c(Hb), Te and Ne parameters for all the pseudo-slits in the rotation analysis module using the PyNeb package (Luridiana et al. 2015) 
# as well as the ionic/elemental abundances and ICFs.
# The routine returns plots for all these parameters as a function of the position angles.
# The results are also provided in an ascii file (PyNeb_output_total_abund_ICFs_per_angles.txt 
# and PyNeb_output_ionic_abund_per_angles.txt and PyNeb_output_Te_and_Ne_per_angles.txt)
# (C) Stavros Akras


import numpy as np
import matplotlib.pyplot as plt
import pyneb as pn
from numpy import inf
import os

from satellite import ionicabundances_angles_script as iaas
from satellite import element_abundances_ICFs_angles_script as eaIas
from satellite import norm_flux_error_script as nfes

global file10

def TeNe(flux_angles,flux_angles_error,ang,line_names,lines_available,param_estimated,param_requered,param_mod_name,param_model_values,par_plotname,par_plotymin,par_plotymax):

    class Te:
        NIISII=[]
        OISII=[]
        OIISII=[]
        OIIISII=[]
        SIIISII=[]
        OIIOII=[]
        NIIOII=[]
        OIOII=[]
        OIIIClIII=[]
        SIIIClIII=[]
        OIIIArVI=[]
        SIIIArVI=[]
        NIIClIII=[]

    class Ne:
        NIISII=[]
        OISII=[]
        OIISII=[]
        OIIISII=[]
        SIIISII=[]
        OIIOII=[]
        NIIOII=[]
        OIOII=[]
        OIIIClIII=[]
        SIIIClIII=[]
        OIIIArVI=[]
        SIIIArVI=[]
        NIIClIII=[]
        
    class Te_error:
        NIISII=[]
        OISII=[]
        OIISII=[]
        OIIISII=[]
        SIIISII=[]
        OIIOII=[]
        NIIOII=[]
        OIOII=[]
        OIIIClIII=[]
        SIIIClIII=[]
        OIIIArVI=[]
        SIIIArVI=[]
        NIIClIII=[]
        
    class Ne_error:
        NIISII=[]
        OISII=[]
        OIISII=[]
        OIIISII=[]
        SIIISII=[]
        OIIOII=[]
        NIIOII=[]
        OIOII=[]
        OIIIClIII=[]
        SIIIClIII=[]
        OIIIArVI=[]
        SIIIArVI=[]
        NIIClIII=[]
        
    class Te_fake:
        NIISII=[]
        OISII=[]
        OIISII=[]
        OIIISII=[]
        SIIISII=[]
        OIIOII=[]
        NIIOII=[]
        OIOII=[]
        OIIIClIII=[]
        SIIIClIII=[]
        OIIIArVI=[]
        SIIIArVI=[]
        NIIClIII=[]
        
    class Ne_fake:
        NIISII=[]
        OISII=[]
        OIISII=[]
        OIIISII=[]
        SIIISII=[]
        OIIOII=[]
        NIIOII=[]
        OIOII=[]
        OIIIClIII=[]
        SIIIClIII=[]
        OIIIArVI=[]
        SIIIArVI=[]
        NIIClIII=[]    
    
    class ion_abun:
        HeIa=[]
        HeIb=[]
        HeIIa=[]
        HeIIb=[]
        NI=[]
        NIIa=[]
        NIIb=[]
        NIIc=[]
        OIa=[]
        OIb=[]
        OIc=[]
        OIIa=[]
        OIIb=[]
        OIIc=[]
        OIId=[]
        OIIIa=[]
        OIIIb=[]
        OIIIc=[]
        SIIa=[]
        SIIb=[]
        SIIIa=[]
        SIIIb=[]
        NeIIIa=[]
        NeIIIb=[]
        ArIII=[]
        ArIVa=[]
        ArIVb=[]
        ClIIIa=[]
        ClIIIb=[]
    
    class ion_abun_error:
        HeIa=[]
        HeIb=[]
        HeIIa=[]
        HeIIb=[]
        NI=[]
        NIIa=[]
        NIIb=[]
        NIIc=[]
        OIa=[]
        OIb=[]
        OIc=[]
        OIIa=[]
        OIIb=[]
        OIIc=[]
        OIId=[]
        OIIIa=[]
        OIIIb=[]
        OIIIc=[]
        SIIa=[]
        SIIb=[]
        SIIIa=[]
        SIIIb=[]
        NeIIIa=[]
        NeIIIb=[]
        ArIII=[]
        ArIVa=[]
        ArIVb=[]
        ClIIIa=[]
        ClIIIb=[]
    
    class elem_abun_KB:
        He=[]
        N=[]
        O=[]
        S=[]
        Ne=[]
        Ar=[]
        Cl=[]
    
    class elem_abun_KB_ratio:
        NO=[]
        SO=[]
        NeO=[]
        ArO=[]
        ClO=[]
    
    class elem_ICF_KB:
        He=[]
        N=[]
        O=[]
        S=[]
        Ne=[]
        Ar=[]
        Cl=[]
        
    class elem_abun_DIMS:
        He=[]
        N=[]
        O=[]
        S=[]
        Ne=[]
        Ar=[]
        Cl=[]
    
    class elem_abun_DIMS_ratio:
        NO=[]
        SO=[]
        NeO=[]
        ArO=[]
        ClO=[]
    
    class elem_ICF_DIMS:
        He=[]
        N=[]
        O=[]
        S=[]
        Ne=[]
        Ar=[]
        Cl=[]
        
        
    class elem_abun_KB_error:
        He=[]
        N=[]
        O=[]
        S=[]
        Ne=[]
        Ar=[]
        Cl=[]
    
    class elem_abun_KB_ratio_error:
        NO=[]
        SO=[]
        NeO=[]
        ArO=[]
        ClO=[]
    
    class elem_ICF_KB_error:
        He=[]
        N=[]
        O=[]
        S=[]
        Ne=[]
        Ar=[]
        Cl=[]
        
    class elem_abun_DIMS_error:
        He=[]
        N=[]
        O=[]
        S=[]
        Ne=[]
        Ar=[]
        Cl=[]
    
    class elem_abun_DIMS_ratio_error:
        NO=[]
        SO=[]
        NeO=[]
        ArO=[]
        ClO=[]
    
    class elem_ICF_DIMS_error:
        He=[]
        N=[]
        O=[]
        S=[]
        Ne=[]
        Ar=[]
        Cl=[]    
        
    ####################
    ### for all the angle in the rotation angle task###
    ####################
    pn.log_.open_file('my_log.txt')
    file10 = open('general_output_file.txt','a')

    
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

    
    file6 = open('PyNeb_output_total_abund_ICFs_per_angles.txt','w')
    print("ang",'He      ', 'N_KB    ','ICF(N)_KB','N_DI    ','ICF(N)_DI', 'O_KB    ','ICF(O)_KB', 'O_DI    ','ICF(O)_DI','S_KB    ', 'ICF(S)_KB','S_DI    ','ICF(S)_DI','Ne_KB  ','ICF(Ne)_KB','Ne_DI  ','ICF(Ne)_DI','Cl_KB  ','ICF(Cl)_KB','Cl_DI  ','ICF(Cl)_DI','Ar_KB ','ICF(Ar)_KB','Ar_DI   ','ICF(Ar)_DI        ','NO_ratio_KB','NO_ratio_DIMS','SO_ratio_KB','SO_ratio_DIMS','NeO_ratio_KB','NeO_ratio_DIMS','ClO_ratio_KB','ClO_ratio_DIMS','ArO_ratio_KB','ArO_ratio_DIMS', file=file6)      
    
    file5 = open('PyNeb_output_ionic_abund_per_angles.txt','w')
    print("ang",'He+(5876)','He+(6678)','He++(4686)','He++(5412)','N0(5199)','N+(5755)','N+(6548)','N+(6584)','O0(5577)','O0(6300)','O0(6363)','O+(3727)','O+(3729)','O+(7320)','O+(7330)','O+++(4363)','O+++(4959)','O+++(5007)','S+(6716)  ','S+(6731)','S++(6312)','S++(9069)','Ne++(3868)','Ne++(3968)','Cl++(5518)', 'Cl++(5538)','Ar++(7136)','Ar+++(4712)', 'Ar+++(4740)', file=file5)  

    
    file4 = open('PyNeb_output_Te_and_Ne_per_angles.txt','w')
    print("ang","Te_NII","Ne_SII ","Te_OI","Ne_SII","Te_OII","Ne_SII","Te_OIII","Ne_SII","Te_SIII","Ne_SII","Te_OII","Ne_OII","Te_NII","Ne_OII","Te_OI ","Ne_OII","Te_OIII","Ne_ClIII","Te_SIII","Ne_ClIII","Te_OIII","Ne_ArVI","Te_SIII","Ne_ArVI","Te_NII","Ne_ClIII", file=file4)
    
    for i in range(0,len(ang)):
        Te.NIISII.append(0),Ne.NIISII.append(0),Te.OISII.append(0),Ne.OISII.append(0),Te.OIISII.append(0),Ne.OIISII.append(0),Te.OIIISII.append(0),Ne.OIIISII.append(0),Te.SIIISII.append(0),Ne.SIIISII.append(0),Te.OIIOII.append(0),Ne.OIIOII.append(0),Te.NIIOII.append(0),Ne.NIIOII.append(0),Te.OIOII.append(0),Ne.OIOII.append(0),Te.OIIIClIII.append(0),Ne.OIIIClIII.append(0),Te.SIIIClIII.append(0),Ne.SIIIClIII.append(0),Te.OIIIArVI.append(0),Ne.OIIIArVI.append(0),Te.SIIIArVI.append(0),Ne.SIIIArVI.append(0),Te.NIIClIII.append(0),Ne.NIIClIII.append(0)
        
        Te_error.NIISII.append(0),Ne_error.NIISII.append(0),Te_error.OISII.append(0),Ne_error.OISII.append(0),Te_error.OIISII.append(0),Ne_error.OIISII.append(0),Te_error.OIIISII.append(0),Ne_error.OIIISII.append(0),Te_error.SIIISII.append(0),Ne_error.SIIISII.append(0),Te_error.OIIOII.append(0),Ne_error.OIIOII.append(0),Te_error.NIIOII.append(0),Ne_error.NIIOII.append(0),Te_error.OIOII.append(0),Ne_error.OIOII.append(0),Te_error.OIIIClIII.append(0),Ne_error.OIIIClIII.append(0),Te_error.SIIIClIII.append(0),Ne_error.SIIIClIII.append(0),Te_error.OIIIArVI.append(0),Ne_error.OIIIArVI.append(0),Te_error.SIIIArVI.append(0),Ne_error.SIIIArVI.append(0),Te_error.NIIClIII.append(0),Ne_error.NIIClIII.append(0)
        
        print("######################################"," slit PA  ", ang[i] ," ##############################", file=file4)
        print("######################################"," slit PA  ", ang[i] ," ##############################", file=file5)
        print("######################################"," slit PA  ", ang[i] ," ##############################", file=file6)
        
            
        file3 = open('test.dat','w')
        print("LINE test err", file=file3)
    
        error=nfes.flux_Error(flux_angles.Ha_6563[i],flux_angles_error.Ha_6563[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("H1r_6563A", '{:03.3f}'.format(flux_angles.Ha_6563[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "H1r_6563Ae", '{:03.3f}'.format(error)

        error=nfes.flux_Error(flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("H1r_4861A", '{:03.3f}'.format(flux_angles.Hb_4861[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "H1r_4861Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.Hg_4340[i],flux_angles_error.Hg_4340[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("H1r_4341A", '{:03.3f}'.format(flux_angles.Hg_4340[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "H1r_4341Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.Hd_4101[i],flux_angles_error.Hd_4101[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("H1r_4102A", '{:03.3f}'.format(flux_angles.Hd_4101[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "H1r_4102Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.HeIIa_4686[i],flux_angles_error.HeIIa_4686[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("He2r_4686A", '{:03.3f}'.format(flux_angles.HeIIa_4686[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "He2r_4686Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.HeIIb_5412[i],flux_angles_error.HeIIb_5412[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("He2r_5411A", '{:03.3f}'.format(flux_angles.HeIIb_5412[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "He2r_5411Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.HeIa_5876[i],flux_angles_error.HeIa_5876[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("He1r_5876A", '{:03.3f}'.format(flux_angles.HeIa_5876[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "He1r_5876Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.HeIb_6678[i],flux_angles_error.HeIb_6678[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("He1r_6678A", '{:03.3f}'.format(flux_angles.HeIb_6678[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "He1r_6678Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.SIIa_6716[i],flux_angles_error.SIIa_6716[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("S2_6716A", '{:03.3f}'.format(flux_angles.SIIa_6716[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "S2_6716Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.SIIb_6731[i],flux_angles_error.SIIb_6731[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("S2_6731A", '{:03.3f}'.format(flux_angles.SIIb_6731[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "S2_6731Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.SIIIa_6312[i],flux_angles_error.SIIIa_6312[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("S3_6312A", '{:03.3f}'.format(flux_angles.SIIIa_6312[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "S3_6312Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.SIIIb_9069[i],flux_angles_error.SIIIb_9069[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("S3_9069A", '{:03.3f}'.format(flux_angles.SIIIb_9069[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "S3_9069Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.NI_5199[i],flux_angles_error.NI_5199[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("N1_5200A", '{:03.3f}'.format(flux_angles.NI_5199[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "N1_5200Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.NIIa_5755[i],flux_angles_error.NIIa_5755[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("N2_5755A", '{:03.3f}'.format(flux_angles.NIIa_5755[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "N2_5755Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.NIIb_6548[i],flux_angles_error.NIIb_6548[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("N2_6548A", '{:03.3f}'.format(flux_angles.NIIb_6548[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "N2_6548Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.NIIc_6584[i],flux_angles_error.NIIc_6584[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("N2_6584A", '{:03.3f}'.format(flux_angles.NIIc_6584[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "N2_6584Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIa_5577[i],flux_angles_error.OIa_5577[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O1_5577A", '{:03.3f}'.format(flux_angles.OIa_5577[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O1_5577Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIb_6300[i],flux_angles_error.OIb_6300[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O1_6300A", '{:03.3f}'.format(flux_angles.OIb_6300[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O1_6300Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIc_6363[i],flux_angles_error.OIc_6363[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O1_6364A", '{:03.3f}'.format(flux_angles.OIc_6363[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O1_6364Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIIa_3727[i],flux_angles_error.OIIa_3727[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O2_3726A", '{:03.3f}'.format(flux_angles.OIIa_3727[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O2_3726Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIIb_3729[i],flux_angles_error.OIIb_3729[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O2_3729A", '{:03.3f}'.format(flux_angles.OIIb_3729[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O2_3729Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIIc_7320[i],flux_angles_error.OIIc_7320[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O2_7319A+", '{:03.3f}'.format(flux_angles.OIIc_7320[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O2_7319A+e", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIId_7330[i],flux_angles_error.OIId_7330[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O2_7330A+", '{:03.3f}'.format(flux_angles.OIId_7330[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O2_7330A+e", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIIIa_4363[i],flux_angles_error.OIIIa_4363[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O3_4363A", '{:03.3f}'.format(flux_angles.OIIIa_4363[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O3_4363Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIIIb_4959[i],flux_angles_error.OIIIb_4959[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O3_4959A", '{:03.3f}'.format(flux_angles.OIIIb_4959[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O3_4959Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.OIIIc_5007[i],flux_angles_error.OIIIc_5007[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("O3_5007A", '{:03.3f}'.format(flux_angles.OIIIc_5007[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "O3_5007Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.ClIIIa_5517[i],flux_angles_error.ClIIIa_5517[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])        
        print("Cl3_5518A", '{:03.3f}'.format(flux_angles.ClIIIa_5517[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "Cl3_5518Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.ClIIIb_5538[i],flux_angles_error.ClIIIb_5538[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("Cl3_5538A", '{:03.3f}'.format(flux_angles.ClIIIb_5538[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)        
#        print >>file3, "Cl3_5538Ae", '{:03.3f}'.format(error)        
        
        error=nfes.flux_Error(flux_angles.NeIIIa_3868[i],flux_angles_error.NeIIIa_3868[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("Ne3_3869A", '{:03.3f}'.format(flux_angles.NeIIIa_3868[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "Ne3_3869Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.NeIIIb_3967[i],flux_angles_error.NeIIIb_3967[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("Ne3_3968A", '{:03.3f}'.format(flux_angles.NeIIIb_3967[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)                
#        print >>file3, "Ne3_3968Ae", '{:03.3f}'.format(error)                
        
        error=nfes.flux_Error(flux_angles.ArIII_7136[i],flux_angles_error.ArIII_7136[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("Ar3_7136A", '{:03.3f}'.format(flux_angles.ArIII_7136[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "Ar3_7136Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.ArIVa_4712[i],flux_angles_error.ArIVa_4712[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("Ar4_4711A", '{:03.3f}'.format(flux_angles.ArIVa_4712[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "Ar4_4711Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.ArIVb_4740[i],flux_angles_error.ArIVb_4740[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("Ar4_4740A", '{:03.3f}'.format(flux_angles.ArIVb_4740[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "Ar4_4740Ae", '{:03.3f}'.format(error)
        
        error=nfes.flux_Error(flux_angles.CI_8727[i],flux_angles_error.CI_8727[i],flux_angles.Hb_4861[i],flux_angles_error.Hb_4861[i])
        print("C1_8728A", '{:03.3f}'.format(flux_angles.CI_8727[i]*100/flux_angles.Hb_4861[i]),'{:03.3f}'.format(error), file=file3)
#        print >>file3, "C1_8728Ae", '{:03.3f}'.format(error)
#        print >>file3, "C2_????A", '{:03.3f}'.format(flux_angles.CII_6461[i]*100/flux_angles.Hb_4861[i])
        
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
        obs_fake.readData(obs_data, fileFormat='lines_in_rows_err_cols', errIsRelative=False)#, err_default=0.05)
        obs.readData(obs_data, fileFormat='lines_in_rows_err_cols', errIsRelative=False)#, err_default=0.05)
        
        # Add a number of "fake" observations aroung the right one, to make Monte Carlo statistics
        index_MC=param_mod_name.index("number_fake_obs_MC")
        obs_fake.addMonteCarloObs(N = param_model_values[index_MC])

        obs.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)
        obs_fake.def_EBV(label1="H1r_6563A", label2="H1r_4861A", r_theo=2.85)

        print("######################################################################################", file=file10)
        print("######################################################################################", file=file10)
#        print , "PyNeb: cH(b)=",obs.extinction.cHbeta[0],"+-",obs.extinction.cHbeta[1]," angle=",ang[i]
        print("PyNeb: cH(b) = {:.2f} +/- {:.2f}".format(obs.extinction.cHbeta[0],np.std(obs_fake.extinction.cHbeta)), file=file10)
        print("PyNeb: E(B-V) = {:.2f} +/- {:.2f}".format(obs.extinction.E_BV[0],np.std(obs_fake.extinction.E_BV)), file=file10)
#        print >>file10, "PyNeb: E(B-V)=", obs.extinction.E_BV[0],"+-",obs.extinction.E_BV[1]," angle=",ang[i]
    

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
#        print(obs.getIntens()['H1r_6563A'])
#        obs.printIntens(returnObs=True)
#        print obs_fake.printIntens()
        
        ########################
        ### Density diagnostics
        ########################
        index_line1=line_names.index("S2_6716s")
        index_line2=line_names.index("S2_6731s")
        index_Ne_SII="no"
        if lines_available[index_line1]=="yes" and  lines_available[index_line2]=="yes":
                diags.addDiag([
                    '[SII] 6731/6716',
                    ])    
                diags.addClabel('[SII] 6731/6716', '[SII]a')
                index_Ne_SII="yes"

        index_line1=line_names.index("O2_3727s")
        index_line2=line_names.index("O2_3729s")
        index_Ne_OII="no"
        if lines_available[index_line1]=="yes" and  lines_available[index_line2]=="yes":
                diags.addDiag([
                    '[OII] 3726/3729',
                    ])    
                index_Ne_OII="yes"

        index_line1=line_names.index("Cl3_5517s")
        index_line2=line_names.index("Cl3_5538s")
        index_Ne_ClIII="no"
        if lines_available[index_line1]=="yes" and  lines_available[index_line2]=="yes":
                diags.addDiag([
                    '[ClIII] 5538/5518',
                    ])    
                index_Ne_ClIII="yes"

        index_line1=line_names.index("Ar4_4712s")
        index_line2=line_names.index("Ar4_4740s")
        index_Ne_ArVI="no"
        if lines_available[index_line1]=="yes" and  lines_available[index_line2]=="yes":
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
        if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes":
            diags.addDiag([
                    '[NII] 5755/6584+',
                    ])    
            index_Te_NII="yes"

        index_line1=line_names.index("O3_4363s")
        index_line2=line_names.index("O3_4959s")
        index_line3=line_names.index("O3_5007s")
        index_Te_OIII="no"
        if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes":
            diags.addDiag([
                '[OIII] 4363/5007+',
                ])      
            index_Te_OIII="yes"

        index_line1=line_names.index("O1_6300s")
        index_line2=line_names.index("O1_6363s")
        index_line3=line_names.index("O1_5577s")
        index_Te_OI="no"
        index_Te_OI2="no"
        if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes":
            diags.addDiag([
                '[OI] 5577/6300+',
                ]) 
            index_Te_OI="yes"
        elif lines_available[index_line1]=="yes" and lines_available[index_line2]=="no" and lines_available[index_line3]=="yes":
            diags.addDiag([
                '[OI] 5577/6300',
                ])
            index_Te_OI2="yes"

        index_line1=line_names.index("O2_3727s")
        index_line2=line_names.index("O2_3729s")
        index_line3=line_names.index("O2_7320s")
        index_line4=line_names.index("O2_7330s")
        index_Te_OII="no"
        if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes" and lines_available[index_line4]=="yes":
            diags.addDiag([
                '[OII] 3727+/7325+',
                ])
            index_Te_OII="yes"
    
        index_line1=line_names.index("S3_6312s")
        index_line2=line_names.index("S3_9069s")
        #index_line3=line_names.index("S3_9592s")
        index_Te_SIII="no"
        if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes":
            diags.addDiag([
                '[SIII] 6312/9069',
                ])
            index_Te_SIII="yes"
#       elif lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes" and lines_available[index_line3]=="yes":
#           diags.addDiag([
#                '[SIII] 6312/9200+',
#           ])  
#           index_Te_SIII="yes"



#        # Create the emission maps to be compared to the observation data (some overkill here)
#        emisgrids = pn.getEmisGridDict(atom_list=diags.getUniqueAtoms(), den_max=1e6)

        ### Plot
        # Create the contour plot as the intersection of tem-den emission maps with dereddened line ratios
#        diags.plot(emisgrids, obs)

        # Place the title
#       plt.title(title)

        # Display the plot
#        plt.show()

        print("######################################################################################")
        #The observed ratio can be automatically extracted from an Observation object named obs:
        if (index_Te_NII=="yes" and index_Ne_SII=="yes"):
            a, b = diags.getCrossTemDen('[NII] 5755/6584+', '[SII] 6731/6716', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+', '[SII] 6731/6716', obs=obs)
            Te_fake.NIISII=a
            Ne_fake.NIISII=b
            
            Te.NIISII[i]=a1
            Ne.NIISII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.NIISII[i]=np.std(np.nan_to_num(a))
            Ne_error.NIISII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: low zones (NII-SII): Te_NII  = {0:5.0f}, Ne_SII = {1:7.1f}  cm-1'.format(Te.NIISII[i], Ne.NIISII[i])," angle=",ang[i], file=file10)
    
        if (index_Te_OI=="yes" and index_Te_OI2=="no" and index_Ne_SII=="yes"):
            a, b = diags.getCrossTemDen('[OI] 5577/6300+', '[SII] 6731/6716', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OI] 5577/6300+', '[SII] 6731/6716', obs=obs)
            Te_fake.OISII=a
            Ne_fake.OISII=b
            
            Te.OISII[i]=a1
            Ne.OISII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OISII[i]=np.std(np.nan_to_num(a))
            Ne_error.OISII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: low zones (OI-SII): Te_OI  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'.format(Te.OISII[i],Ne.OISII[i])," angle=",ang[i], file=file10) 
        
        if (index_Te_OI=="no" and index_Te_OI2=="yes" and index_Ne_SII=="yes"):
            a, b = diags.getCrossTemDen('[OI] 5577/6300', '[SII] 6731/6716', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OI] 5577/6300', '[SII] 6731/6716', obs=obs)
            Te_fake.OISII=a
            Ne_fake.OISII=b
            
            Te.OISII[i]=a1
            Ne.OISII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OISII[i]=np.std(np.nan_to_num(a))
            Ne_error.OISII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: low zones (OI-SII): Te_OI  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'.format(Te.OISII[i],Ne.OISII[i])," angle=",ang[i], file=file10)
    
        if (index_Te_OII=="yes" and index_Ne_SII=="yes"):
            a, b = diags.getCrossTemDen('[OII] 3727+/7325+', '[SII] 6731/6716', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OII] 3727+/7325+', '[SII] 6731/6716', obs=obs)
            Te_fake.OIISII=a
            Ne_fake.OIISII=b
            
            Te.OIISII[i]=a1
            Ne.OIISII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OIISII[i]=np.std(np.nan_to_num(a))
            Ne_error.OIISII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: low zones (OII-SII): Te_OII  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'.format(Te.OIISII[i], Ne.OIISII[i])," angle=",ang[i], file=file10)
    
        if (index_Te_OIII=="yes" and index_Ne_SII=="yes"):
            a, b = diags.getCrossTemDen('[OIII] 4363/5007+', '[SII] 6731/6716', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+', '[SII] 6731/6716', obs=obs)
            Te_fake.OIIISII=a
            Ne_fake.OIIISII=b
            
            Te.OIIISII[i]=a1
            Ne.OIIISII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OIIISII[i]=np.std(np.nan_to_num(a))
            Ne_error.OIIISII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: medium zones (OIII-SII): Te_OIII  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'.format(Te.OIIISII[i], Ne.OIIISII[i])," angle=",ang[i], file=file10)

        if (index_Te_SIII=="yes" and index_Ne_SII=="yes"):
            a, b = diags.getCrossTemDen('[SIII] 6312/9069', '[SII] 6731/6716', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069', '[SII] 6731/6716', obs=obs)
            Te_fake.SIIISII=a
            Ne_fake.SIIISII=b
            
            Te.SIIISII[i]=a1
            Ne.SIIISII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.SIIISII[i]=np.std(np.nan_to_num(a))
            Ne_error.SIIISII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: medium zones (SIII-SII): Te_SIII  = {0:5.0f},Ne_SII = {1:7.1f}  cm-1'.format(Te.SIIISII[i], Ne.SIIISII[i])," angle=",ang[i], file=file10)  
        
        if (index_Te_OII=="yes" and index_Ne_OII=="yes"):
            a, b = diags.getCrossTemDen('[OII] 3727+/7325+', '[OII] 3726/3729', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OII] 3727+/7325+', '[OII] 3726/3729', obs=obs)
            Te_fake.OIIOII=a
            Ne_fake.OIIOII=b
            
            Te.OIIOII[i]=a1
            Ne.OIIOII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OIIOII[i]=np.std(np.nan_to_num(a))
            Ne_error.OIIOII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: low zones (OII-OII): Te_OII  = {0:5.0f},Ne_OII = {1:7.1f}  cm-1'.format(Te.OIIOII[i], Ne.OIIOII[i])," angle=",ang[i], file=file10)
    
        if (index_Te_NII=="yes" and index_Ne_OII=="yes"):
            a, b = diags.getCrossTemDen('[NII] 5755/6584+', '[OII] 3726/3729', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+', '[OII] 3726/3729', obs=obs)
            Te_fake.NIIOII=a
            Ne_fake.NIIOII=b
            
            Te.NIIOII[i]=a1
            Ne.NIIOII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.NIIOII[i]=np.std(np.nan_to_num(a))
            Ne_error.NIIOII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: low zones (NII-OII): Te_NII  = {0:5.0f},Ne_OII = {1:7.1f}  cm-1'.format(Te.NIIOII[i],Ne.NIIOII[i])," angle=",ang[i], file=file10)
                  
        if (index_Te_OI=="yes" and index_Te_OI2=="no" and index_Ne_OII=="yes"):
            a, b = diags.getCrossTemDen('[OI] 5577/6300+', '[OII] 3726/3729', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OI] 5577/6300+', '[OII] 3726/3729', obs=obs)
            Te_fake.OIOII=a
            Ne_fake.OIOII=b
                        
            Te.OIOII[i]=a1
            Ne.OIOII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OIOII[i]=np.std(np.nan_to_num(a))
            Ne_error.OIOII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: low zones (OI-OII): Te_OI  = {0:5.0f},Ne_OII = {1:7.1f}  cm-1'.format(Te.OIOII[i], Ne.OIOII[i])," angle=",ang[i], file=file10)
        
        if (index_Te_OI=="no" and index_Te_OI2=="yes" and index_Ne_OII=="yes"):
            a, b = diags.getCrossTemDen('[OI] 5577/6300', '[OII] 3726/3729', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OI] 5577/6300', '[OII] 3726/3729', obs=obs)
            Te_fake.OIOII=a
            Ne_fake.OIOII=b
                        
            Te.OIOII[i]=a1
            Ne.OIOII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OIOII[i]=np.std(np.nan_to_num(a))
            Ne_error.OIOII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: low zones (OI-OII): Te_OI  = {0:5.0f},Ne_OII = {1:7.1f}  cm-1'.format(Te.OIOII[i], Ne.OIOII[i])," angle=",ang[i], file=file10)
    
        if (index_Te_OIII=="yes" and index_Ne_ClIII=="yes"):
            a, b = diags.getCrossTemDen('[OIII] 4363/5007+', '[ClIII] 5538/5518', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+', '[ClIII] 5538/5518', obs=obs)
            Te_fake.OIIIClIII=a
            Ne_fake.OIIIClIII=b
            
            Te.OIIIClIII[i]=a1
            Ne.OIIIClIII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OIIIClIII[i]=np.std(np.nan_to_num(a))
            Ne_error.OIIIClIII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: medium zones (OIII-ClIII): Te_OIII  = {0:5.0f},Ne_ClIII = {1:7.1f}  cm-1'.format(Te.OIIIClIII[i],Ne.OIIIClIII[i])," angle=",ang[i], file=file10)          
                  
        if (index_Te_SIII=="yes" and index_Ne_ClIII=="yes"):
            a, b = diags.getCrossTemDen('[SIII] 6312/9069', '[ClIII] 5538/5518', obs=obs_fake)         
            a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069', '[ClIII] 5538/5518', obs=obs)
            Te_fake.SIIIClIII=a
            Ne_fake.SIIIClIII=b
            
            Te.SIIIClIII[i]=a1
            Ne.SIIIClIII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.SIIIClIII[i]=np.std(np.nan_to_num(a))
            Ne_error.SIIIClIII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: medium zones (SIII-ClIII): Te_SIII  = {0:5.0f},Ne_ClIII = {1:7.1f}  cm-1'.format(Te.SIIIClIII[i], Ne.SIIIClIII[i])," angle=",ang[i], file=file10)
    
        if (index_Te_OIII=="yes" and index_Ne_ArVI=="yes"):
            a, b = diags.getCrossTemDen('[OIII] 4363/5007+', '[ArIV] 4740/4711', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[OIII] 4363/5007+', '[ArIV] 4740/4711', obs=obs)
            Te_fake.OIIIArVI=a
            Ne_fake.OIIIArVI=b
            
            Te.OIIIArVI[i]=a1
            Ne.OIIIArVI[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.OIIIArVI[i]=np.std(np.nan_to_num(a))
            Ne_error.OIIIArVI[i]=np.std(np.nan_to_num(b))
            print('PyNeb: medium zones (OIII-ArVI): Te_OIII  = {0:5.0f},Ne_ArVI = {1:7.1f}  cm-1'.format(Te.OIIIArVI[i], Ne.OIIIArVI[i])," angle=",ang[i], file=file10)          
        
        if (index_Te_SIII=="yes" and index_Ne_ArVI=="yes"):
            a, b = diags.getCrossTemDen('[SIII] 6312/9069', '[ArIV] 4740/4711', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[SIII] 6312/9069', '[ArIV] 4740/4711', obs=obs)
            Te_fake.NIISII=a
            Ne_fake.NIISII=b
            
            Te.SIIIArVI[i]=a1
            Ne.SIIIArVI[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.SIIIArVI[i]=np.std(np.nan_to_num(a))
            Ne_error.SIIIArVI[i]=np.std(np.nan_to_num(b))
            print('PyNeb: medium zones (SIII-ArVI): Te_SIII  = {0:5.0f},Ne_ArVI = {1:7.1f}  cm-1'.format(Te.SIIIArVI[i],Ne.SIIIArVI[i])," angle=",ang[i], file=file10)
       
        if (index_Te_NII=="yes" and index_Ne_ClIII=="yes"):
            a, b = diags.getCrossTemDen('[NII] 5755/6584+', '[ClIII] 5538/5518', obs=obs_fake)
            a1, b1 = diags.getCrossTemDen('[NII] 5755/6584+', '[ClIII] 5538/5518', obs=obs)
            Te_fake.NIIClIII=a
            Ne_fake.NIIClIII=b
            
            Te.NIIClIII[i]=a1
            Ne.NIIClIII[i]=b1
            
            a[a == inf] =  0
            a=np.nan_to_num(a)
            a = [i2 for i2 in a if i2 != 0]
            a = [i2 for i2 in a if i2 < a1*3 and i2>a1*0.3]

            b[b == inf] =  0
            b=np.nan_to_num(b)
            b = [i2 for i2 in b if i2 != 0]
            b = [i2 for i2 in b if i2 < b1*3 and i2>b1*0.3]
            
            Te_error.NIIClIII[i]=np.std(np.nan_to_num(a))
            Ne_error.NIIClIII[i]=np.std(np.nan_to_num(b))
            print('PyNeb: medium zones (NII-ClIII): Te_NII  = {0:5.0f},Ne_ClIII = {1:7.1f}  cm-1'.format(Te.NIIClIII[i],Ne.NIIClIII[i])," angle=",ang[i], file=file10)

                  
        print('{:03.0f}'.format(ang[i]), '{:06.0f}'.format(Te.NIISII[i]), '{:06.0f}'.format(Ne.NIISII[i]),'{:06.0f}'.format(Te.OISII[i]),'{:06.0f}'.format(Ne.OISII[i]),'{:06.0f}'.format(Te.OIISII[i]),'{:06.0f}'.format(Ne.OIISII[i]),'{:06.0f}'.format(Te.OIIISII[i]),'{:07.0f}'.format(Ne.OIIISII[i]),'{:07.0f}'.format(Te.SIIISII[i]),'{:06.0f}'.format(Ne.SIIISII[i]),'{:06.0f}'.format(Te.OIIOII[i]),'{:06.0f}'.format(Ne.OIIOII[i]),'{:06.0f}'.format(Te.NIIOII[i]),'{:06.0f}'.format(Ne.NIIOII[i]),'{:06.0f}'.format(Te.OIOII[i]),'{:06.0f}'.format(Ne.OIOII[i]),'{:07.0f}'.format(Te.OIIIClIII[i]),'{:08.0f}'.format(Ne.OIIIClIII[i]),'{:07.0f}'.format(Te.SIIIClIII[i]),'{:08.0f}'.format(Ne.SIIIClIII[i]),'{:07.0f}'.format(Te.OIIIArVI[i]),'{:07.0f}'.format(Ne.OIIIArVI[i]),'{:07.0f}'.format(Te.SIIIArVI[i]),'{:07.0f}'.format(Ne.SIIIArVI[i]),'{:07.0f}'.format(Te.NIIClIII[i]),'{:07.0f}'.format(Ne.NIIClIII[i]), file=file4)

        print('{:03.0f}'.format(ang[i]), '{:06.0f}'.format(Te_error.NIISII[i]), '{:06.0f}'.format(Ne_error.NIISII[i]),'{:06.0f}'.format(Te_error.OISII[i]),'{:06.0f}'.format(Ne_error.OISII[i]),'{:06.0f}'.format(Te_error.OIISII[i]),'{:06.0f}'.format(Ne_error.OIISII[i]),'{:06.0f}'.format(Te_error.OIIISII[i]),'{:07.0f}'.format(Ne_error.OIIISII[i]),'{:07.0f}'.format(Te_error.SIIISII[i]),'{:06.0f}'.format(Ne_error.SIIISII[i]),'{:06.0f}'.format(Te_error.OIIOII[i]),'{:06.0f}'.format(Ne_error.OIIOII[i]),'{:06.0f}'.format(Te_error.NIIOII[i]),'{:06.0f}'.format(Ne_error.NIIOII[i]),'{:06.0f}'.format(Te_error.OIOII[i]),'{:06.0f}'.format(Ne_error.OIOII[i]),'{:07.0f}'.format(Te_error.OIIIClIII[i]),'{:08.0f}'.format(Ne_error.OIIIClIII[i]),'{:07.0f}'.format(Te_error.SIIIClIII[i]),'{:08.0f}'.format(Ne_error.SIIIClIII[i]),'{:07.0f}'.format(Te_error.OIIIArVI[i]),'{:07.0f}'.format(Ne_error.OIIIArVI[i]),'{:07.0f}'.format(Te_error.SIIIArVI[i]),'{:07.0f}'.format(Ne_error.SIIIArVI[i]),'{:07.0f}'.format(Te_error.NIIClIII[i]),'{:07.0f}'.format(Ne_error.NIIClIII[i]), file=file4)
        
        ###########################
        ### calculate, print and save ionic abundances
        ##########################
        ion_abun,ion_abun_error,atom_abun_fake=iaas.ionicabundances(obs,obs_fake,Te,Te_error,Te_fake,Ne,Ne_error,Ne_fake,i,ang[i],line_names,lines_available,ion_abun,ion_abun_error,param_estimated,param_requered,param_mod_name)
        
        print('{:03.0f}'.format(ang[i]),'{:5.3e}'.format(ion_abun.HeIa[i]), '{:5.3e}'.format(ion_abun.HeIb[i]),'{:5.4e}'.format(ion_abun.HeIIa[i]),'{:5.4e}'.format(ion_abun.HeIIb[i]),'{:5.2e}'.format(ion_abun.NI[i]),'{:5.2e}'.format(ion_abun.NIIa[i]),'{:5.2e}'.format(ion_abun.NIIb[i]),'{:5.2e}'.format(ion_abun.NIIc[i]),'{:5.2e}'.format(ion_abun.OIa[i]),'{:5.2e}'.format(ion_abun.OIb[i]),'{:5.2e}'.format(ion_abun.OIc[i]),'{:5.2e}'.format(ion_abun.OIIa[i]),'{:5.2e}'.format(ion_abun.OIIb[i]),'{:5.2e}'.format(ion_abun.OIIc[i]),'{:5.2e}'.format(ion_abun.OIId[i]),'{:5.4e}'.format(ion_abun.OIIIa[i]),'{:5.4e}'.format(ion_abun.OIIIb[i]),'{:5.4e}'.format(ion_abun.OIIIc[i]),'{:5.2e}'.format(ion_abun.SIIa[i]),'{:5.2e}'.format(ion_abun.SIIb[i]),'{:5.3e}'.format(ion_abun.SIIIa[i]),'{:5.3e}'.format(ion_abun.SIIIb[i]),'{:5.4e}'.format(ion_abun.NeIIIa[i]),'{:5.4e}'.format(ion_abun.NeIIIb[i]),'{:5.4e}'.format(ion_abun.ClIIIa[i]),'{:5.4e}'.format(ion_abun.ClIIIb[i]),'{:5.4e}'.format(ion_abun.ArIII[i]),'{:5.5e}'.format(ion_abun.ArIVa[i]),'{:5.5e}'.format(ion_abun.ArIVb[i]), file=file5)
 
        print('{:03.0f}'.format(ang[i]),'{:5.3e}'.format(ion_abun_error.HeIa[i]), '{:5.3e}'.format(ion_abun_error.HeIb[i]),'{:5.4e}'.format(ion_abun_error.HeIIa[i]),'{:5.4e}'.format(ion_abun_error.HeIIb[i]),'{:5.2e}'.format(ion_abun_error.NI[i]),'{:5.2e}'.format(ion_abun_error.NIIa[i]),'{:5.2e}'.format(ion_abun_error.NIIb[i]),'{:5.2e}'.format(ion_abun_error.NIIc[i]),'{:5.2e}'.format(ion_abun_error.OIa[i]),'{:5.2e}'.format(ion_abun_error.OIb[i]),'{:5.2e}'.format(ion_abun_error.OIc[i]),'{:5.2e}'.format(ion_abun_error.OIIa[i]),'{:5.2e}'.format(ion_abun_error.OIIb[i]),'{:5.2e}'.format(ion_abun_error.OIIc[i]),'{:5.2e}'.format(ion_abun_error.OIId[i]),'{:5.4e}'.format(ion_abun_error.OIIIa[i]),'{:5.4e}'.format(ion_abun_error.OIIIb[i]),'{:5.4e}'.format(ion_abun_error.OIIIc[i]),'{:5.2e}'.format(ion_abun_error.SIIa[i]),'{:5.2e}'.format(ion_abun_error.SIIb[i]),'{:5.3e}'.format(ion_abun_error.SIIIa[i]),'{:5.3e}'.format(ion_abun_error.SIIIb[i]),'{:5.4e}'.format(ion_abun_error.NeIIIa[i]),'{:5.4e}'.format(ion_abun_error.NeIIIb[i]),'{:5.4e}'.format(ion_abun_error.ClIIIa[i]),'{:5.4e}'.format(ion_abun_error.ClIIIb[i]),'{:5.4e}'.format(ion_abun_error.ArIII[i]),'{:5.5e}'.format(ion_abun_error.ArIVa[i]),'{:5.5e}'.format(ion_abun_error.ArIVb[i]), file=file5)
        
        ###########################
        ### calculate, print and save total elemental abundances and ICFs
        ##########################
        elem_abun_KB,elem_ICF_KB,elem_abun_DIMS,elem_ICF_DIMS,elem_abun_KB_error,elem_ICF_KB_error,elem_abun_DIMS_error,elem_ICF_DIMS_error,elem_abun_KB_ratio,elem_abun_KB_ratio_error,elem_abun_DIMS_ratio,elem_abun_DIMS_ratio_error=eaIas.total_abundances_ICFs(ion_abun,ion_abun_error,atom_abun_fake,i,line_names,lines_available,elem_abun_KB,elem_ICF_KB,elem_abun_DIMS,elem_ICF_DIMS,elem_abun_KB_error,elem_ICF_KB_error,elem_abun_DIMS_error,elem_ICF_DIMS_error,elem_abun_KB_ratio,elem_abun_KB_ratio_error,elem_abun_DIMS_ratio,elem_abun_DIMS_ratio_error,param_mod_name)
        
        print('{:03.0f}'.format(ang[i]),'{:5.2e}'.format(elem_abun_KB.He[i]), '{:5.2e}'.format(elem_abun_KB.N[i]),'{:09.5f}'.format(elem_ICF_KB.N[i]),'{:5.2e}'.format(elem_abun_DIMS.N[i]),'{:09.5f}'.format(elem_ICF_DIMS.N[i]   ),'{:5.2e}'.format(elem_abun_KB.O[i]),'{:09.5f}'.format(elem_ICF_KB.O[i]),'{:5.2e}'.format(elem_abun_DIMS.O[i]),'{:09.5f}'.format(elem_ICF_DIMS.O[i]),'{:5.2e}'.format(elem_abun_KB.S[i]),'{:09.5f}'.format(elem_ICF_KB.S[i]),'{:5.2e}'.format(elem_abun_DIMS.S[i]),'{:09.5f}'.format(elem_ICF_DIMS.S[i]),'{:5.2e}'.format(elem_abun_KB.Ne[i]),'{:09.5f}'.format(elem_ICF_KB.Ne[i]),'{:5.2e}'.format(elem_abun_DIMS.Ne[i]),'{:09.5f}'.format(elem_ICF_DIMS.Ne[i]),'{:5.2e}'.format(elem_abun_KB.Cl[i]),'{:09.5f}'.format(elem_ICF_KB.Cl[i]),'{:5.2e}'.format(elem_abun_DIMS.Cl[i]),'{:09.5f}'.format(elem_ICF_DIMS.Cl[i]),'{:5.2e}'.format(elem_abun_KB.Ar[i]),'{:09.5f}'.format(elem_ICF_KB.Ar[i]),'{:5.2e}'.format(elem_abun_DIMS.Ar[i]),'{:09.5f}'.format(elem_ICF_DIMS.Ar[i]),'{:9.2e}'.format(elem_abun_KB_ratio.NO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio.NO[i]),'{:9.2e}'.format(elem_abun_KB_ratio.SO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio.SO[i]), '{:9.2e}'.format(elem_abun_KB_ratio.NeO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio.NeO[i]),'{:9.2e}'.format(elem_abun_KB_ratio.ClO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio.ClO[i]),'{:9.2e}'.format(elem_abun_KB_ratio.ClO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio.ClO[i]), file=file6)
 
        print('{:03.0f}'.format(ang[i]), '{:5.2e}'.format(elem_abun_KB_error.He[i]), '{:5.2e}'.format(elem_abun_KB_error.N[i]),'{:09.5f}'.format(elem_ICF_KB_error.N[i]),'{:5.2e}'.format(elem_abun_DIMS_error.N[i]),'{:09.5f}'.format(elem_ICF_DIMS_error.N[i]),'{:5.2e}'.format(elem_abun_KB_error.O[i]),'{:09.5f}'.format(elem_ICF_KB_error.O[i]),'{:5.2e}'.format(elem_abun_DIMS_error.O[i]),'{:09.5f}'.format(elem_ICF_DIMS_error.O[i]),'{:5.2e}'.format(elem_abun_KB_error.S[i]),'{:09.5f}'.format(elem_ICF_KB_error.S[i]),'{:5.2e}'.format(elem_abun_DIMS_error.S[i]),'{:09.5f}'.format(elem_ICF_DIMS_error.S[i]),'{:5.2e}'.format(elem_abun_KB_error.Ne[i]),'{:09.5f}'.format(elem_ICF_KB_error.Ne[i]),'{:5.2e}'.format(elem_abun_DIMS_error.Ne[i]),'{:09.5f}'.format(elem_ICF_DIMS_error.Ne[i]),'{:5.2e}'.format(elem_abun_KB_error.Cl[i]),'{:09.5f}'.format(elem_ICF_KB_error.Cl[i]),'{:5.2e}'.format(elem_abun_DIMS_error.Cl[i]),'{:09.5f}'.format(elem_ICF_DIMS_error.Cl[i]),'{:5.2e}'.format(elem_abun_KB_error.Ar[i]),'{:09.5f}'.format(elem_ICF_KB_error.Ar[i]),'{:5.2e}'.format(elem_abun_DIMS_error.Ar[i]),'{:09.5f}'.format(elem_ICF_DIMS_error.Ar[i]),'{:9.2e}'.format(elem_abun_KB_ratio_error.NO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio_error.NO[i]),'{:9.2e}'.format(elem_abun_KB_ratio_error.SO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio_error.SO[i]), '{:9.2e}'.format(elem_abun_KB_ratio_error.NeO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio_error.NeO[i]),'{:9.2e}'.format(elem_abun_KB_ratio_error.ClO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio_error.ClO[i]),'{:9.2e}'.format(elem_abun_KB_ratio_error.ClO[i]),'{:10.2e}'.format(elem_abun_DIMS_ratio_error.ClO[i]), file=file6)
        
            
        os.remove("test.dat")
    file4.close() 
    file5.close() 
    file10.close()

    pn.log_.close_file()


#############################################################################################

# Te and Ne plot    
    index_line1=param_estimated.index("Te(NII6548_84)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(NII6548_84)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(NII6548_84)_Ne(SII6716_31)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        fig,ax1=plt.subplots(figsize=(12,7))
        plt.errorbar(ang,Te.NIISII,yerr=Te_error.NIISII,label="Te ([N II])")
        plt.errorbar(ang,Ne.NIISII,yerr=Ne_error.NIISII,label="Ne ([S II])")
#        plt.plot(ang,Te.NIISII,label="Te ([N II])")
#        plt.plot(ang,Ne.NIISII,label="Ne ([S II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_Ne(SII6716_31)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_Ne(SII6716_31)_both_angles.pdf")
#        plt.show()

# Te plot    
    if param_requered[index_line1]=="yes":
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(NII6548_84)_only_Ne(SII6716_31)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        fig,ax1=plt.subplots(figsize=(12,7))
        plt.errorbar(ang,Te.NIISII,yerr=Te_error.NIISII,label="Te ([N II])") 
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_only_Ne(SII6716_31)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_only_Ne(SII6716_31)_angles.pdf")

#Ne plot
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(SII6716_31)_only_Te(NII6548_84)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        fig,ax1=plt.subplots(figsize=(12,7))
        plt.errorbar(ang,Ne.NIISII,yerr=Ne_error.NIISII,label="Ne ([S II])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(NII6548_84)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(NII6548_84)_angles.pdf")


#############################################################################################
# Te and Ne plot    
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OI6300_63)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OI6300_63)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OI6300_63)_Ne(SII6716_31)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OISII,yerr=Te_error.OISII,label="Te ([O I])")
        plt.errorbar(ang,Ne.OISII,yerr=Ne_error.OISII,label="Ne ([S II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(OI6300_63)_Ne(SII6716_31)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OI6300_63)_Ne(SII6716_31)_both_angles.pdf")
#        plt.show()

#Te plot    
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OI6300_63)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OI6300_63)")        
    if param_requered[index_line1]=="yes":
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OI6300_63)_only_Ne(SII6716_31)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OISII,yerr=Te_error.OISII,label="Te ([O I])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(OI6300_63)_only_Ne(SII6716_31)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OI6300_63)_only_Ne(SII6716_31)_angles.pdf")
 
#Ne plot    
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OI6300_63)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OI6300_63)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(SII6716_31)_only_Te(OI6300_63)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.OISII,yerr=Ne_error.OISII,label="Ne ([S II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(OI6300_63)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(OI6300_63)_angles.pdf")
       
    
#############################################################################################
# Te and Ne plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OII3727_29_7320_30)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OII3727_29_7320_30)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OII3727_29_7320_30)_Ne(SII6716_31)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIISII,yerr=Te_error.OIISII,label="Te ([O II])")
        plt.errorbar(ang,Ne.OIISII,yerr=Ne_error.OIISII,label="Ne ([S II])")  
#        plt.plot(ang,Te.OIISII,label="Te ([O II])")
#        plt.plot(ang,Ne.OIISII,label="Ne ([S II])") 
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(OII3727_29_7320_30)_Ne(SII6716_31)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OII3727_29_7320_30)_Ne(SII6716_31)_both_angles.pdf")
#        plt.show() 

# Te plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OII3727_29_7320_30)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OII3727_29_7320_30)")        
    if param_requered[index_line1]=="yes":
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OII3727_29_7320_30)_only_Ne(SII6716_31)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIISII,yerr=Te_error.OIISII,label="Te ([O II])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(OII3727_29_7320_30)_only_Ne(SII6716_31)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OII3727_29_7320_30)_only_Ne(SII6716_31)_angles.pdf")

# Ne plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OII3727_29_7320_30)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OII3727_29_7320_30)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(SII6716_31)_only_Te(OII3727_29_7320_30)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.OIISII,yerr=Ne_error.OIISII,label="Ne ([S II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(OII3727_29_7320_30)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(OII3727_29_7320_30)_angles.pdf")
    

#############################################################################################
# Te and Ne plot         
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OIII4959_5007)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OIII4959_5007)_Ne(SII6716_31)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIIISII,yerr=Te_error.OIIISII,label="Te ([O III])")
        plt.errorbar(ang,Ne.OIIISII,yerr=Ne_error.OIIISII,label="Ne ([S II])")  
#        plt.plot(ang,Te.OIIISII,label="Te ([O III])")
#        plt.plot(ang,Ne.OIIISII,label="Ne ([S II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_Ne(SII6716_31)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_Ne(SII6716_31)_both_angles.pdf")
#        plt.show() 

# Te plot         
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OIII4959_5007)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OIII4959_5007)_only_Ne(SII6716_31)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIIISII,yerr=Te_error.OIIISII,label="Te ([O III])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_only_Ne(SII6716_31)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_only_Ne(SII6716_31)_angles.pdf")

# Ne plot         
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(OIII4959_5007)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(SII6716_31)_only_Te(OIII4959_5007)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.OIIISII,yerr=Ne_error.OIIISII,label="Ne ([S II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(OIII4959_5007)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(OIII4959_5007)_angles.pdf")



#############################################################################################
# Te and Ne plot       
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(SIII6312_9069)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(SIII6312_9069)_Ne(SII6716_31)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.SIIISII,yerr=Te_error.SIIISII,label="Te ([S III])")
        plt.errorbar(ang,Ne.SIIISII,yerr=Ne_error.SIIISII,label="Ne ([S II])")  
#        plt.plot(ang,Te.OISII,label="Te ([S III])")
#        plt.plot(ang,Ne.OISII,label="Ne ([S II])") 
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_Ne(SII6716_31)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_Ne(SII6716_31)_both_angles.pdf")
#        plt.show() 
    
# Te plot       
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(SIII6312_9069)")        
    if param_requered[index_line1]=="yes":
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(SIII6312_9069)_only_Ne(SII6716_31)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.SIIISII,yerr=Te_error.SIIISII,label="Te ([S III])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_only_Ne(SII6716_31)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_only_Ne(SII6716_31)_angles.pdf")
 
# Ne plot       
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(SII6716_31)")
    index_line2=param_estimated.index("Ne(SII6716_31)_Te(SIII6312_9069)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(SII6716_31)_only_Te(SIII6312_9069)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.SIIISII,yerr=Ne_error.SIIISII,label="Ne ([S II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(SIII6312_9069)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(SII6716_31)_only_Te(SIII6312_9069)_angles.pdf")
   
   
#############################################################################################
# Te and Ne plot       
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OII3727_29_7320_30)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(OII3727_29_7320_30)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OII3727_29_7320_30)_Ne(OII3727_29)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIIOII,yerr=Te_error.OIIOII,label="Te ([O II])")
        plt.errorbar(ang,Ne.OIIOII,yerr=Ne_error.OIIOII,label="Ne ([O II])")  
#        plt.plot(ang,Te.OISII,label="Te ([O II])")
#        plt.plot(ang,Ne.OISII,label="Ne ([O II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(OII3727_29_7320_30)_Ne(OII3727_29)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OII3727_29_7320_30)_Ne(OII3727_29)_both_angles.pdf")
#        plt.show() 
    
# Te plot       
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OII3727_29_7320_30)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(OII3727_29_7320_30)")        
    if param_requered[index_line1]=="yes":
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OII3727_29_7320_30)_only_Ne(OII3727_29)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIIOII,yerr=Te_error.OIIOII,label="Te ([O II])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(OII3727_29_7320_30)_only_Ne(OII3727_29)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OII3727_29_7320_30)_only_Ne(OII3727_29)_angles.pdf")

# Ne plot       
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OII3727_29_7320_30)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(OII3727_29_7320_30)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(OII3727_29)_only_Te(OII3727_29_7320_30)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.OIIOII,yerr=Ne_error.OIIOII,label="Ne ([O II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(OII3727_29)_only_Te(OII3727_29_7320_30)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(OII3727_29)_only_Te(OII3727_29_7320_30)_angles.pdf")

    
#############################################################################################
# Te and Ne plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(NII6548_84)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(NII6548_84)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(NII6548_84)_Ne(OII3727_29)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.NIIOII,yerr=Te_error.NIIOII,label="Te ([N II])")
        plt.errorbar(ang,Ne.NIIOII,yerr=Ne_error.NIIOII,label="Ne ([O II])")  
#        plt.plot(ang,Te.OISII,label="Te ([N II])")
#        plt.plot(ang,Ne.OISII,label="Ne ([O II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_Ne(OII3727_29)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_Ne(OII3727_29)_both_angles.pdf")
#        plt.show() 

# Te plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(NII6548_84)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(NII6548_84)")        
    if param_requered[index_line1]=="yes" :
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(NII6548_84)_only_Ne(OII3727_29)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.NIIOII,yerr=Te_error.NIIOII,label="Te ([N II])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_only_Ne(OII3727_29)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_only_Ne(OII3727_29)_angles.pdf")

# Ne plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(NII6548_84)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(NII6548_84)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(OII3727_29)_only_Te(NII6548_84)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.NIIOII,yerr=Ne_error.NIIOII,label="Ne ([O II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(OII3727_29)_only_Te(NII6548_84)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(OII3727_29)_only_Te(NII6548_84)_angles.pdf")

#############################################################################################
# Te and Ne plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OI6300_63)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(OI6300_63)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OI6300_63)_Ne(OII3727_29)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIOII,yerr=Te_error.OIOII,label="Te ([O I])")
        plt.errorbar(ang,Ne.OIOII,yerr=Ne_error.OIOII,label="Ne ([O II])")  
#        plt.plot(ang,Te.OISII,label="Te ([O I])")
#        plt.plot(ang,Ne.OISII,label="Ne ([O II])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(OI6300_63)_Ne(OII3727_29)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OI6300_63)_Ne(OII3727_29)_both_angles.pdf")
#        plt.show() 

# Te plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OI6300_63)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(OI6300_63)")        
    if param_requered[index_line1]=="yes" :
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OI6300_63)_only_Ne(OII3727_29)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIOII,yerr=Te_error.OIOII,label="Te ([O I])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(OI6300_63)_only_Ne(OII3727_29)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OI6300_63)_only_Ne(OII3727_29)_angles.pdf")

# Ne plot     
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OI6300_63)_Ne(OII3727_29)")
    index_line2=param_estimated.index("Ne(OII3727_29)_Te(OI6300_63)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(OII3727_29)_only_Te(OI6300_63)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.OIOII,yerr=Ne_error.OIOII,label="Ne ([O II])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(OII3727_29)_only_Te(OI6300_63)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(OII3727_29)_only_Te(OI6300_63)_angles.pdf")



#############################################################################################
# Te and Ne plot 
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(OIII4959_5007)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OIII4959_5007)_Ne(ClIII5517_38)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIIIClIII,yerr=Te_error.OIIIClIII,label="Te ([O III])")
        plt.errorbar(ang,Ne.OIIIClIII,yerr=Ne_error.OIIIClIII,label="Ne ([Cl III])")  
#        plt.plot(ang,Te.OISII,label="Te ([O III])")
#        plt.plot(ang,Ne.OISII,label="Ne ([Cl III])") 
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_Ne(ClIII5517_38)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_Ne(ClIII5517_38)_both_angles.pdf")
#        plt.show() 

#############################################################################################
# Te plot 
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(OIII4959_5007)")        
    if param_requered[index_line1]=="yes":
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OIII4959_5007)_only_Ne(ClIII5517_38)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIIIClIII,yerr=Te_error.OIIIClIII,label="Te ([O III])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_only_Ne(ClIII5517_38)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_only_Ne(ClIII5517_38)_angles.pdf")

#############################################################################################
# Ne plot 
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(OIII4959_5007)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(ClIII5517_38)_only_Te(OIII4959_5007)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.OIIIClIII,yerr=Ne_error.OIIIClIII)  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(ClIII5517_38)_only_Te(OIII4959_5007)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(ClIII5517_38)_only_Te(OIII4959_5007)_angles.pdf")

#############################################################################################
# Te and Ne plot 
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(SIII6312_9069)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(SIII6312_9069)_Ne(ClIII5517_38)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.SIIIClIII,yerr=Te_error.SIIIClIII,label="Te ([S III])")
        plt.errorbar(ang,Ne.SIIIClIII,yerr=Ne_error.SIIIClIII,label="Ne ([Cl III])")  
#        plt.plot(ang,Te.OISII,label="Te ([S III])")
#        plt.plot(ang,Ne.OISII,label="Ne ([Cl III])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_Ne(ClIII5517_38)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_Ne(ClIII5517_38)_both_angles.pdf")
#        plt.show() 

# Te plot 
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(SIII6312_9069)")        
    if param_requered[index_line1]=="yes" :
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(SIII6312_9069)_only_Ne(ClIII5517_38)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.SIIIClIII,yerr=Te_error.SIIIClIII,label="Te ([S III])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_only_Ne(ClIII5517_38)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_only_Ne(ClIII5517_38)_angles.pdf")

# Ne plot 
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(SIII6312_9069)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(ClIII5517_38)_only_Te(SIII6312_9069)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.SIIIClIII,yerr=Ne_error.SIIIClIII,label="Ne ([Cl III])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(ClIII5517_38)_only_Te(SIII6312_9069)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(ClIII5517_38)_only_Te(SIII6312_9069)_angles.pdf")


#############################################################################################
# Te and Ne plot
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(ArVI4712_40)")
    index_line2=param_estimated.index("Ne(ArVI4712_40)_Te(OIII4959_5007)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OIII4959_5007)_Ne(ArVI4712_40)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIIIArVI,yerr=Te_error.OIIIArVI,label="Te ([O III])")
        plt.errorbar(ang,Ne.OIIIArVI,yerr=Ne_error.OIIIArVI,label="Ne ([Ar VI])")  
#        plt.plot(ang,Te.OISII,label="Te ([O III])")
#        plt.plot(ang,Ne.OISII,label="Ne ([Ar VI])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_Ne(ArVI4712_40)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_Ne(ArVI4712_40)_both_angles.pdf")
#        plt.show() 

# Te  plot
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(ArVI4712_40)")
    index_line2=param_estimated.index("Ne(ArVI4712_40)_Te(OIII4959_5007)")        
    if param_requered[index_line1]=="yes" :
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(OIII4959_5007)_only_Ne(ArVI4712_40)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.OIIIArVI,yerr=Te_error.OIIIArVI,label="Te ([O III])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_only_Ne(ArVI4712_40)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(OIII4959_5007)_only_Ne(ArVI4712_40)_angles.pdf")

# Ne plot
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(OIII4959_5007)_Ne(ArVI4712_40)")
    index_line2=param_estimated.index("Ne(ArVI4712_40)_Te(OIII4959_5007)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(ArVI4712_40)_only_Te(OIII4959_5007)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.OIIIArVI,yerr=Ne_error.OIIIArVI,label="Ne ([Ar VI])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(ArVI4712_40)_only_Te(OIII4959_5007)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(ArVI4712_40)_only_Te(OIII4959_5007)_angles.pdf")

#############################################################################################
# Te and Ne plot    
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(ArVI4712_40)")
    index_line2=param_estimated.index("Ne(ArVI4712_40)_Te(SIII6312_9069)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(SIII6312_9069)_Ne(ArVI4712_40)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.SIIIArVI,yerr=Te_error.SIIIArVI,label="Te ([S III])")
        plt.errorbar(ang,Ne.SIIIArVI,yerr=Ne_error.SIIIArVI,label="Ne ([Ar VI])")  
#        plt.plot(ang,Te.OISII,label="Te ([S III])")
#        plt.plot(ang,Ne.OISII,label="Ne ([Ar VI])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_Ne(ArVI4712_40)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_Ne(ArVI4712_40)_both_angles.pdf")
#        plt.show()

# Te plot    
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(ArVI4712_40)")
    index_line2=param_estimated.index("Ne(ArVI4712_40)_Te(SIII6312_9069)")        
    if param_requered[index_line1]=="yes" :
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(SIII6312_9069)_only_Ne(ArVI4712_40)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.SIIIArVI,yerr=Te_error.SIIIArVI,label="Te ([S III])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_only_Ne(ArVI4712_40)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(SIII6312_9069)_only_Ne(ArVI4712_40)_angles.pdf")

# Ne plot    
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(SIII6312_9069)_Ne(ArVI4712_40)")
    index_line2=param_estimated.index("Ne(ArVI4712_40)_Te(SIII6312_9069)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(ArVI4712_40)_only_Te(SIII6312_9069)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.SIIIArVI,yerr=Ne_error.SIIIArVI,label="Ne ([Ar VI])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Ne(ArVI4712_40)_only_Te(SIII6312_9069)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(ArVI4712_40)_only_Te(SIII6312_9069)_angles.pdf")
 
 
 
#############################################################################################
# Te and Ne plot
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(NII6548_84)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(NII6548_84)")        
    if param_requered[index_line1]=="yes" or param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(NII6548_84)_Ne(ClIII5517_38)_both")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.NIIClIII,yerr=Te_error.NIIClIII,label="Te ([N II])")
        plt.errorbar(ang,Ne.NIIClIII,yerr=Ne_error.NIIClIII,label="Ne ([Cl III])")  
#        plt.plot(ang,Te.OISII,label="Te ([S III])")
#        plt.plot(ang,Ne.OISII,label="Ne ([Ar VI])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K)/Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_Ne(ClIII5517_38)_both_angles.png")
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_Ne(ClIII5517_38)_both _angles.pdf")
#        plt.show()

# Te  plot
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(NII6548_84)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(NII6548_84)")        
    if param_requered[index_line1]=="yes":
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Te(NII6548_84)_only_Ne(ClIII5517_38)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Te.NIIClIII,yerr=Te_error.NIIClIII,label="Te ([N II])")
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Te (K) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_only_Ne(ClIII5517_38)_angles.png")
        fig.savefig("output_angles_plots/fig_Te(NII6548_84)_only_Ne(ClIII5517_38)_angles.pdf")

# Ne plot
    fig,ax1=plt.subplots(figsize=(12,7))
    index_line1=param_estimated.index("Te(NII6548_84)_Ne(ClIII5517_38)")
    index_line2=param_estimated.index("Ne(ClIII5517_38)_Te(NII6548_84)")        
    if param_requered[index_line2]=='yes':
        
        # Read the necessary parameters (x or y axis) for the plots
        index_plotname=par_plotname.index("Ne(ClIII5517_38)_only_Te(NII6548_84)")
        ymin=float(par_plotymin[index_plotname])
        ymax=float(par_plotymax[index_plotname])
        
        plt.errorbar(ang,Ne.NIIClIII,yerr=Ne_error.NIIClIII,label="Ne ([Cl III])")  
        plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
        plt.ylabel("Ne (cm$^{-3}$) ",fontsize=21,labelpad=6)
        ax1.tick_params(axis='x', labelsize=21, pad=6)
        ax1.tick_params(axis='y', labelsize=21, pad=6)
        ax1.yaxis.get_offset_text().set_fontsize(18)
        plt.ylim(ymin,ymax)
        plt.legend(loc='upper right',ncol=2,prop={'size':18})
        fig.savefig("output_angles_plots/fig_Ne(ClIII5517_38)_only_Te(NII6548_84)_angles.png")
        fig.savefig("output_angles_plots/fig_Ne(ClIII5517_38)_only_Te(NII6548_84)_angles.pdf")


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

    
    
    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB.He,yerr=elem_abun_KB_error.He,label="He (KB)")
    plt.errorbar(ang,elem_abun_DIMS.He,yerr=elem_abun_DIMS_error.He,label="He (DIMS)")
    index_line1=line_names.index("HeI_5876s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.HeIa,yerr=ion_abun_error.HeIa,label="He+ 5876")
    index_line1=line_names.index("HeI_6678s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.HeIb,yerr=ion_abun_error.HeIb,label="He+ 6678")
    index_line1=line_names.index("HeII_4686s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.HeIIa,yerr=ion_abun_error.HeIIa,label="He++ 4686")
    index_line1=line_names.index("HeII_5412s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.HeIIb,yerr=ion_abun_error.HeIIb,label="He++ 5412")
    
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("He_abundances")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ionic and Total abudnance (X/H+) ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/He_abundances_angles.png")
    fig.savefig("output_angles_plots/He_abundances_angles.pdf")
#    plt.show()
    
#############################################################################################
    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB.O,yerr=elem_abun_KB_error.O,label="O (KB)")
    plt.errorbar(ang,elem_abun_DIMS.O,yerr=elem_abun_DIMS_error.O,label="O (DIMS)")
    index_line1=line_names.index("O3_4363s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.OIIIa,yerr=ion_abun_error.OIIIa,label="O++ 4363")
    index_line1=line_names.index("O3_5007s")
    index_line2=line_names.index("O3_4959s")

    if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes": 
        plt.errorbar(ang,ion_abun.OIIIc,yerr=ion_abun_error.OIIIc,label="O++ 5007")
    if lines_available[index_line1]=="yes" and lines_available[index_line2]=="no": 
        plt.errorbar(ang,ion_abun.OIIIc,yerr=ion_abun_error.OIIIc,label="O++ 5007")
    if lines_available[index_line1]=="no" and lines_available[index_line2]=="yes": 
        plt.errorbar(ang,ion_abun.OIIIb,yerr=ion_abun_error.OIIIb,label="O++ 4959")
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("O_abundances_OIII")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ionic and Total abudnance (X/H+) ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/O_abundances_angles_a.png")
    fig.savefig("output_angles_plots/O_abundances_angles_a.pdf")
#    plt.show()
    

#############################################################################################    
    fig,ax2=plt.subplots(figsize=(12,7))
#    plt.errorbar(ang,elem_abun_KB.O,yerr=elem_abun_KB_error.O,label="O (KB)")
#    plt.errorbar(ang,elem_abun_DIMS.O,yerr=elem_abun_DIMS_error.O,label="O (DIMS)")
    index_line1=line_names.index("O2_3727s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.OIIa,yerr=ion_abun_error.OIIa,label="O+ 3727")
    index_line1=line_names.index("O2_3729s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.OIIb,yerr=ion_abun_error.OIIb,label="O+ 3729")
    index_line1=line_names.index("O2_7320s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.OIIc,yerr=ion_abun_error.OIIc,label="O+ 7320")
    index_line1=line_names.index("O2_7330s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.OIId,yerr=ion_abun_error.OIId,label="O+ 7330")
    index_line1=line_names.index("O1_5577s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.OIa,yerr=ion_abun_error.OIa,label="Oo 5755")    
    index_line1=line_names.index("O1_6300s")
    index_line2=line_names.index("O1_6363s")
    if lines_available[index_line1]=="yes" and lines_available[index_line2]=="yes": 
        plt.errorbar(ang,ion_abun.OIb,yerr=ion_abun_error.OIb,label="Oo 6300")
    if lines_available[index_line1]=="yes" and lines_available[index_line2]=="no": 
        plt.errorbar(ang,ion_abun.OIb,yerr=ion_abun_error.OIb,label="Oo 6300")
    if lines_available[index_line1]=="no" and lines_available[index_line2]=="yes": 
        plt.errorbar(ang,ion_abun.OIc,yerr=ion_abun_error.OIc,label="Oo 6363")

    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("O_abundances_OIIOI")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])

    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ionic and Total abudnance (X/H+) ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/O_abundances_angles_b.png")
    fig.savefig("output_angles_plots/O_abundances_angles_b.pdf")
#    plt.show()
    


#############################################################################################    
    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB.N,yerr=elem_abun_KB_error.N,label="N (KB)")
    plt.errorbar(ang,elem_abun_DIMS.N,yerr=elem_abun_DIMS_error.N,label="N (DIMS)")
    index_line1=line_names.index("N1_5199s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.NI,yerr=ion_abun_error.NI,label="No 5199")
    index_line1=line_names.index("N2_5755s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.NIIa,yerr=ion_abun_error.NIIa,label="N+ 5755")
    index_line1=line_names.index("N2_6548s")
    index_line2=line_names.index("N2_6583s")
    if lines_available[index_line1]=="yes"and lines_available[index_line2]=="yes": 
        plt.errorbar(ang,ion_abun.NIIc,yerr=ion_abun_error.NIIc,label="N+ 6584")
    if lines_available[index_line1]=="no"and lines_available[index_line2]=="yes": 
        plt.errorbar(ang,ion_abun.NIIc,yerr=ion_abun_error.NIIc,label="N+ 6584")
    if lines_available[index_line1]=="yes"and lines_available[index_line2]=="no":
        plt.errorbar(ang,ion_abun.NIIb,yerr=ion_abun_error.NIIb,label="N+ 6548")
    

    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("N_abundances")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ionic and Total abudnance (X/H+) ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3, prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/N_abundances_angles.png")
    fig.savefig("output_angles_plots/N_abundances_angles.pdf")
    

#############################################################################################        
    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB.S,yerr=elem_abun_KB_error.S,label="S (KB)")
    plt.errorbar(ang,elem_abun_DIMS.S,yerr=elem_abun_DIMS_error.S,label="S (DIMS)")
    index_line1=line_names.index("S2_6716s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.SIIa,yerr=ion_abun_error.SIIa,label="S+ 6716")
    index_line1=line_names.index("S2_6731s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.SIIb,yerr=ion_abun_error.SIIb,label="S+ 6731")
    index_line1=line_names.index("S3_6312s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.SIIIa,yerr=ion_abun_error.SIIIa,label="S++ 6312")
    index_line1=line_names.index("S3_9069s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.SIIIb,yerr=ion_abun_error.SIIIb,label="S++ 9069")
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("S_abundances")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ionic and Total abudnance (X/H+) ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.ylim(ymin,ymax)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    fig.savefig("output_angles_plots/S_abundances_angles.png")
    fig.savefig("output_angles_plots/S_abundances_angles.pdf")
#    plt.show()
    
 
#############################################################################################    
    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB.Ne,yerr=elem_abun_KB_error.Ne,label="Ne (KB)")
    plt.errorbar(ang,elem_abun_DIMS.Ne,yerr=elem_abun_DIMS_error.Ne,label="Ne (DIMS)")
    index_line1=line_names.index("Ne3_3868s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.NeIIIa,yerr=ion_abun_error.NeIIIa,label="Ne++ 3868")
    index_line1=line_names.index("Ne3_3967s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.NeIIIb,yerr=ion_abun_error.NeIIIb,label="Ne++ 3967")

    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("Ne_abundances")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])

    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ionic and Total abudnance (X/H+) ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=2,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/Ne_abundances_angles.png")
    fig.savefig("output_angles_plots/Ne_abundances_angles.pdf")
    

#############################################################################################        
    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB.Cl,yerr=elem_abun_KB_error.Cl,label="Cl (KB)")
    plt.errorbar(ang,elem_abun_DIMS.Cl,yerr=elem_abun_DIMS_error.Cl,label="Cl (DIMS)")
    index_line1=line_names.index("Cl3_5517s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.ClIIIa,yerr=ion_abun_error.ClIIIa,label="Cl++ 5517 ")
    index_line1=line_names.index("Cl3_5538s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.ClIIIb,yerr=ion_abun_error.ClIIIb,label="Cl++ 5538")
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("Cl_abundances")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ionic and Total abudnance (X/H+) ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=2,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/Cl_abundances_angles.png")
    fig.savefig("output_angles_plots/Cl_abundances_angles.pdf")
#    plt.show()

#############################################################################################    
    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB.Ar,yerr=elem_abun_KB_error.Ar,label="Ar (KB)")
    plt.errorbar(ang,elem_abun_DIMS.Ar,yerr=elem_abun_DIMS_error.Ar,label="Ar (DIMS)")
    index_line1=line_names.index("Ar3_7136s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.ArIII,yerr=ion_abun_error.ArIII,label="Ar++ 7136")
    index_line1=line_names.index("Ar4_4712s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.ArIVa,yerr=ion_abun_error.ArIVa,label="Ar+++ 4712")
    index_line1=line_names.index("Ar4_4740s")
    if lines_available[index_line1]=="yes": plt.errorbar(ang,ion_abun.ArIVb,yerr=ion_abun_error.ArIVb,label="Ar+++ 4740")
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("Ar_abundances")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ionic and Total abudnance (X/H+) ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(-3,4),axis='both')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=2,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/Ar_abundances_angles.png")
    fig.savefig("output_angles_plots/Ar_abundances_angles.pdf")
#    plt.show()



#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################


    fig,ax2=plt.subplots(figsize=(12,7))

    plt.errorbar(ang,elem_abun_KB_ratio.NO,yerr=elem_abun_KB_ratio_error.NO,label="N/O (KB)")
    plt.errorbar(ang,elem_abun_DIMS_ratio.NO,yerr=elem_abun_DIMS_ratio_error.NO,label="N/O (DIMS)")

    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("NO_abundances_ratio")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("N/O abundance ratio ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/NO_abundances_ratio_angles.png")
    fig.savefig("output_angles_plots/NO_abundances_ratio_angles.pdf")
    
#############################################################################################

    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB_ratio.SO,yerr=elem_abun_KB_ratio_error.SO,label="S/O (KB)")
    plt.errorbar(ang,elem_abun_DIMS_ratio.SO,yerr=elem_abun_DIMS_ratio_error.SO,label="S/O (DIMS)")

    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("SO_abundances_ratio")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("S/O abundance ratio ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/SO_abundances_ratio_angles.png")
    fig.savefig("output_angles_plots/SO_abundances_ratio_angles.pdf")

#############################################################################################

    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB_ratio.NeO,yerr=elem_abun_KB_ratio_error.NeO,label="Ne/O (KB)")
    plt.errorbar(ang,elem_abun_DIMS_ratio.NeO,yerr=elem_abun_DIMS_ratio_error.NeO,label="Ne/O (DIMS)")

    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("NeO_abundances_ratio")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ne/O abundance ratio ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/NeO_abundances_ratio_angles.png")
    fig.savefig("output_angles_plots/NeO_abundances_ratio_angles.pdf")
#############################################################################################

    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB_ratio.ClO,yerr=elem_abun_KB_ratio_error.ClO,label="Cl/O (KB)")
    plt.errorbar(ang,elem_abun_DIMS_ratio.ClO,yerr=elem_abun_DIMS_ratio_error.ClO,label="Cl/O (DIMS)")

    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("ClO_abundances_ratio")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Cl/O abundance ratio ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/ClO_abundances_ratio_angles.png")
    fig.savefig("output_angles_plots/ClO_abundances_ratio_angles.pdf")

#############################################################################################

    fig,ax2=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_abun_KB_ratio.ArO,yerr=elem_abun_KB_ratio_error.ArO,label="Ar/O (KB)")
    plt.errorbar(ang,elem_abun_DIMS_ratio.ArO,yerr=elem_abun_DIMS_ratio_error.ArO,label="Ar/O (DIMS)")
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("ArO_abundances_ratio")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("Ar/O abundance ratio ",fontsize=21,labelpad=6)
    ax2.ticklabel_format(style='sci',scilimits=(0,0),axis='y')
    ax2.tick_params(axis='x', labelsize=21, pad=6)
    ax2.tick_params(axis='y', labelsize=21, pad=6)
    ax2.yaxis.get_offset_text().set_fontsize(18)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/ArO_abundances_ratio_angles.png")
    fig.savefig("output_angles_plots/ArO_abundances_ratio_angles.pdf")


#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################


    fig,ax1=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_ICF_KB.N,yerr=elem_ICF_KB_error.N,label="N_ICF (KB)")
    plt.errorbar(ang,elem_ICF_DIMS.N,yerr=elem_ICF_DIMS_error.N,label="N_ICF (DIMS)")
#    plt.errorbar(ang,elem_ICF_KB.S,yerr=elem_ICF_KB_error.S,label="S_ICF (KB)")
#    plt.errorbar(ang,elem_ICF_DIMS.S,yerr=elem_ICF_DIMS_error.S,label="S_ICF (DIMS)")

    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("N_ICF")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])


    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("ICFs ",fontsize=21,labelpad=6)
    ax1.tick_params(axis='x', labelsize=21, pad=6)
    ax1.tick_params(axis='y', labelsize=21, pad=6)
    plt.legend(loc='upper right',ncol=2,prop={'size':18})
    plt.ylim(ymin,ymax)
    fig.savefig("output_angles_plots/N_ICFs_angle.png")
    fig.savefig("output_angles_plots/N_ICFs_angle.pdf")

#############################################################################################        
    fig,ax1=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_ICF_KB.O,yerr=elem_ICF_KB_error.O,label="O_ICF (KB)")
    plt.errorbar(ang,elem_ICF_DIMS.O,yerr=elem_ICF_DIMS_error.O,label="O_ICF (DIMS)")
    plt.errorbar(ang,elem_ICF_KB.S,yerr=elem_ICF_KB_error.S,label="S_ICF (KB)")
    plt.errorbar(ang,elem_ICF_DIMS.S,yerr=elem_ICF_DIMS_error.S,label="S_ICF (DIMS)")
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("O_S_ICF")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("ICFs ",fontsize=21,labelpad=6)
    ax1.tick_params(axis='x', labelsize=21, pad=6)
    ax1.tick_params(axis='y', labelsize=21, pad=6)
    plt.ylim(ymin,ymax)
    plt.legend(loc='upper right',ncol=2,prop={'size':18})
    fig.savefig("output_angles_plots/O_S_ICFs_angle.png")
    fig.savefig("output_angles_plots/O_S_ICFs_angle.pdf")
#       plt.show()

#############################################################################################    
    fig,ax1=plt.subplots(figsize=(12,7))
    plt.errorbar(ang,elem_ICF_KB.Cl,yerr=elem_ICF_KB_error.Cl,label="Cl_ICF (KB)")
    plt.errorbar(ang,elem_ICF_DIMS.Cl,yerr=elem_ICF_DIMS_error.Cl,label="Cl_ICF (DIMS)")
    plt.errorbar(ang,elem_ICF_KB.Ne,yerr=elem_ICF_KB_error.Ne,label="Ne_ICF (KB)")
    plt.errorbar(ang,elem_ICF_DIMS.Ne,yerr=elem_ICF_DIMS_error.Ne,label="Ne_ICF (DIMS)")
    plt.errorbar(ang,elem_ICF_KB.Ar,yerr=elem_ICF_KB_error.Ar,label="Ar_ICF (KB)")
    plt.errorbar(ang,elem_ICF_DIMS.Ar,yerr=elem_ICF_DIMS_error.Ar,label="Ar_ICF (DIMS)")
    
    # Read the necessary parameters (x or y axis) for the plots
    index_plotname=par_plotname.index("Cl_Ne_Ar_ICF")
    ymin=float(par_plotymin[index_plotname])
    ymax=float(par_plotymax[index_plotname])
    
    plt.xlabel("Position angle of simulated slits (degrees)",fontsize=21,labelpad=6)
    plt.ylabel("ICFs ",fontsize=21,labelpad=6)
    ax1.tick_params(axis='x', labelsize=21, pad=6)
    ax1.tick_params(axis='y', labelsize=21, pad=6)
    plt.ylim(ymin,ymax)
    plt.legend(loc='upper right',ncol=3,prop={'size':18})
    fig.savefig("output_angles_plots/Ne_Cl_Ar_ICFs_angle.png")
    fig.savefig("output_angles_plots/Ne_Cl_Ar_ICFs_angle.pdf")



    return Te,Ne     

        
