# ionicabundances_2D_script.py:
# Calculates the ionic abundances for each spaxel in 2D maps.
# Te and Ne used for the ionic abundances computations are provided by the user.
# If both parameters have zero values, then the routine computes the averabe value 
# from the rotation analysis module. 
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
import matplotlib.patches as patches
import matplotlib.font_manager as font_manager
import random
import math
import numpy as np
import seaborn as sns
import pyneb as pn
import os, fnmatch
import warnings
from numpy import inf


from satellite import diagnotic_diagrams_script as dds
from satellite import mean_value_script as mvs


def ionicabundances(obs,obs_fake,Te,Te_error,Te_fake,Ne,Ne_error,Ne_fake,i,j,line_names,lines_available,ion_abun,ion_abun_error,param_estimated,param_requered,Te_PA,Ne_PA,param_mod_name):

    
    file10 = open('general_output_file.txt','a')
    warnings.filterwarnings("ignore")
    
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
    
##############################################################
    He1 = pn.RecAtom('He',1)
    
    index_name=param_estimated.index("He+_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

    
    index1=index2=0
    index_line=line_names.index("HeI_5876s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        He1a_abund_fake = He1.getIonAbundance(int_ratio=obs_fake.getIntens(0)['He1r_5876A'], tem=Teinput_fake, den=Neinput_fake, to_eval='L(5876)', Hbeta=100.0)
        He1a_abund = He1.getIonAbundance(int_ratio=obs.getIntens(0)['He1r_5876A'], tem=Teinput, den=Neinput, to_eval='L(5876)', Hbeta=100.0)
        index1=1
        
        He1a_abund_fake1=np.nan_to_num(He1a_abund_fake)
        He1a_abund_fake1[He1a_abund_fake1 == inf] =  0

        He1a_abund_fake[He1a_abund_fake == inf] =  0
        He1a_abund_fake=np.nan_to_num(He1a_abund_fake)
        He1a_abund_fake = [i2 for i2 in He1a_abund_fake if i2 != 0]
        He1a_abund_fake = [i2 for i2 in He1a_abund_fake if i2 < He1a_abund*3 and i2>He1a_abund*0.3]
        He1a_abund_fake=np.std(np.nan_to_num(He1a_abund_fake), ddof=1)
    else:
        He1a_abund_fake=0.0
        He1a_abund=[0.0,0.0]
    
    index_line=line_names.index("HeI_6678s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        He1b_abund_fake = He1.getIonAbundance(int_ratio=obs_fake.getIntens(0)['He1r_6678A'], tem=Teinput_fake, den=Neinput_fake, to_eval='L(6678)', Hbeta=100.0)
        He1b_abund = He1.getIonAbundance(int_ratio=obs.getIntens(0)['He1r_6678A'], tem=Teinput, den=Neinput, to_eval='L(6678)', Hbeta=100.0)
        index2=1
        
        He1b_abund_fake2=np.nan_to_num(He1b_abund_fake)
        He1b_abund_fake2[He1b_abund_fake2 == inf] =  0

        He1b_abund_fake[He1b_abund_fake == inf] =  0
        He1b_abund_fake=np.nan_to_num(He1b_abund_fake)
        He1b_abund_fake = [i2 for i2 in He1b_abund_fake if i2 != 0]
        He1b_abund_fake = [i2 for i2 in He1b_abund_fake if i2 < He1b_abund*3 and i2>He1b_abund*0.3]
        He1b_abund_fake=np.std(np.nan_to_num(He1b_abund_fake), ddof=1)
    else:
        He1b_abund_fake=0.0
        He1b_abund=[0.0,0.0]
    
    #He1a_abund_fake1=np.nan_to_num(He1a_abund_fake)
    #He1b_abund_fake2=np.nan_to_num(He1b_abund_fake)
    
    #He1a_abund_fake=np.std(np.nan_to_num(He1a_abund_fake), ddof=1)
    #He1b_abund_fake=np.std(np.nan_to_num(He1b_abund_fake), ddof=1)

    He1_abund_fake=0
    He1_abund_fake2 =0
    He1_abund=0
    
    if (index1==1 and index2==1):
        He1_abund_fake=np.sqrt(He1a_abund_fake**2+He1b_abund_fake**2)/2.
        He1_abund_fake2=(He1a_abund_fake1+He1b_abund_fake2)/2.
        He1_abund=(He1a_abund+He1b_abund)/2.
    elif (index1==1 and index2==0):
        He1_abund_fake=(He1a_abund_fake)
        He1_abund_fake2=(He1a_abund_fake1)
        He1_abund=(He1a_abund)
    elif (index1==0 and index2==1):
        He1_abund_fake=(He1b_abund_fake)
        He1_abund_fake2=(He1b_abund_fake2)
        He1_abund=(He1b_abund)


    
    if (He1a_abund[0]<1 and He1b_abund[0]<1):
        if(He1a_abund[0]>=0  and He1b_abund[0]>=0 ):
            ion_abun.HeIa[i,j]=He1a_abund[0]
            ion_abun.HeIb[i,j]=He1b_abund[0]
            if np.isnan(He1a_abund_fake)==True: He1a_abund_fake=0
            ion_abun_error.HeIa[i]=He1a_abund_fake
            if np.isnan(He1b_abund_fake)==True: He1b_abund_fake=0
            ion_abun_error.HeIb[i]=He1b_abund_fake
    
#    print >>file10, 'He+(5786)/H+ = {:5.2e}'.format(ion_abun.HeIa[i,j]),'He+(6678)/H+ = {:5.2e}'.format(ion_abun.HeIb[i,j])
#    print >>file10, 'error He+(5786)/H+ = {:5.2e}'.format(ion_abun_error.HeIa[i,j]),'He+(6678)/H+ = {:5.2e}'.format(ion_abun_error.HeIb[i,j])  
#############################################################

#############################################################
    He2 = pn.RecAtom('He',2)
    
    index_name=param_estimated.index("He++_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII
        
    index1=index2=0
    index_line=line_names.index("HeII_4686s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        He2a_abund_fake = He2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['He2r_4686A'], tem=Teinput_fake, den=Neinput_fake, to_eval='L(4686)', Hbeta=100.0)
        He2a_abund = He2.getIonAbundance(int_ratio=obs.getIntens(0)['He2r_4686A'], tem=Teinput, den=Neinput, to_eval='L(4686)', Hbeta=100.0)
        index1=1
        
        He2a_abund_fake1=np.nan_to_num(He2a_abund_fake)
        He2a_abund_fake1[He2a_abund_fake1 == inf] =  0
    
        He2a_abund_fake[He2a_abund_fake == inf] =  0
        He2a_abund_fake=np.nan_to_num(He2a_abund_fake)
        He2a_abund_fake = [i2 for i2 in He2a_abund_fake if i2 != 0]
        He2a_abund_fake = [i2 for i2 in He2a_abund_fake if i2 < He2a_abund*3 and i2>He2a_abund*0.3]
        He2a_abund_fake=np.std(np.nan_to_num(He2a_abund_fake), ddof=1)
    else:
        He2a_abund_fake=0.0
        He2a_abund=[0.0,0.0]
    
    index_line=line_names.index("HeII_5412s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        He2b_abund_fake = He2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['He2r_5411A'], tem=Teinput_fake, den=Neinput_fake, to_eval='L(5411)', Hbeta=100.0)
        He2b_abund = He2.getIonAbundance(int_ratio=obs.getIntens(0)['He2r_5411A'], tem=Teinput, den=Neinput, to_eval='L(5411)', Hbeta=100.0)
        index2=1
        
        He2b_abund_fake2=np.nan_to_num(He2b_abund_fake)
        He2b_abund_fake2[He2b_abund_fake2 == inf] =  0
    
        He2b_abund_fake[He2b_abund_fake == inf] =  0
        He2b_abund_fake=np.nan_to_num(He2b_abund_fake)
        He2b_abund_fake = [i2 for i2 in He2b_abund_fake if i2 != 0]
        He2b_abund_fake = [i2 for i2 in He2b_abund_fake if i2 < He2b_abund*3 and i2>He2b_abund*0.3]
        He2b_abund_fake=np.std(np.nan_to_num(He2b_abund_fake), ddof=1)
    else:
        He2b_abund_fake=0.0
        He2b_abund=[0.0,0.0]
    
    #He2a_abund_fake1=np.nan_to_num(He2a_abund_fake)
    #He2b_abund_fake2=np.nan_to_num(He2b_abund_fake)
        
    #He2a_abund_fake=np.std(np.nan_to_num(He2a_abund_fake), ddof=1)
    #He2b_abund_fake=np.std(np.nan_to_num(He2b_abund_fake), ddof=1)
    
    He2_abund_fake=0
    He2_abund_fake2=0
    He2_abund=0
    
    if (index1==1 and index2==1):
        He2_abund_fake=np.sqrt(He2a_abund_fake**2+He2b_abund_fake**2)/2.
        He2_abund_fake2=(He2a_abund_fake1+He2b_abund_fake2)/2.
        He2_abund=(He2a_abund+He2b_abund)/2.
    elif (index1==1 and index2==0):
        He2_abund_fake=(He2a_abund_fake)
        He2_abund_fake2=(He2a_abund_fake1)
        He2_abund=(He2a_abund)
    elif (index1==0 and index2==1):
        He2_abund_fake=(He2b_abund_fake)
        He2_abund_fake2=(He2b_abund_fake2)
        He2_abund=(He2b_abund)
        

    if (He2a_abund[0]<1 and He2b_abund[0]<1):
        if(He2a_abund[0]>=0  and He2b_abund[0]>=0 ):
            ion_abun.HeIIa[i,j]=He2a_abund[0]
            ion_abun.HeIIb[i,j]=He2b_abund[0]
            if np.isnan(He2a_abund_fake)==True: He2a_abund_fake=0
            ion_abun_error.HeIIa[i]=He2a_abund_fake
            if np.isnan(He2b_abund_fake)==True: He2b_abund_fake=0
            ion_abun_error.HeIIb[i]=He2b_abund_fake
    
#    print >>file10, 'He+2(4686)/H+ = {:5.2e}'.format(ion_abun.HeIIa[i,j]),'He+2(5412)/H+ = {:5.2e}'.format(ion_abun.HeIIb[i,j])
#    print >>file10, 'He+2(4686)/H+ = {:5.2e}'.format(ion_abun_error.HeIIa[i,j]),'He+2(5412)/H+ = {:5.2e}'.format(ion_abun_error.HeIIb[i,j])
#############################################################


#############################################################
    O1 = pn.Atom('O',1)

    index_name=param_estimated.index("Oo_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

    
    index1=index2=index3=0    
    index_line=line_names.index("O1_5577s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O1a_abund_fake = O1.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O1_5577A'], tem=Teinput_fake, den=Neinput_fake, to_eval='L(5577)', Hbeta=100.0)
        O1a_abund = O1.getIonAbundance(int_ratio=obs.getIntens(0)['O1_5577A'], tem=Teinput, den=Neinput, to_eval='L(5577)', Hbeta=100.0)
        index1=1
        
        O1a_abund_fake1=np.nan_to_num(O1a_abund_fake)
        O1a_abund_fake1[O1a_abund_fake1 == inf] =  0

        O1a_abund_fake[O1a_abund_fake == inf] =  0
        O1a_abund_fake=np.nan_to_num(O1a_abund_fake)
        O1a_abund_fake = [i2 for i2 in O1a_abund_fake if i2 != 0]
        O1a_abund_fake = [i2 for i2 in O1a_abund_fake if i2 < O1a_abund*3 and i2>O1a_abund*0.3]
        O1a_abund_fake=np.std(np.nan_to_num(O1a_abund_fake), ddof=1)    
    else:
        O1a_abund_fake=0.0
        O1a_abund=[0.0,0.0]
        
    index1=index2=0    
    index_line=line_names.index("O1_6300s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O1b_abund_fake = O1.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O1_6300A'], tem=Teinput_fake, den=Neinput_fake, to_eval='L(6300)', Hbeta=100.0)
        O1b_abund = O1.getIonAbundance(int_ratio=obs.getIntens(0)['O1_6300A'], tem=Teinput, den=Neinput, to_eval='L(6300)', Hbeta=100.0)
        index2=1
        
        O1b_abund_fake2=np.nan_to_num(O1b_abund_fake)
        O1b_abund_fake2[O1b_abund_fake2 == inf] =  0
    
        O1b_abund_fake[O1b_abund_fake == inf] =  0
        O1b_abund_fake=np.nan_to_num(O1b_abund_fake)
        O1b_abund_fake = [i2 for i2 in O1b_abund_fake if i2 != 0]
        O1b_abund_fake = [i2 for i2 in O1b_abund_fake if i2 < O1b_abund*3 and i2>O1b_abund*0.3]
        O1b_abund_fake=np.std(np.nan_to_num(O1b_abund_fake), ddof=1)
    else:
        O1b_abund_fake=0.0
        O1b_abund=[0.0,0.0]

    index_line=line_names.index("O1_6363s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O1c_abund_fake = O1.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O1_6364A'], tem=Teinput_fake, den=Neinput_fake, to_eval='L(6364)', Hbeta=100.0)
        O1c_abund = O1.getIonAbundance(int_ratio=obs.getIntens(0)['O1_6364A'], tem=Teinput, den=Neinput, to_eval='L(6364)', Hbeta=100.0)
        index3=1
        
        O1c_abund_fake3=np.nan_to_num(O1c_abund_fake)
        O1c_abund_fake3[O1c_abund_fake3 == inf] =  0
    
        O1c_abund_fake[O1c_abund_fake == inf] =  0
        O1c_abund_fake=np.nan_to_num(O1c_abund_fake)
        O1c_abund_fake = [i2 for i2 in O1c_abund_fake if i2 != 0]
        O1c_abund_fake = [i2 for i2 in O1c_abund_fake if i2 < O1c_abund*3 and i2>O1c_abund*0.3]
        O1c_abund_fake=np.std(np.nan_to_num(O1c_abund_fake), ddof=1)
    else:
        O1c_abund_fake=0.0
        O1c_abund=[0.0,0.0]


    #O1a_abund_fake1=np.nan_to_num(O1a_abund_fake)
    #O1b_abund_fake2=np.nan_to_num(O1b_abund_fake)
    #O1c_abund_fake3=np.nan_to_num(O1c_abund_fake)
    
    #O1a_abund_fake=np.std(np.nan_to_num(O1a_abund_fake), ddof=1)
    #O1b_abund_fake=np.std(np.nan_to_num(O1b_abund_fake), ddof=1)
    #O1c_abund_fake=np.std(np.nan_to_num(O1c_abund_fake), ddof=1)
    
    O1_abund_fake=0
    O1_abund_fake2=0
    O1_abund=0

    if (index1==1 and index2==0 and index3==0):
        O1_abund_fake=O1a_abund_fake
        O1_abund_fake2=O1a_abund_fake1
        O1_abund=O1a_abund
    elif (index1==0 and index2==1 and index3==0):
        O1_abund_fake=O1b_abund_fake
        O1_abund_fake2=O1b_abund_fake2
        O1_abund=O1b_abund
    elif (index1==0 and index2==1 and index3==1):
        O1_abund_fake=np.sqrt(O1b_abund_fake**2+O1c_abund_fake**2)/2.
        O1_abund_fake2=(O1b_abund_fake2+O1c_abund_fake3)/2.
        O1_abund=(O1b_abund+O1c_abund)/2.
    elif (index1==1 and index2==1 and index3==1):
        O1_abund_fake=np.sqrt(O1a_abund_fake**2+O1b_abund_fake**2+O1c_abund_fake**2)/3.
        O1_abund_fake2=(O1a_abund_fake1+O1b_abund_fake2+O1c_abund_fake3)/3.
        O1_abund=(O1a_abund+O1b_abund+O1c_abund)/3.
    
    


    if (O1a_abund[0]<0.01 and O1b_abund[0]<0.01 and O1c_abund[0]<0.01):
        if(O1a_abund[0]>=0  and O1b_abund[0]>=0  and O1c_abund[0]>=0 ):
            ion_abun.OIa[i,j]=O1a_abund[0]
            ion_abun.OIb[i,j]=O1b_abund[0]
            ion_abun.OIc[i,j]=O1c_abund[0]
            if np.isnan(O1a_abund_fake)==True: O1a_abund_fake=0
            ion_abun_error.OIa[i]=O1a_abund_fake
            if np.isnan(O1b_abund_fake)==True: O1b_abund_fake=0
            ion_abun_error.OIb[i]=O1b_abund_fake
            if np.isnan(O1c_abund_fake)==True: O1c_abund_fake=0
            ion_abun_error.OIc[i]=O1c_abund_fake
    
#    print >>file10, 'Oo(5577)/H+ = {:5.2e}'.format(ion_abun.OIa[i,j]),'Oo(6300)/H+ = {:5.2e}'.format(ion_abun.OIb[i,j]),'Oo(6363)/H+ = {:5.2e}'.format(ion_abun.OIc[i,j])
#    print >>file10, 'Oo(5577)/H+ = {:5.2e}'.format(ion_abun_error.OIa[i,j]),'Oo(6300)/H+ = {:5.2e}'.format(ion_abun_error.OIb[i,j]),'Oo(6363)/H+ = {:5.2e}'.format(ion_abun_error.OIc[i,j])
#############################################################

#############################################################
    O2 = pn.Atom('O',2)
    
    index_name=param_estimated.index("O+_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

            
            
    index1=index2=index3=index4=0    
    index_line=line_names.index("O2_3727s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O2a_abund_fake = O2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O2_3726A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(3726)', Hbeta=100.0)
        O2a_abund = O2.getIonAbundance(int_ratio=obs.getIntens(0)['O2_3726A'],tem=Teinput, den=Neinput, to_eval='L(3726)', Hbeta=100.0)
        index1=1
        
        O2a_abund_fake1=np.nan_to_num(O2a_abund_fake)
        O2a_abund_fake1[O2a_abund_fake1 == inf] =  0

        O2a_abund_fake[O2a_abund_fake == inf] =  0
        O2a_abund_fake=np.nan_to_num(O2a_abund_fake)
        O2a_abund_fake = [i2 for i2 in O2a_abund_fake if i2 != 0]
        O2a_abund_fake = [i2 for i2 in O2a_abund_fake if i2 < O2a_abund*3 and i2>O2a_abund*0.3]
        O2a_abund_fake=np.std(np.nan_to_num(O2a_abund_fake), ddof=1)
    else:
        O2a_abund_fake=0.0
        O2a_abund=[0.0,0.0]
        
    index_line=line_names.index("O2_3729s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O2b_abund_fake = O2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O2_3729A'], tem=Teinput_fake, den=Neinput_fake, to_eval='L(3729)', Hbeta=100.0)
        O2b_abund = O2.getIonAbundance(int_ratio=obs.getIntens(0)['O2_3729A'], tem=Teinput, den=Neinput, to_eval='L(3729)', Hbeta=100.0)
        index2=1
        
        O2b_abund_fake2=np.nan_to_num(O2b_abund_fake)
        O2b_abund_fake2[O2b_abund_fake2 == inf] =  0

        O2b_abund_fake[O2b_abund_fake == inf] =  0
        O2b_abund_fake=np.nan_to_num(O2b_abund_fake)
        O2b_abund_fake = [i2 for i2 in O2b_abund_fake if i2 != 0]
        O2b_abund_fake = [i2 for i2 in O2b_abund_fake if i2 < O2b_abund*3 and i2>O2b_abund*0.3]
        O2b_abund_fake=np.std(np.nan_to_num(O2b_abund_fake), ddof=1)
    else:
        O2b_abund_fake=0.0
        O2b_abund=[0.0,0.0]
        
    index_line=line_names.index("O2_7320s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O2c_abund_fake = O2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O2_7319A+'],tem=Teinput_fake, den=Neinput_fake, to_eval='I(4,2)+I(5,2)', Hbeta=100.0)
        O2c_abund = O2.getIonAbundance(int_ratio=obs.getIntens(0)['O2_7319A+'],tem=Teinput, den=Neinput, to_eval='I(4,2)+I(5,2)', Hbeta=100.0)
        index3=1
        
        O2c_abund_fake3=np.nan_to_num(O2c_abund_fake)
        O2c_abund_fake3[O2c_abund_fake3 == inf] =  0

        O2c_abund_fake[O2c_abund_fake == inf] =  0
        O2c_abund_fake=np.nan_to_num(O2c_abund_fake)
        O2c_abund_fake = [i2 for i2 in O2c_abund_fake if i2 != 0]
        O2c_abund_fake = [i2 for i2 in O2c_abund_fake if i2 < O2c_abund*3 and i2>O2c_abund*0.3]
        O2c_abund_fake=np.std(np.nan_to_num(O2c_abund_fake), ddof=1)
    else:
        O2c_abund_fake=0.0
        O2c_abund=[0.0,0.0]
        
    index_line=line_names.index("O2_7330s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O2d_abund_fake = O2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O2_7330A+'],tem=Teinput_fake, den=Neinput_fake, to_eval='I(4,3)+I(5,3)', Hbeta=100.0)
        O2d_abund = O2.getIonAbundance(int_ratio=obs.getIntens(0)['O2_7330A+'],tem=Teinput, den=Neinput, to_eval='I(4,3)+I(5,3)', Hbeta=100.0)
        index4=1
        
        O2d_abund_fake4=np.nan_to_num(O2d_abund_fake)
        O2d_abund_fake4[O2d_abund_fake4 == inf] =  0

        O2d_abund_fake[O2d_abund_fake == inf] =  0
        O2d_abund_fake=np.nan_to_num(O2d_abund_fake)
        O2d_abund_fake = [i2 for i2 in O2d_abund_fake if i2 != 0]
        O2d_abund_fake = [i2 for i2 in O2d_abund_fake if i2 < O2d_abund*3 and i2>O2d_abund*0.3]
        O2d_abund_fake=np.std(np.nan_to_num(O2d_abund_fake), ddof=1)
    else:
        O2d_abund_fake=0.0
        O2d_abund=[0.0,0.0]
    
    #O2a_abund_fake1=np.nan_to_num(O2a_abund_fake)
    #O2b_abund_fake2=np.nan_to_num(O2b_abund_fake)
    #O2c_abund_fake3=np.nan_to_num(O2c_abund_fake)
    #O2d_abund_fake4=np.nan_to_num(O2d_abund_fake)
    
    #O2a_abund_fake=np.std(np.nan_to_num(O2a_abund_fake), ddof=1)
    #O2b_abund_fake=np.std(np.nan_to_num(O2b_abund_fake), ddof=1)
    #O2c_abund_fake=np.std(np.nan_to_num(O2c_abund_fake), ddof=1)
    #O2d_abund_fake=np.std(np.nan_to_num(O2d_abund_fake), ddof=1)
    
    O2_abund_fake=0    
    O2_abund_fake2=0    
    O2_abund=0

    if (index1==1 and index2==0 and index3==0 and index4==0):
        O2_abund_fake=O2a_abund_fake
        O2_abund_fake2=O2a_abund_fake1
        O2_abund=O2a_abund
    elif (index1==0 and index2==1 and index3==0 and index4==0):
        O2_abund_fake=O2a_abund_fake
        O2_abund_fake2=O2a_abund_fake2
        O2_abund=O2a_abund
    elif (index1==1 and index2==1 and index3==0 and index4==0):
        O2_abund_fake=np.sqrt(O2a_abund_fake**2+O2b_abund_fake**2)/2.
        O2_abund_fake2=(O2a_abund_fake1+O2b_abund_fake2)/2.
        O2_abund=(O2a_abund+O2b_abund)/2.
    elif (index1==0 and index2==0 and index3==1 and index4==0):
        O2_abund_fake=O2c_abund_fake
        O2_abund_fake2=O2c_abund_fake3
        O2_abund=O2c_abund
    elif (index1==0 and index2==0 and index3==0 and index4==1):
        O2_abund_fake=O2d_abund_fake
        O2_abund_fake2=O2d_abund_fake4
        O2_abund=O2d_abund
    elif (index1==0 and index2==0 and index3==1 and index4==1):
        O2_abund_fake=np.sqrt(O2c_abund_fake**2+O2d_abund_fake**2)/2.
        O2_abund_fake2=(O2c_abund_fake3+O2d_abund_fake4)/2.
        O2_abund=(O2c_abund+O2d_abund)/2.
    elif (index1==1 and index2==1 and index3==1 and index4==1):    
        O2_abund_fake=np.sqrt(O2a_abund_fake**2+ O2b_abund_fake**2+ O2c_abund_fake**2+ O2d_abund_fake**2)/4.
        O2_abund_fake2=(O2a_abund_fake1+ O2b_abund_fake2+ O2c_abund_fake3+ O2d_abund_fake4)/4.
        O2_abund=(O2a_abund+ O2b_abund+ O2c_abund+ O2d_abund)/4.
    
    
    
    if (O2a_abund[0]<0.01 and O2b_abund[0]<0.01 and O2c_abund[0]<0.01 and O2d_abund[0]<0.01):
        if(O2a_abund[0]>=0  and O2b_abund[0]>=0  and O2c_abund[0]>=0  and O2d_abund[0]>=0 ):         
            ion_abun.OIIa[i,j]=O2a_abund[0]
            ion_abun.OIIb[i,j]=O2b_abund[0]
            ion_abun.OIIc[i,j]=O2c_abund[0]
            ion_abun.OIId[i,j]=O2d_abund[0]
            if np.isnan(O2a_abund_fake)==True: O2a_abund_fake=0
            ion_abun_error.OIIa[i]=O2a_abund_fake
            if np.isnan(O2b_abund_fake)==True: O2b_abund_fake=0
            ion_abun_error.OIIb[i]=O2b_abund_fake
            if np.isnan(O2c_abund_fake)==True: O2c_abund_fake=0
            ion_abun_error.OIIc[i]=O2c_abund_fake
            if np.isnan(O2d_abund_fake)==True: O2d_abund_fake=0
            ion_abun_error.OIId[i]=O2d_abund_fake
    
#    print >>file10, 'O+(3727)/H+ = {:5.2e}'.format(ion_abun.OIIa[i,j]),'O+(3729)/H+ = {:5.2e}'.format(ion_abun.OIIb[i,j]),'O+(7320)/H+ = {:5.2e}'.format(ion_abun.OIIc[i,j]),'O+(7330)/H+ = {:5.2e}'.format(ion_abun.OIId[i,j])
#    print >>file10, 'O+(3727)/H+ = {:5.2e}'.format(ion_abun_error.OIIa[i,j]),'O+(3729)/H+ = {:5.2e}'.format(ion_abun_error.OIIb[i,j]),'O+(7320)/H+ = {:5.2e}'.format(ion_abun_error.OIIc[i,j]),'O+(7330)/H+ = {:5.2e}'.format(ion_abun_error.OIId[i,j])

#############################################################

#############################################################
    O3 = pn.Atom('O',3)

    index_name=param_estimated.index("O++_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

    
    index1=index2=index3=0    
    index_line=line_names.index("O3_4363s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O3a_abund_fake = O3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O3_4363A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(4363)', Hbeta=100.0)
        O3a_abund = O3.getIonAbundance(int_ratio=obs.getIntens(0)['O3_4363A'],tem=Teinput, den=Neinput, to_eval='L(4363)', Hbeta=100.0)
        index1=1
        
        O3a_abund_fake1=np.nan_to_num(O3a_abund_fake)
        O3a_abund_fake1[O3a_abund_fake1 == inf] =  0

        O3a_abund_fake[O3a_abund_fake == inf] =  0
        O3a_abund_fake=np.nan_to_num(O3a_abund_fake)
        O3a_abund_fake = [i2 for i2 in O3a_abund_fake if i2 != 0]
        O3a_abund_fake = [i2 for i2 in O3a_abund_fake if i2 < O3a_abund*3 and i2>O3a_abund*0.3]
        O3a_abund_fake=np.std(np.nan_to_num(O3a_abund_fake), ddof=1)
    else:
        O3a_abund_fake=0.0
        O3a_abund=[0.0,0.0]
        
    index_line=line_names.index("O3_4959s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O3b_abund_fake = O3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O3_4959A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(4959)', Hbeta=100.0)
        O3b_abund = O3.getIonAbundance(int_ratio=obs.getIntens(0)['O3_4959A'],tem=Teinput, den=Neinput, to_eval='L(4959)', Hbeta=100.0)
        index2=1
        
        O3b_abund_fake2=np.nan_to_num(O3b_abund_fake)
        O3b_abund_fake2[O3b_abund_fake2 == inf] =  0

        O3b_abund_fake[O3b_abund_fake == inf] =  0
        O3b_abund_fake=np.nan_to_num(O3b_abund_fake)
        O3b_abund_fake = [i2 for i2 in O3b_abund_fake if i2 != 0]
        O3b_abund_fake = [i2 for i2 in O3b_abund_fake if i2 < O3b_abund*3 and i2>O3b_abund*0.3]
        O3b_abund_fake=np.std(np.nan_to_num(O3b_abund_fake), ddof=1)
    else:
        O3b_abund_fake=0.0
        O3b_abund=[0.0,0.0]
        
    index_line=line_names.index("O3_5007s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        O3c_abund_fake = O3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['O3_5007A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(5007)', Hbeta=100.0)
        O3c_abund = O3.getIonAbundance(int_ratio=obs.getIntens(0)['O3_5007A'],tem=Teinput, den=Neinput, to_eval='L(5007)', Hbeta=100.0)
        index3=1
        
        O3c_abund_fake3=np.nan_to_num(O3c_abund_fake)
        O3c_abund_fake3[O3c_abund_fake3 == inf] =  0

        O3c_abund_fake[O3c_abund_fake == inf] =  0
        O3c_abund_fake=np.nan_to_num(O3c_abund_fake)
        O3c_abund_fake = [i2 for i2 in O3c_abund_fake if i2 != 0]
        O3c_abund_fake = [i2 for i2 in O3c_abund_fake if i2 < O3c_abund*3 and i2>O3c_abund*0.3]
        O3c_abund_fake=np.std(np.nan_to_num(O3c_abund_fake), ddof=1)
    else:
        O3c_abund_fake=0.0
        O3c_abund=[0.0,0.0]


    #O3a_abund_fake1=np.nan_to_num(O3a_abund_fake)
    #O3b_abund_fake2=np.nan_to_num(O3b_abund_fake)
    #O3c_abund_fake3=np.nan_to_num(O3c_abund_fake)

    #O3a_abund_fake=np.std(np.nan_to_num(O3a_abund_fake), ddof=1)
    #O3b_abund_fake=np.std(np.nan_to_num(O3b_abund_fake), ddof=1)
    #O3c_abund_fake=np.std(np.nan_to_num(O3c_abund_fake), ddof=1)
    
    
    O3_abund_fake=0 
    O3_abund_fake2=0 
    O3_abund=0 

    if (index1==0 and index2==0 and index3==1):
        O3_abund_fake=O3c_abund_fake
        O3_abund_fake2=O3c_abund_fake3
        O3_abund=O3c_abund
    elif (index1==0 and index2==1 and index3==0 ):
        O3_abund_fake=O3b_abund_fake
        O3_abund_fake2=O3b_abund_fake2
        O3_abund=O3b_abund
    elif (index1==1 and index2==1 and index3==1 ):
        O3_abund_fake=np.sqrt(O3a_abund_fake**2+O3b_abund_fake**2+O3c_abund_fake**2)/3.
        O3_abund_fake2=(O3a_abund_fake1+O3b_abund_fake2+O3c_abund_fake3)/3.
        O3_abund=(O3a_abund+O3b_abund+O3c_abund)/3.
    elif (index1==0 and index2==1 and index3==1 ):
        O3_abund_fake=np.sqrt(O3b_abund_fake**2+O3c_abund_fake**2)/2.
        O3_abund_fake2=(O3b_abund_fake2+O3c_abund_fake3)/2.
        O3_abund=(O3b_abund+O3c_abund)/2.


    
    if (O3a_abund[0]<0.01 and O3b_abund[0]<0.01 and O3c_abund[0]<0.01):
        if (O3a_abund[0]>=0  and O3b_abund[0]>=0  and O3c_abund[0]>=0 ):  
            ion_abun.OIIIa[i,j]=O3a_abund[0]
            ion_abun.OIIIb[i,j]=O3b_abund[0]
            ion_abun.OIIIc[i,j]=O3c_abund[0]
            if np.isnan(O3a_abund_fake)==True: O3a_abund_fake=0
            ion_abun_error.OIIIa[i]=O3a_abund_fake
            if np.isnan(O3b_abund_fake)==True: O3b_abund_fake=0
            ion_abun_error.OIIIb[i]=O3b_abund_fake
            if np.isnan(O3c_abund_fake)==True: O3c_abund_fake=0
            ion_abun_error.OIIIc[i]=O3c_abund_fake
    
#    print >>file10, 'O+2(4363)/H+ = {:5.2e}'.format(ion_abun.OIIIa[i,j]),'O+2(4959)/H+ = {:5.2e}'.format(ion_abun.OIIIb[i,j]),'O+2(5007)/H+ = {:5.2e}'.format(ion_abun.OIIIc[i,j])
#    print >>file10, 'O+2(4363)/H+ = {:5.2e}'.format(ion_abun_error.OIIIa[i,j]),'O+2(4959)/H+ = {:5.2e}'.format(ion_abun_error.OIIIb[i,j]),'O+2(5007)/H+ = {:5.2e}'.format(ion_abun_error.OIIIc[i,j])
#############################################################

#############################################################
    N1 = pn.Atom('N',1)
    
    index_name=param_estimated.index("No_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

        
#    index_line=line_names.index("N1_5199s")
#    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
#        N1a_abund = N1.getIonAbundance(int_ratio=obs.getIntens(0)['N1_5198A'],tem=Teinput, den=Neinput, to_eval='L(5198)', Hbeta=100.0)
#    else:
#        N1a_abund=[0.0,0.0]
    
    index1=0
    index_line=line_names.index("N1_5199s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        N1b_abund_fake = N1.getIonAbundance(int_ratio=obs_fake.getIntens(0)['N1_5200A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(5200)', Hbeta=100.0)
        N1b_abund = N1.getIonAbundance(int_ratio=obs.getIntens(0)['N1_5200A'],tem=Teinput, den=Neinput, to_eval='L(5200)', Hbeta=100.0)
        index1=1
        
        N1b_abund_fake1=np.nan_to_num(N1b_abund_fake)
        N1b_abund_fake1[N1b_abund_fake1 == inf] =  0

        N1b_abund_fake[N1b_abund_fake == inf] =  0
        N1b_abund_fake=np.nan_to_num(N1b_abund_fake)
        N1b_abund_fake = [i2 for i2 in N1b_abund_fake if i2 != 0]
        N1b_abund_fake = [i2 for i2 in N1b_abund_fake if i2 < N1b_abund*3 and i2>N1b_abund*0.3]
        N1b_abund_fake=np.std(np.nan_to_num(N1b_abund_fake), ddof=1)
    else:
        N1b_abund_fake=0.0
        N1b_abund=[0.0,0.0]
    
    #N1b_abund_fake1=np.nan_to_num(N1b_abund_fake)
    #N1b_abund_fake=np.std(np.nan_to_num(N1b_abund_fake), ddof=1)
    
    N1_abund_fake=0 
    N1_abund_fake2=0 
    N1_abund=0 
    if (index1==1):
        N1_abund_fake=N1b_abund_fake
        N1_abund_fake2=N1b_abund_fake1
        N1_abund=N1b_abund
    
    
    if (N1b_abund[0]<0.001):
        if(N1b_abund[0]>=0 ):
#           ion_abun.NIa[i,j]=N1a_abund
            ion_abun.NI[i,j]=N1b_abund[0]
            if np.isnan(N1b_abund_fake)==True: N1b_abund_fake=0
            ion_abun_error.NI[i]=N1b_abund_fake

#    print('No(5199)/H+ = {:5.2e}'.format(N1a_abund),'No(5200)/H+ = {:5.2e}'.format(N1b_abund))
#    print >>file10, 'No(5200)/H+ = {:5.2e}'.format(ion_abun.NI[i,j])
#    print >>file10, 'No(5200)/H+ = {:5.2e}'.format(ion_abun_error.NI[i,j])

#############################################################


#############################################################
    N2 = pn.Atom('N',2)
    
    index_name=param_estimated.index("N+_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

        
    index1=index2=index3=0    
    index_line=line_names.index("N2_5755s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        N2a_abund_fake = N2.getIonAbundance(int_ratio=obs_fake.getIntens()['N2_5755A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(5755)',   Hbeta=100.0)
        N2a_abund = N2.getIonAbundance(int_ratio=obs.getIntens()['N2_5755A'],tem=Teinput, den=Neinput, to_eval='L(5755)',   Hbeta=100.0)
        index1=1
        
        N2a_abund_fake1=np.nan_to_num(N2a_abund_fake)
        N2a_abund_fake1[N2a_abund_fake1 == inf] =  0

        N2a_abund_fake[N2a_abund_fake == inf] =  0
        N2a_abund_fake=np.nan_to_num(N2a_abund_fake)
        N2a_abund_fake = [i2 for i2 in N2a_abund_fake if i2 != 0]
        N2a_abund_fake = [i2 for i2 in N2a_abund_fake if i2 < N2a_abund*3 and i2>N2a_abund*0.3]
        N2a_abund_fake=np.std(np.nan_to_num(N2a_abund_fake), ddof=1)
    else:
        N2a_abund_fake=0.0
        N2a_abund=[0.0,0.0]
        
    index_line=line_names.index("N2_6548s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        N2b_abund_fake = N2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['N2_6548A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(6548)', Hbeta=100.0)
        N2b_abund = N2.getIonAbundance(int_ratio=obs.getIntens(0)['N2_6548A'],tem=Teinput, den=Neinput, to_eval='L(6548)', Hbeta=100.0)
        index2=1
        
        N2b_abund_fake2=np.nan_to_num(N2b_abund_fake)
        N2b_abund_fake2[N2b_abund_fake2 == inf] =  0

        N2b_abund_fake[N2b_abund_fake == inf] =  0
        N2b_abund_fake=np.nan_to_num(N2b_abund_fake)
        N2b_abund_fake = [i2 for i2 in N2b_abund_fake if i2 != 0]
        N2b_abund_fake = [i2 for i2 in N2b_abund_fake if i2 < N2b_abund*3 and i2>N2b_abund*0.3]
        N2b_abund_fake=np.std(np.nan_to_num(N2b_abund_fake), ddof=1)
    else:
        N2b_abund_fake=0.0
        N2b_abund=[0.0,0.0]
    
    index_line=line_names.index("N2_6583s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        N2c_abund_fake = N2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['N2_6584A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(6584)', Hbeta=100.0)      
        N2c_abund = N2.getIonAbundance(int_ratio=obs.getIntens(0)['N2_6584A'],tem=Teinput, den=Neinput, to_eval='L(6584)', Hbeta=100.0)
        index3=1
        
        N2c_abund_fake3=np.nan_to_num(N2c_abund_fake)
        N2c_abund_fake3[N2c_abund_fake3 == inf] =  0

        N2c_abund_fake[N2c_abund_fake == inf] =  0
        N2c_abund_fake=np.nan_to_num(N2c_abund_fake)
        N2c_abund_fake = [i2 for i2 in N2c_abund_fake if i2 != 0]
        N2c_abund_fake = [i2 for i2 in N2c_abund_fake if i2 < N2c_abund*3 and i2>N2c_abund*0.3]
        N2c_abund_fake=np.std(np.nan_to_num(N2c_abund_fake), ddof=1)
    else:
        N2c_abund_fake=0.0
        N2c_abund=[0.0,0.0]
        
    #N2a_abund_fake1=np.nan_to_num(N2a_abund_fake)
    #N2b_abund_fake2=np.nan_to_num(N2b_abund_fake)
    #N2c_abund_fake3=np.nan_to_num(N2c_abund_fake)
    
    #N2a_abund_fake=np.std(np.nan_to_num(N2a_abund_fake), ddof=1)
    #N2b_abund_fake=np.std(np.nan_to_num(N2b_abund_fake), ddof=1)
    #N2c_abund_fake=np.std(np.nan_to_num(N2c_abund_fake), ddof=1)
    
    N2_abund_fake=0 
    N2_abund_fake2=0 
    N2_abund=0 
    
    if (index1==0 and index2==0 and index3==1):
        N2_abund_fake=N2c_abund_fake
        N2_abund_fake2=N2c_abund_fake3
        N2_abund=N2c_abund
    elif (index1==0 and index2==1 and index3==0 ):
        N2_abund_fake=N2b_abund_fake
        N2_abund_fake2=N2b_abund_fake2
        N2_abund=N2b_abund
    elif (index1==1 and index2==1 and index3==1 ):
        N2_abund_fake=np.sqrt(N2a_abund_fake**2+N2b_abund_fake**2+N2c_abund_fake**2)/3.
        N2_abund_fake2=(N2a_abund_fake1+N2b_abund_fake2+N2c_abund_fake3)/3.
        N2_abund=(N2a_abund+N2b_abund+N2c_abund)/3.
    elif (index1==0 and index2==1 and index3==1 ):
        N2_abund_fake=np.sqrt(N2b_abund_fake**2+N2c_abund_fake**2)/2.     
        N2_abund_fake2=(N2b_abund_fake2+N2c_abund_fake3)/2.     
        N2_abund=(N2b_abund+N2c_abund)/2.
 
    if (N2a_abund[0]<0.001 and N2b_abund[0]<0.001 and N2c_abund[0]<0.001):
        if(N2a_abund[0]>=0  and N2b_abund[0]>=0  and N2c_abund[0]>=0 ): 
            ion_abun.NIIa[i,j]=N2a_abund[0]
            ion_abun.NIIb[i,j]=N2b_abund[0]
            ion_abun.NIIc[i,j]=N2c_abund[0]
            if np.isnan(N2a_abund_fake)==True: N2a_abund_fake=0
            ion_abun_error.NIIa[i]=N2a_abund_fake
            if np.isnan(N2b_abund_fake)==True: N2b_abund_fake=0
            ion_abun_error.NIIb[i]=N2b_abund_fake
            if np.isnan(N2c_abund_fake)==True: N2c_abund_fake=0
            ion_abun_error.NIIc[i]=N2c_abund_fake

#    print >>file10, 'N+(5755)/H+ = {:5.2e}'.format(ion_abun.NIIa[i,j]),'N+(6548)/H+ = {:5.2e}'.format(ion_abun.NIIb[i,j]),'N+(6584)/H+ = {:5.2e}'.format(ion_abun.NIIc[i,j])
#    print >>file10, 'N+(5755)/H+ = {:5.2e}'.format(ion_abun_error.NIIa[i,j]),'N+(6548)/H+ = {:5.2e}'.format(ion_abun_error.NIIb[i,j]),'N+(6584)/H+ = {:5.2e}'.format(ion_abun_error.NIIc[i,j])
#############################################################
     
#############################################################
    S2 = pn.Atom('S',2)
    
    index_name=param_estimated.index("S+_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

    
    index1=index2=0    
    index_line=line_names.index("S2_6716s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        S2a_abund_fake = S2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['S2_6716A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(6717)', Hbeta=100.0)
        S2a_abund = S2.getIonAbundance(int_ratio=obs.getIntens(0)['S2_6716A'],tem=Teinput, den=Neinput, to_eval='L(6717)', Hbeta=100.0)
        index1=1
        
        S2a_abund_fake1=np.nan_to_num(S2a_abund_fake)
        S2a_abund_fake1[S2a_abund_fake1 == inf] =  0

        S2a_abund_fake[S2a_abund_fake == inf] =  0
        S2a_abund_fake=np.nan_to_num(S2a_abund_fake)
        S2a_abund_fake = [i2 for i2 in S2a_abund_fake if i2 != 0]
        S2a_abund_fake = [i2 for i2 in S2a_abund_fake if i2 < S2a_abund*3 and i2>S2a_abund*0.3]
        S2a_abund_fake=np.std(np.nan_to_num(S2a_abund_fake), ddof=1)
    else:
        S2a_abund_fake=0.0
        S2a_abund=[0.0,0.0]
        
    index_line=line_names.index("S2_6731s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        S2b_abund_fake = S2.getIonAbundance(int_ratio=obs_fake.getIntens(0)['S2_6731A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(6731)', Hbeta=100.0)
        S2b_abund = S2.getIonAbundance(int_ratio=obs.getIntens(0)['S2_6731A'],tem=Teinput, den=Neinput, to_eval='L(6731)', Hbeta=100.0)
        index2=1
        
        S2b_abund_fake2=np.nan_to_num(S2b_abund_fake)
        S2b_abund_fake2[S2b_abund_fake2 == inf] =  0

        S2b_abund_fake[S2b_abund_fake == inf] =  0
        S2b_abund_fake=np.nan_to_num(S2b_abund_fake)
        S2b_abund_fake = [i2 for i2 in S2b_abund_fake if i2 != 0]
        S2b_abund_fake = [i2 for i2 in S2b_abund_fake if i2 < S2b_abund*3 and i2>S2b_abund*0.3]
        S2b_abund_fake=np.std(np.nan_to_num(S2b_abund_fake), ddof=1)
    else:
        S2b_abund_fake=0.0
        S2b_abund=[0.0,0.0]
    
    
    #S2a_abund_fake1=np.nan_to_num(S2a_abund_fake)
    #S2b_abund_fake2=np.nan_to_num(S2b_abund_fake)
    
    #S2a_abund_fake=np.std(np.nan_to_num(S2a_abund_fake), ddof=1)
    #S2b_abund_fake=np.std(np.nan_to_num(S2b_abund_fake), ddof=1)
    
    S2_abund_fake=0 
    S2_abund_fake2=0 
    S2_abund=0 
    if (index1==1 and index2==0):
        S2_abund_fake=S2a_abund_fake
        S2_abund_fake2=S2a_abund_fake1
        S2_abund=S2a_abund
    elif (index1==0 and index2==1):
        S2_abund_fake=S2b_abund_fake
        S2_abund_fake2=S2b_abund_fake2
        S2_abund=S2b_abund
    elif (index1==1 and index2==1):
        S2_abund_fake=np.sqrt(S2a_abund_fake**2+S2b_abund_fake**2)/2.
        S2_abund_fake2=(S2a_abund_fake1+S2b_abund_fake2)/2.
        S2_abund=(S2a_abund+S2b_abund)/2.
    

    if (S2a_abund[0]<0.01 and S2b_abund[0]<0.01):
        if(S2a_abund[0]>=0  and S2b_abund[0]>=0 ): 
            ion_abun.SIIa[i,j]=S2a_abund[0]
            ion_abun.SIIb[i,j]=S2b_abund[0]
            if np.isnan(S2a_abund_fake)==True: S2a_abund_fake=0
            ion_abun_error.SIIa[i]=S2a_abund_fake
            if np.isnan(S2b_abund_fake)==True: S2b_abund_fake=0
            ion_abun_error.SIIb[i]=S2b_abund_fake
    
#    print >>file10, 'S+(6717)/H+ = {:5.2e}'.format(ion_abun.SIIa[i,j]),'S+(6731)/H+ = {:5.2e}'.format(ion_abun.SIIb[i,j])
#    print >>file10, 'S+(6717)/H+ = {:5.2e}'.format(ion_abun_error.SIIa[i,j]),'S+(6731)/H+ = {:5.2e}'.format(ion_abun_error.SIIb[i,j])
#############################################################     
    
#############################################################
    S3 = pn.Atom('S',3)
    
    index_name=param_estimated.index("S++_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j]
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

    
    index1=index2=0    
    index_line=line_names.index("S3_6312s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        S3a_abund_fake = S3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['S3_6312A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(6312)', Hbeta=100.0)
        S3a_abund = S3.getIonAbundance(int_ratio=obs.getIntens(0)['S3_6312A'],tem=Teinput, den=Neinput, to_eval='L(6312)', Hbeta=100.0)
        index1=1
        
        S3a_abund_fake1=np.nan_to_num(S3a_abund_fake)
        S3a_abund_fake1[S3a_abund_fake1 == inf] =  0

        S3a_abund_fake[S3a_abund_fake == inf] =  0
        S3a_abund_fake=np.nan_to_num(S3a_abund_fake)
        S3a_abund_fake = [i2 for i2 in S3a_abund_fake if i2 != 0]
        S3a_abund_fake = [i2 for i2 in S3a_abund_fake if i2 < S3a_abund*3 and i2>S3a_abund*0.3]
        S3a_abund_fake=np.std(np.nan_to_num(S3a_abund_fake), ddof=1)
    else:
        S3a_abund_fake=0.0
        S3a_abund=[0.0,0.0]
        
    index_line=line_names.index("S3_9069s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        S3b_abund_fake = S3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['S3_9069A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(9069)', Hbeta=100.0)
        S3b_abund = S3.getIonAbundance(int_ratio=obs.getIntens(0)['S3_9069A'],tem=Teinput, den=Neinput, to_eval='L(9069)', Hbeta=100.0)
        index2=1
        
        S3b_abund_fake2=np.nan_to_num(S3b_abund_fake)
        S3b_abund_fake2[S3b_abund_fake2 == inf] =  0

        S3b_abund_fake[S3b_abund_fake == inf] =  0
        S3b_abund_fake=np.nan_to_num(S3b_abund_fake)
        S3b_abund_fake = [i2 for i2 in S3b_abund_fake if i2 != 0]
        S3b_abund_fake = [i2 for i2 in S3b_abund_fake if i2 < S3b_abund*3 and i2>S3b_abund*0.3] 
        S3b_abund_fake=np.std(np.nan_to_num(S3b_abund_fake), ddof=1)
    else:
        S3b_abund_fake=0.0
        S3b_abund=[0.0,0.0]

#    index_line=line_names.index("S3_9532s")
#    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
#        S3c_abund = S3.getIonAbundance(int_ratio=obs.getIntens(0)['S3_9531A'],tem=Teinput, den=Neinput, to_eval='L(9531)', Hbeta=100.0)
#    else:
#        S3c_abund=[0.0,0.0]
    
    #S3a_abund_fake1=np.nan_to_num(S3a_abund_fake)
    #S3b_abund_fake2=np.nan_to_num(S3b_abund_fake)
    
    #S3a_abund_fake=np.std(np.nan_to_num(S3a_abund_fake), ddof=1)
    #S3b_abund_fake=np.std(np.nan_to_num(S3b_abund_fake), ddof=1)
    
    S3_abund_fake=0 
    S3_abund_fake2=0 
    S3_abund=0 
    if (index1==1 and index2==0):
        S3_abund_fake=S3a_abund_fake
        S3_abund_fake2=S3a_abund_fake1
        S3_abund=S3a_abund
    if (index1==0 and index2==1):
        S3_abund_fake=S3b_abund_fake
        S3_abund_fake2=S3b_abund_fake2
        S3_abund=S3b_abund
    elif (index1==1 and index2==1):
        S3_abund_fake=np.sqrt(S3a_abund_fake**2+S3b_abund_fake**2)/2.
        S3_abund_fake2=(S3a_abund_fake1+S3b_abund_fake2)/2.
        S3_abund=(S3a_abund+S3b_abund)/2.
    
    if (S3a_abund[0]<0.01 and S3b_abund[0]<0.01):
        if(S3a_abund[0]>=0  and S3b_abund[0]>=0 ):    
            ion_abun.SIIIa[i,j]=S3a_abund[0]
            ion_abun.SIIIb[i,j]=S3b_abund[0]
            if np.isnan(S3a_abund_fake)==True: S3a_abund_fake=0
            ion_abun_error.SIIIa[i]=S3a_abund_fake
            if np.isnan(S3b_abund_fake)==True: S3b_abund_fake=0
            ion_abun_error.SIIIb[i]=S3b_abund_fake
                                             
#    print('S+2(6312)/H+ = {:5.2e}'.format(S3a_abund),'S+2(9068)/H+ = {:5.2e}'.format(S3b_abund),'S+2(9532)/H+ = {:5.2e}'.format(S3c_abund))
#    print >>file10, 'S+2(6312)/H+ = {:5.2e}'.format(ion_abun.SIIIa[i,j]),'S+2(9068)/H+ = {:5.2e}'.format(ion_abun.SIIIb[i,j])
#    print >>file10, 'S+2(6312)/H+ = {:5.2e}'.format(ion_abun_error.SIIIa[i,j]),'S+2(9068)/H+ = {:5.2e}'.format(ion_abun_error.SIIIb[i,j])

#############################################################   

   
#############################################################
    Ne3 = pn.Atom('Ne',3)
    
    index_name=param_estimated.index("Ne++_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j]
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

            
    index1=index2=0
    index_line=line_names.index("Ne3_3868s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        Ne3a_abund_fake = Ne3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['Ne3_3869A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(3869)', Hbeta=100.0)
        Ne3a_abund = Ne3.getIonAbundance(int_ratio=obs.getIntens(0)['Ne3_3869A'],tem=Teinput, den=Neinput, to_eval='L(3869)', Hbeta=100.0)
        index1=1
        
        Ne3a_abund_fake1=np.nan_to_num(Ne3a_abund_fake)
        Ne3a_abund_fake1[Ne3a_abund_fake1 == inf] =  0

        Ne3a_abund_fake[Ne3a_abund_fake == inf] =  0
        Ne3a_abund_fake=np.nan_to_num(Ne3a_abund_fake)
        Ne3a_abund_fake = [i2 for i2 in Ne3a_abund_fake if i2 != 0]
        Ne3a_abund_fake = [i2 for i2 in Ne3a_abund_fake if i2 < Ne3a_abund*3 and i2>Ne3a_abund*0.3]
        Ne3a_abund_fake=np.std(np.nan_to_num(Ne3a_abund_fake), ddof=1)
    else:
        Ne3a_abund_fake=0.0
        Ne3a_abund=[0.0,0.0]
        
    index_line=line_names.index("Ne3_3967s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        Ne3b_abund_fake = Ne3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['Ne3_3968A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(3968)', Hbeta=100.0)
        Ne3b_abund = Ne3.getIonAbundance(int_ratio=obs.getIntens(0)['Ne3_3968A'],tem=Teinput, den=Neinput, to_eval='L(3968)', Hbeta=100.0)
        index2=1
        
        Ne3b_abund_fake2=np.nan_to_num(Ne3b_abund_fake)
        Ne3b_abund_fake2[Ne3b_abund_fake2 == inf] =  0

        Ne3b_abund_fake[Ne3b_abund_fake == inf] =  0
        Ne3b_abund_fake=np.nan_to_num(Ne3b_abund_fake)
        Ne3b_abund_fake = [i2 for i2 in Ne3b_abund_fake if i2 != 0]
        Ne3b_abund_fake = [i2 for i2 in Ne3b_abund_fake if i2 < Ne3b_abund*3 and i2>Ne3b_abund*0.3]
        Ne3b_abund_fake=np.std(np.nan_to_num(Ne3b_abund_fake), ddof=1)
    else:
        Ne3b_abund_fake=0.0
        Ne3b_abund=[0.0,0.0]

    #Ne3a_abund_fake1=np.nan_to_num(Ne3a_abund_fake)
    #Ne3b_abund_fake2=np.nan_to_num(Ne3b_abund_fake)
    
    #Ne3a_abund_fake=np.std(np.nan_to_num(Ne3a_abund_fake), ddof=1)
    #Ne3b_abund_fake=np.std(np.nan_to_num(Ne3b_abund_fake), ddof=1)

    Ne3_abund_fake=0 
    Ne3_abund_fake2=0 
    Ne3_abund=0 
    if (index1==1 and index2==0):
        Ne3_abund_fake=Ne3a_abund_fake
        Ne3_abund_fake2=Ne3a_abund_fake1
        Ne3_abund=Ne3a_abund
    elif (index1==0 and index2==1):
        Ne3_abund_fake=Ne3b_abund_fake
        Ne3_abund_fake2=Ne3b_abund_fake2
        Ne3_abund=Ne3b_abund
    elif (index1==1 and index2==1):
        Ne3_abund_fake=np.sqrt(Ne3a_abund_fake**2+Ne3b_abund_fake**2)/2.
        Ne3_abund_fake2=(Ne3a_abund_fake1+Ne3b_abund_fake2)/2.
        Ne3_abund=(Ne3a_abund+Ne3b_abund)/2.

    if (Ne3a_abund[0]<0.01 and Ne3b_abund[0]<0.01):
        if(Ne3a_abund[0]>=0  and Ne3b_abund[0]>=0 ):  
            ion_abun.NeIIIa[i,j]=Ne3a_abund[0]
            ion_abun.NeIIIb[i,j]=Ne3b_abund[0]     
            if np.isnan(Ne3a_abund_fake)==True: Ne3a_abund_fake=0
            ion_abun_error.NeIIIa[i]=Ne3a_abund_fake
            if np.isnan(Ne3b_abund_fake)==True: Ne3b_abund_fake=0
            ion_abun_error.NeIIIb[i]=Ne3b_abund_fake
    
#    print >>file10, 'Ne+2(3869)/H+ = {:5.2e}'.format(ion_abun.NeIIIa[i,j]),'Ne+2(3968)/H+ = {:5.2e}'.format(ion_abun.NeIIIb[i,j])
#    print >>file10, 'Ne+2(3869)/H+ = {:5.2e}'.format(ion_abun_error.NeIIIa[i,j]),'Ne+2(3968)/H+ = {:5.2e}'.format(ion_abun_error.NeIIIb[i,j])

#############################################################   


#############################################################
    Ar3 = pn.Atom('Ar',3)
    
    index_name=param_estimated.index("Ar++_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

    
    index1=0    
    index_line=line_names.index("Ar3_7136s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        Ar3a_abund_fake = Ar3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['Ar3_7136A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(7136)', Hbeta=100.0)
        Ar3a_abund = Ar3.getIonAbundance(int_ratio=obs.getIntens(0)['Ar3_7136A'],tem=Teinput, den=Neinput, to_eval='L(7136)', Hbeta=100.0)
        index1=1
        
        Ar3a_abund_fake1=np.nan_to_num(Ar3a_abund_fake)
        Ar3a_abund_fake1[Ar3a_abund_fake1 == inf] =  0

        Ar3a_abund_fake[Ar3a_abund_fake == inf] =  0
        Ar3a_abund_fake=np.nan_to_num(Ar3a_abund_fake)
        Ar3a_abund_fake = [i2 for i2 in Ar3a_abund_fake if i2 != 0]
        Ar3a_abund_fake = [i2 for i2 in Ar3a_abund_fake if i2 < Ar3a_abund*3 and i2>Ar3a_abund*0.3]
        Ar3a_abund_fake=np.std(np.nan_to_num(Ar3a_abund_fake), ddof=1)
    else:
        Ar3a_abund_fake=0.0
        Ar3a_abund=[0.0,0.0]

    #Ar3a_abund_fake1=np.nan_to_num(Ar3a_abund_fake)
    #Ar3a_abund_fake=np.std(np.nan_to_num(Ar3a_abund_fake), ddof=1)

    Ar3_abund_fake=0 
    Ar3_abund_fake2=0 
    Ar3_abund=0 
    if (index1==1):
        Ar3_abund_fake=Ar3a_abund_fake
        Ar3_abund_fake2=Ar3a_abund_fake1
        Ar3_abund=Ar3a_abund


    if (Ar3a_abund[0]<0.01):
        if(Ar3a_abund[0]>=0 ):
            ion_abun.ArIII[i,j]=Ar3a_abund[0]
            if np.isnan(Ar3a_abund_fake)==True: Ar3a_abund_fake=0
            ion_abun_error.ArIII[i]=Ar3a_abund_fake
    
#    print >>file10, 'Ar+2(7136)/H+ = {:5.2e}'.format(ion_abun.ArIII[i,j])
#    print >>file10, 'Ar+2(7136)/H+ = {:5.2e}'.format(ion_abun_error.ArIII[i,j])
#############################################################  

#############################################################
    Ar4 = pn.Atom('Ar',4)
        
    index_name=param_estimated.index("Ar+++_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

    
    index1=index2=0
    index_line=line_names.index("Ar4_4712s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        Ar4a_abund_fake = Ar4.getIonAbundance(int_ratio=obs_fake.getIntens(0)['Ar4_4711A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(4711)', Hbeta=100.0)
        Ar4a_abund = Ar4.getIonAbundance(int_ratio=obs.getIntens(0)['Ar4_4711A'],tem=Teinput, den=Neinput, to_eval='L(4711)', Hbeta=100.0)
        index1=1
        
        Ar4a_abund_fake1=np.nan_to_num(Ar4a_abund_fake)
        Ar4a_abund_fake1[Ar4a_abund_fake1 == inf] =  0

        Ar4a_abund_fake[Ar4a_abund_fake == inf] =  0
        Ar4a_abund_fake=np.nan_to_num(Ar4a_abund_fake)
        Ar4a_abund_fake = [i2 for i2 in Ar4a_abund_fake if i2 != 0]
        Ar4a_abund_fake = [i2 for i2 in Ar4a_abund_fake if i2 < Ar4a_abund*3 and i2>Ar4a_abund*0.3]
        Ar4a_abund_fake=np.std(np.nan_to_num(Ar4a_abund_fake), ddof=1)
    else:
        Ar4a_abund_fake=0.0
        Ar4a_abund=[0.0,0.0]
        
    index_line=line_names.index("Ar4_4740s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        Ar4b_abund_fake = Ar4.getIonAbundance(int_ratio=obs_fake.getIntens(0)['Ar4_4740A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(4740)', Hbeta=100.0)
        Ar4b_abund = Ar4.getIonAbundance(int_ratio=obs.getIntens(0)['Ar4_4740A'],tem=Teinput, den=Neinput, to_eval='L(4740)', Hbeta=100.0)
        index2=1
        
        Ar4b_abund_fake2=np.nan_to_num(Ar4b_abund_fake)
        Ar4b_abund_fake2[Ar4b_abund_fake2 == inf] =  0

        Ar4b_abund_fake[Ar4b_abund_fake == inf] =  0
        Ar4b_abund_fake=np.nan_to_num(Ar4b_abund_fake)
        Ar4b_abund_fake = [i2 for i2 in Ar4b_abund_fake if i2 != 0]
        Ar4b_abund_fake = [i2 for i2 in Ar4b_abund_fake if i2 < Ar4b_abund*3 and i2>Ar4b_abund*0.3]
        Ar4b_abund_fake=np.std(np.nan_to_num(Ar4b_abund_fake), ddof=1)
    else:
        Ar4b_abund_fake=0.0
        Ar4b_abund=[0.0,0.0]

    #Ar4a_abund_fake1=np.nan_to_num(Ar4a_abund_fake)
    #Ar4b_abund_fake2=np.nan_to_num(Ar4b_abund_fake)
    
    #Ar4a_abund_fake=np.std(np.nan_to_num(Ar4a_abund_fake), ddof=1)
    #Ar4b_abund_fake=np.std(np.nan_to_num(Ar4b_abund_fake), ddof=1)
    
    Ar4_abund_fake=0
    Ar4_abund_fake2=0
    Ar4_abund=0

    if (index2==1 and index2==0):
        Ar4_abund_fake=Ar4a_abund_fake
        Ar4_abund_fake2=Ar4a_abund_fake1
        Ar4_abund=Ar4a_abund
    elif (index2==0 and index2==1):
        Ar4_abund_fake=Ar4b_abund_fake
        Ar4_abund_fake2=Ar4b_abund_fake2
        Ar4_abund=Ar4b_abund
    elif (index1==1 and index2==1):
        Ar4_abund_fake=np.sqrt(Ar4a_abund_fake**2+Ar4b_abund_fake**2)/2.
        Ar4_abund_fake2=(Ar4a_abund_fake1+Ar4b_abund_fake2)/2.
        Ar4_abund=(Ar4a_abund+Ar4b_abund)/2.
    
    if (Ar4a_abund[0]<0.01 and Ar4b_abund[0]<0.01):
        if(Ar4a_abund[0]>=0  and Ar4b_abund[0]>=0 ):
            ion_abun.ArIVa[i,j]=Ar4a_abund[0]
            ion_abun.ArIVb[i,j]=Ar4b_abund[0]
            if np.isnan(Ar4a_abund_fake)==True: Ar4a_abund_fake=0
            ion_abun_error.ArIVa[i]=Ar4a_abund_fake
            if np.isnan(Ar4b_abund_fake)==True: Ar4b_abund_fake=0
            ion_abun_error.ArIVb[i]=Ar4b_abund_fake
    
#    print >>file10, 'Ar+3(4712)/H+ = {:5.2e}'.format(ion_abun.ArIVa[i,j]),'Ar+3(4740)/H+ = {:5.2e}'.format(ion_abun.ArIVb[i,j])
#    print >>file10, 'Ar+3(4712)/H+ = {:5.2e}'.format(ion_abun_error.ArIVa[i,j]),'Ar+3(4740)/H+ = {:5.2e}'.format(ion_abun_error.ArIVb[i,j])
#############################################################  

#############################################################
    Cl3 = pn.Atom('Cl',3)
    
    index_name=param_estimated.index("Cl++_Te_Ne")
    if (param_requered[index_name]=='TeNIINeSII'):
        Teinput=Te.NIISII[i,j]
        Neinput=Ne.NIISII[i,j]
        
        Teinput_fake=Te_fake.NIISII
        Neinput_fake=Ne_fake.NIISII
    if (param_requered[index_name]=='TeOINeSII'):
        Teinput=Te.OISII[i,j]
        Neinput=Ne.OISII[i,j]
        
        Teinput_fake=Te_fake.OISII
        Neinput_fake=Ne_fake.OISII
    if (param_requered[index_name]=='TeOIINeSII'):
        Teinput=Te.OIISII[i,j]
        Neinput=Ne.OIISII[i,j]
        
        Teinput_fake=Te_fake.OIISII
        Neinput_fake=Ne_fake.OIISII
    if (param_requered[index_name]=='TeOIIINeSII'):
        Teinput=Te.OIIISII[i,j]
        Neinput=Ne.OIIISII[i,j]
        
        Teinput_fake=Te_fake.OIIISII
        Neinput_fake=Ne_fake.OIIISII
    if (param_requered[index_name]=='TeSIIINeSII'):
        Teinput=Te.SIIISII[i,j]
        Neinput=Ne.SIIISII[i,j] 
        
        Teinput_fake=Te_fake.SIIISII
        Neinput_fake=Ne_fake.SIIISII
    if (param_requered[index_name]=='TeOIINeOII'):
        Teinput=Te.OIIOII[i,j]
        Neinput=Ne.OIIOII[i,j]
        
        Teinput_fake=Te_fake.OIIOII
        Neinput_fake=Ne_fake.OIIOII
    if (param_requered[index_name]=='TeNIINeOII'):
        Teinput=Te.NIIOII[i,j]
        Neinput=Ne.NIIOII[i,j]
        
        Teinput_fake=Te_fake.NIIOII
        Neinput_fake=Ne_fake.NIIOII
    if (param_requered[index_name]=='TeOINeOII'):
        Teinput=Te.OIOII[i,j]
        Neinput=Ne.OIOII[i,j]
        
        Teinput_fake=Te_fake.OIOII
        Neinput_fake=Ne_fake.OIOII
    if (param_requered[index_name]=='TeOIIINeClIII'):
        Teinput=Te.OIIIClIII[i,j]
        Neinput=Ne.OIIIClIII[i,j]
        
        Teinput_fake=Te_fake.OIIIClIII
        Neinput_fake=Ne_fake.OIIIClIII
    if (param_requered[index_name]=='TeSIIINeClIII'):
        Teinput=Te.SIIIClIII[i,j]
        Neinput=Ne.SIIIClIII[i,j]
        
        Teinput_fake=Te_fake.SIIIClIII
        Neinput_fake=Ne_fake.SIIIClIII
    if (param_requered[index_name]=='TeOIIINeArVI'):
        Teinput=Te.OIIIArVI[i,j]
        Neinput=Ne.OIIIArVI[i,j] 
        
        Teinput_fake=Te_fake.OIIIArVI
        Neinput_fake=Ne_fake.OIIIArVI
    if (param_requered[index_name]=='TeSIIINeArVI'):
        Teinput=Te.SIIIArVI[i,j]
        Neinput=Ne.SIIIArVI[i,j]
        
        Teinput_fake=Te_fake.SIIIArVI
        Neinput_fake=Ne_fake.SIIIArVI
    if (param_requered[index_name]=='TeNIINeClIII'):
        Teinput=Te.NIIClIII[i,j]
        Neinput=Ne.NIIClIII[i,j]
        
        Teinput_fake=Te_fake.NIIClIII
        Neinput_fake=Ne_fake.NIIClIII

    
    index1=index2=0
    index_line=line_names.index("Cl3_5517s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        Cl3a_abund_fake = Cl3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['Cl3_5518A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(5518)', Hbeta=100.0)
        Cl3a_abund = Cl3.getIonAbundance(int_ratio=obs.getIntens(0)['Cl3_5518A'],tem=Teinput, den=Neinput, to_eval='L(5518)', Hbeta=100.0)
        index1=1
        
        Cl3a_abund_fake1=np.nan_to_num(Cl3a_abund_fake)
        Cl3a_abund_fake1[Cl3a_abund_fake1 == inf] =  0

        Cl3a_abund_fake[Cl3a_abund_fake == inf] =  0
        Cl3a_abund_fake=np.nan_to_num(Cl3a_abund_fake)
        Cl3a_abund_fake = [i2 for i2 in Cl3a_abund_fake if i2 != 0]
        Cl3a_abund_fake = [i2 for i2 in Cl3a_abund_fake if i2 < Cl3a_abund*3 and i2>Cl3a_abund*0.3]
        Cl3a_abund_fake=np.std(np.nan_to_num(Cl3a_abund_fake), ddof=1)
    else:
        Cl3a_abund_fake=0.0
        Cl3a_abund=[0.0,0.0]
        
    index_line=line_names.index("Cl3_5538s")
    if lines_available[index_line]=="yes" and Teinput!=0 and Neinput!=0:
        Cl3b_abund_fake = Cl3.getIonAbundance(int_ratio=obs_fake.getIntens(0)['Cl3_5538A'],tem=Teinput_fake, den=Neinput_fake, to_eval='L(5538)', Hbeta=100.0)
        Cl3b_abund = Cl3.getIonAbundance(int_ratio=obs.getIntens(0)['Cl3_5538A'],tem=Teinput, den=Neinput, to_eval='L(5538)', Hbeta=100.0)
        index2=1
        
        Cl3b_abund_fake2=np.nan_to_num(Cl3b_abund_fake)
        Cl3b_abund_fake2[Cl3b_abund_fake2 == inf] =  0

        Cl3b_abund_fake[Cl3b_abund_fake == inf] =  0
        Cl3b_abund_fake=np.nan_to_num(Cl3b_abund_fake)
        Cl3b_abund_fake = [i2 for i2 in Cl3b_abund_fake if i2 != 0]
        Cl3b_abund_fake = [i2 for i2 in Cl3b_abund_fake if i2 < Cl3b_abund*3 and i2>Cl3b_abund*0.3]
        Cl3b_abund_fake=np.std(np.nan_to_num(Cl3b_abund_fake), ddof=1)
    else:
        Cl3b_abund_fake=0.0
        Cl3b_abund=[0.0,0.0]
 
    #Cl3a_abund_fake1=np.nan_to_num(Cl3a_abund_fake)
    #Cl3b_abund_fake2=np.nan_to_num(Cl3b_abund_fake)
  
    #Cl3a_abund_fake=np.std(np.nan_to_num(Cl3a_abund_fake), ddof=1)
    #Cl3b_abund_fake=np.std(np.nan_to_num(Cl3b_abund_fake), ddof=1)
 
    Cl3_abund_fake=0
    Cl3_abund_fake2=0
    Cl3_abund=0
    if (index1==1 and index2==0):
        Cl3_abund_fake=Cl3a_abund_fake
        Cl3_abund_fake2=Cl3a_abund_fake1
        Cl3_abund=Cl3a_abund
    elif (index2==0 and index2==1):
        Cl3_abund_fake=Cl3b_abund_fake
        Cl3_abund_fake2=Cl3b_abund_fake2
        Cl3_abund=Cl3b_abund
    elif (index1==1 and index2==1):
        Cl3_abund_fake=np.sqrt(Cl3a_abund_fake**2+Cl3b_abund_fake**2)/2.
        Cl3_abund_fake2=(Cl3a_abund_fake1+Cl3b_abund_fake2)/2.
        Cl3_abund=(Cl3a_abund+Cl3b_abund)/2.

    if (Cl3a_abund[0]<0.01 and Cl3b_abund[0]<0.01):
        if(Cl3a_abund[0]>=0  and Cl3b_abund[0]>=0 ):
            ion_abun.ClIIIa[i,j]=Cl3a_abund[0]
            ion_abun.ClIIIb[i,j]=Cl3b_abund[0]
            if np.isnan(Cl3a_abund_fake)==True: Cl3a_abund_fake=0
            ion_abun_error.ClIIIa[i]=Cl3a_abund_fake
            if np.isnan(Cl3b_abund_fake)==True: Cl3b_abund_fake=0
            ion_abun_error.ClIIIb[i]=Cl3b_abund_fake
    
#    print >>file10, 'Cl+2(5517)/H+ = {:5.2e}'.format(ion_abun.ClIIIa[i,j]),'Cl+2(5538)/H+ = {:5.2e}'.format(ion_abun.ClIIIb[i,j])
#    print >>file10, 'Cl+2(5517)/H+ = {:5.2e}'.format(ion_abun_error.ClIIIa[i,j]),'Cl+2(5538)/H+ = {:5.2e}'.format(ion_abun_error.ClIIIb[i,j])
#############################################################  
    

    atom_abun_fake = {'He2': He1_abund_fake2, 'He3':He2_abund_fake2,'O2': O2_abund_fake2, 'O3': O3_abund_fake2, 'N2':N2_abund_fake2,'S2': S2_abund_fake2, 'S3': S3_abund_fake2,'Ne3': Ne3_abund_fake2,'Ar3': Ar3_abund_fake2,'Ar4': Ar4_abund_fake2,'Cl3': Cl3_abund_fake2}

    file10.close()
    return ion_abun,ion_abun_error,atom_abun_fake

        
        
