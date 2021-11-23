#!/usr/bin/env python
 
# SATELLITE, the Spectroscopic Analysis Tool for intEgraL fieLd unIt daTacubEs
# (C) 2021- Stavros Akras


# SATELLITE uses several software packages in Python: Matplotlib
# (Hunter 2007), NumPy (van der Walt et al. 2011), SciPy (Virtanen
# et al. 2020), seaborn (Waskom & the seaborn development team
# 2020) and AstroPy (Astropy Collaboration et al. 2013, 2018), 
# as well as the PyNeb package developed by Luridiana V., Morisset C., 
# Shaw R. A., 2015, A&A, 573, A42

# SATELLITE is a free code: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation.

# It is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


import warnings
import datetime

from satellite import read_input_script as read
from satellite import read_input_lines_parameters_script as read_ilps
from satellite import radial_analysis_script as ras
from satellite import rotate_line_fluxes_script as rotlfs
from satellite import specificPA_line_fluxes_script as sPAlfs
from satellite import analysis2D_script as a2Ds
from satellite import my_analysis2D_script as mya2Ds
from satellite import TeNe_angles_script as TeNeangles
from satellite import TeNe_specific_slits_script as TeNeslits
from satellite import TeNe_2D_script as TeNe2Ds
from satellite import TeNe_radial_script as TeNers
from satellite import flux_angles_norm_definition_script as fands
from satellite import define_Te_Ne_for_null_spaxels_script as dTeNenullspaxelss
from satellite import exclude_probematic_errors_script as epes
from satellite import physical_parameters_script as pps

class flux2D:
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

class flux2D_error:
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
 
    
    

global file10

print('=========================================') 
prog_start = datetime.datetime.now()
print("Start Program")
print(prog_start.strftime("%Y-%m-%d %H:%M:%S"))
print('=========================================') 


warnings.filterwarnings("ignore")
file10 = open('general_output_file.txt','w')

#######################################################################################################
# read all input parameters given by the user.
#######################################################################################################
line_names,line_ext_error,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,DD_name,DD_avail,DDxmin,DDxmax,DDymin,DDymax,par_plotname,par_plotymin,par_plotymax=read_ilps.read_input_lines_parameters()

### if new lines must be included in the code, they have to be added in the input file, read them in the read.read_input_script for the flux and error 
### and finally define their tables in the flux classes

#######################################################################################################
# read line flux images  
#######################################################################################################
flux2D.Ha_6563,flux2D.Hb_4861,flux2D.Hg_4340,flux2D.Hd_4101,flux2D.HeIa_5876,flux2D.HeIb_6678,flux2D.HeIIa_4686,flux2D.HeIIb_5412,flux2D.NIIa_5755,flux2D.NIIb_6548,flux2D.NIIc_6584,flux2D.NI_5199,flux2D.OIIIa_4363,flux2D.OIIIb_4959,flux2D.OIIIc_5007,flux2D.OIIa_3727,flux2D.OIIb_3729,flux2D.OIIc_7320,flux2D.OIId_7330,flux2D.OIa_5577,flux2D.OIb_6300,flux2D.OIc_6363,flux2D.SIIa_6716,flux2D.SIIb_6731,flux2D.SIIIa_6312,flux2D.SIIIb_9069,flux2D.ClIIIa_5517,flux2D.ClIIIb_5538,flux2D.ArIII_7136,flux2D.ArIVa_4712,flux2D.ArIVb_4740,flux2D.CI_8727,flux2D.CII_6461,flux2D.NeIIIa_3868,flux2D.NeIIIb_3967,hdr=read.read_input_images(line_names,lines_available,param_model_values,"fluxes")


#######################################################################################################
# read line flux error images  
#######################################################################################################
flux2D_error.Ha_6563,flux2D_error.Hb_4861,flux2D_error.Hg_4340,flux2D_error.Hd_4101,flux2D_error.HeIa_5876,flux2D_error.HeIb_6678,flux2D_error.HeIIa_4686,flux2D_error.HeIIb_5412,flux2D_error.NIIa_5755,flux2D_error.NIIb_6548,flux2D_error.NIIc_6584,flux2D_error.NI_5199,flux2D_error.OIIIa_4363,flux2D_error.OIIIb_4959,flux2D_error.OIIIc_5007,flux2D_error.OIIa_3727,flux2D_error.OIIb_3729,flux2D_error.OIIc_7320,flux2D_error.OIId_7330,flux2D_error.OIa_5577,flux2D_error.OIb_6300,flux2D_error.OIc_6363,flux2D_error.SIIa_6716,flux2D_error.SIIb_6731,flux2D_error.SIIIa_6312,flux2D_error.SIIIb_9069,flux2D_error.ClIIIa_5517,flux2D_error.ClIIIb_5538,flux2D_error.ArIII_7136,flux2D_error.ArIVa_4712,flux2D_error.ArIVb_4740,flux2D_error.CI_8727,flux2D_error.CII_6461,flux2D_error.NeIIIa_3868,flux2D_error.NeIIIb_3967,hdrer=read.read_input_images(line_names,lines_available,param_model_values,"errors")



flux2D,flux2D_error=epes.exclude_problematic_error_values(flux2D,flux2D_error,param_mod_name,param_model_values)  

        
index_size=param_mod_name.index("total_num_pixels_horiz")
index_pixel_scale=param_mod_name.index("pixel_scale")
sizex=param_model_values[index_size]
sizey=param_model_values[index_size]
pixscale=float(param_model_values[index_pixel_scale])*0.01
minx=miny=0
print('x=',sizex,'y=', sizey,'pixel scale=',pixscale, file=file10)   

#######################################################################################################
# calculate line fluxes, normalized lines relative to Hbeta, and a number of line ratios for radial slits with PA from X1 to X2 given by the user.
#######################################################################################################
index=param_estimated.index('angular_analysis_task')
if param_requered[index]=="yes":
    flux_angles,flux_angles_error,flux_angles_norm,flux_angles_norm_error,ratio_angles,ratio_angles_error,angles=rotlfs.rotate_line_fluxes(flux2D,flux2D_error,line_names,line_ext_error,lines_available,param_estimated,param_requered,param_mod_name,param_model_values,par_plotname,par_plotymin,par_plotymax)
    
#######################################################################################################
# use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as 
# Te, Ne, ionic and chemical abudnances,  and ICFs
#######################################################################################################    
    Te_PA, Ne_PA = TeNeangles.TeNe(flux_angles,flux_angles_error,angles,line_names,lines_available,param_estimated,param_requered,param_mod_name,param_model_values,par_plotname,par_plotymin,par_plotymax)



#######################################################################################################
# radial analysis of various emission lines for a specific PA of a slit 
#######################################################################################################
index=param_estimated.index('radial_analysis_task')
if param_requered[index]=="yes":
    flux_radial,flux_radial_error,distance_CS=ras.radial_analysis(flux2D,flux2D_error,line_names,line_ext_error,lines_available,lines_radial,param_mod_name,param_model_values,par_plotname,par_plotymin,par_plotymax)

#######################################################################################################
# use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as 
# Te, Ne, ionic and chemical abudnances,  and ICFs
#######################################################################################################    
    Te_rad, Ne_rad = TeNers.TeNe(flux_radial,flux_radial_error,distance_CS,line_names,lines_available,param_estimated,param_requered,param_mod_name,param_model_values,par_plotname,par_plotymin,par_plotymax)


#######################################################################################################
# calculate line fluxes, normalized lines relative to Hbeta, and a number of line ratios for specific slits with PAs, widths and lengths given by the user.
#######################################################################################################
index=param_estimated.index('specific_slits_analysis_task')
if param_requered[index]=="yes":
    flux_spec_slit,flux_spec_slit_error,flux_spec_slit_norm,flux_spec_slit_norm_error,ratio_spec_slit,ratio_spec_slit_error=sPAlfs.specficPA_line_fluxes(flux2D,flux2D_error,line_names,line_ext_error,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values)

#######################################################################################################
# use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as 
# Te, Ne, ionic and chemical abudnances,  and ICFs
#######################################################################################################
    Te_slits, Ne_slits = TeNeslits.TeNe(flux_spec_slit,flux_spec_slit_error,line_names,lines_available,param_estimated,param_requered,param_mod_name,param_model_values)


#######################################################################################################
# 2D analysis and diagnostic diagrams
#######################################################################################################
index=param_estimated.index('2D_anaysis_task')
if param_requered[index]=="yes":
    
    index=param_estimated.index('angular_analysis_task')
    if param_requered[index]=="no":
        flux_angles_norm,flux_angles_norm_error,angles=fands.flux_ang_norm_definition(param_mod_name,param_model_values)
        Te_PA,Ne_PA=dTeNenullspaxelss.nullTeNe(param_mod_name,param_model_values)
   
    index=param_estimated.index('specific_slits_analysis_task')
    if param_requered[index]=="no":
        flux_spec_slit_norm,flux_spec_slit_norm_error=fands.flux_ang_norm_definition_specificslits(param_mod_name,param_model_values)
        Te_slits,Ne_slits=dTeNenullspaxelss.nullTeNe_specific_slits(param_mod_name,param_model_values)

    a2Ds.analysis2D(flux2D,flux2D_error,flux_angles_norm,angles,line_names,line_ext_error,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,DD_name,DD_avail,DDxmin,DDxmax,DDymin,DDymax,hdr,flux_spec_slit_norm,flux_spec_slit_norm_error)

#######################################################################################################
# use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as 
# Te, Ne, ionic and chemical abudnances,  and ICFs
#######################################################################################################
    Te_2D, Ne_2D = TeNe2Ds.TeNe(flux2D,flux2D_error,line_names,line_ext_error,lines_available,param_mod_name,param_model_values,param_estimated,param_requered,hdr,Te_PA,Ne_PA)

#######################################################################################################
# user can write his/her proper script for a 2D analysis and diagnostic diagrams
#######################################################################################################
index=param_estimated.index('My_2D_anaysis_task')
if param_requered[index]=="yes":
    mya2Ds.analysis2D(flux2D,flux_angles_norm,angles,param_mod_name,param_model_values)

#######################################################################################################
# use the derived fluxes and ratios as input to Pyneb and calculate paramteres such as 
# Te, Ne, ionic and chemical abudnances,  and ICFs
#######################################################################################################
#pps.analysiswithPyNeb(flux2D,flux2D_error,flux_angles,flux_angles_error,flux_spec_slit,flux_spec_slit_error,angles,line_names,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,hdr)

#######################################################################################################
# user can write his/her proper PyNeb script
#######################################################################################################
index=param_estimated.index('My_physical_parameters_task')
if param_requered[index]=="yes":
    pps.analysiswithPyNeb(flux2D,flux_angles,flux_spec_slit,angles,line_names,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,hdr)

file10.close()

prog_end = datetime.datetime.now()
print('=========================================') 
print("End Program")
print(prog_end.strftime("%Y-%m-%d %H:%M:%S"))
print("Elapsed time")
print((prog_end - prog_start))
print('=========================================') 

