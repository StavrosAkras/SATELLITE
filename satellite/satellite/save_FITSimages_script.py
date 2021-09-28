# save_FITSimages_script.py:
# Saves the 2D arrays as FITS image. The name is taken from the input.txt file and 
# only if it is requested by the user (outputs.txt file defined as "yes").
# The routine also searches in the current folder if the image is already saved.
# (C) Stavros Akras


import numpy as np
from astropy.io import fits
import os, fnmatch



def save_image(a1,a2,name,parameter,request):
               
    index_req=parameter.index(name)
    if request[index_req]=="yes":
        file_save_exist=0
        listOfFiles = os.listdir('output_images')
        pattern = name+".fits"
        for entry in listOfFiles:
            if fnmatch.fnmatch(entry, pattern):
                file_save_exist=1
                
        a11=np.nan_to_num(a1)
        index_par=parameter.index(name)
        if parameter[index_par]==name:
            if file_save_exist==1:
                os.remove("output_images/"+name+".fits")
                fits.writeto("output_images/"+name+".fits", a11, a2)
            elif file_save_exist==0:
                fits.writeto("output_images/"+name+".fits", a11, a2)        
        
    return 

def save_image_ionabun(a1,a2,name1,name,parameter,request,line_names,lines_available):
               
    index_req=line_names.index(name1)
    if lines_available[index_req]=="yes":
        file_save_exist=0
        listOfFiles = os.listdir('output_images')
        pattern = name+".fits"
        for entry in listOfFiles:
            if fnmatch.fnmatch(entry, pattern):
                file_save_exist=1

        a11=np.nan_to_num(a1)
        index_par=line_names.index(name1)
        if lines_available[index_par]=="yes":            
            if file_save_exist==1:
                os.remove("output_images/"+name+".fits")
                fits.writeto("output_images/"+name+".fits", a11, a2)
            elif file_save_exist==0:
                fits.writeto("output_images/"+name+".fits", a11, a2)        
        
    return 

def save_image2(a1,a2,name1,name,parameter,request):
               
    index_req=parameter.index(name1)
    if request[index_req]=="yes":
        file_save_exist=0
        listOfFiles = os.listdir('output_images')
        pattern = name+".fits"
        for entry in listOfFiles:
            if fnmatch.fnmatch(entry, pattern):
                file_save_exist=1
        
        a11=np.nan_to_num(a1)
        index_par=parameter.index(name1)
        if parameter[index_par]==name1:
            if file_save_exist==1:
                os.remove("output_images/"+name+".fits")
                fits.writeto("output_images/"+name+".fits", a11, a2)
            elif file_save_exist==0:
                fits.writeto("output_images/"+name+".fits", a11, a2)        
        
    return 
