# radial_line_flux_script.py:
# Computes the flux of an emission line based on its 2D map within a pseudo-slit of a specific position angle, 
# width and length given by the user. The "rotate" routine rotates the 2D images (or 2D array) clockwise. 
# Note: The images has its 0,0 pixel at the upper-left corner and increase from left to right 
# and from top to bottom. This means that it is flipped compared to the ds9 representation. 
# The fits images has its 0,0 pixel at the lower-left corner and it increases from left to right and from bottom to top.
# The fluxes, error and radial distances are return in three 1D arrays.
# (C) Stavros Akras


import numpy as np
from scipy.ndimage.interpolation import rotate

###
### This script/function 


def radialfluxline(flux,flux_err,angle,haa,hbb,widthslit,lengthslit,max_width,pix_sc):
                
       
        flux_rot=rotate(flux,angle,order=5)
        flux_err_rot=rotate(flux_err,angle,order=5)

        pixels_conversion=1 

        ha_rot=rotate(haa,angle,order=5)
        hb_rot=rotate(hbb,angle,order=5)
        size_rotx=flux_rot.shape[0]
        size_roty=flux_rot.shape[1]    
        centerx=(size_rotx)/2
        centery=(size_roty)/2
        
        
        widthslit_new=int(widthslit/pixels_conversion)
        lengthslit_new=int(lengthslit/pixels_conversion)
        

        radial_flux_table=np.zeros((centerx+lengthslit_new))
        radial_flux_table_error=np.zeros((centerx+lengthslit_new))

        arcsec=np.zeros((centerx+lengthslit_new))
        pixel_scale=pix_sc
        
        ###for slits 
        if widthslit<max_width and lengthslit<=max_width:
            
            k=0
            for i in range(centerx,centerx+lengthslit_new): # i anebainei pros ta panw
                sum=sumerror=0.0
                if (widthslit_new % 2)==0:
                    for j in range(centery-widthslit_new/2,centery+widthslit_new/2): # to j einai orizontia
                        if (flux_rot[i,j] ==0):
                            sum=sum
                        if flux_rot[i,j] > 0 and ha_rot[i,j] >= hb_rot[i,j]*2.86  and hb_rot[i,j] > 0:
                            sum=sum+flux_rot[i,j]
                            sumerror=sumerror+flux_err_rot[i,j]**2
                        else:
                            sum=sum+0
                            sumerror=sumerror+0
                            
                    radial_flux_table[k]=sum
                    radial_flux_table_error[k]=np.sqrt(sumerror)
                    arcsec[k]=k*pixel_scale
                    k=k+1
                    
                else:
                    for j in range(centery-(widthslit_new-1)/2,centery+(widthslit_new-1)/2+1): # to j einai orizontia
                        if (flux_rot[i,j] ==0):
                            sum=sum
                        if flux_rot[i,j] > 0 and ha_rot[i,j] >= hb_rot[i,j]*2.86  and hb_rot[i,j] > 0:
                            sum=sum+flux_rot[i,j]
                            sumerror=sumerror+flux_err_rot[i,j]**2
                        else:
                            sum=sum+0
                            sumerror=sumerror+0    

                    radial_flux_table[k]=sum
                    radial_flux_table_error[k]=np.sqrt(sumerror)
                    arcsec[k]=k*pixel_scale
                    k=k+1
        
        ll=0
        radial_flux_table2,radial_flux_table2_error,arcsec2=[],[],[]
        for i in range(0,len(radial_flux_table)):
            if radial_flux_table[i]!=0: 
                radial_flux_table2.append(0)
                radial_flux_table2_error.append(0)
                arcsec2.append(0)
                radial_flux_table2[ll]=radial_flux_table[i]
                radial_flux_table2_error[ll]=radial_flux_table_error[i]
                arcsec2[ll]=arcsec[i]
                ll=ll+1
        
        file10 = open('general_output_file.txt','a')        
        if widthslit>max_width or lengthslit>max_width: print("Sorry, your slit width or/and length are larger that the true size of the image", file=file10)
        
        return radial_flux_table2,radial_flux_table2_error,arcsec2       
                
                


