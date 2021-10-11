# slit_line_flux_script.py:
# Calculates the flux of an emission line along a pseudo-slit at a specific position angle, width and length given by the user (numerical_input.txt). 
# The width should always be an integer number. The pseudo-slit starts from the center of the image/center of the nebula, and it covers only one half of the nebula. 
# For flux calculation, the code sums up the values  of all the pixels except those which have F($_Hb$)$<$0 and/or F($_Ha$)$<$F($_Hb$)*2.86 (negative or unrealistic c(Hb)).
# Morevoer, the routine rotates the entire image/table by a given angle, then calculates the new size (x,y) of the rotated imaged and the center of the new image. 
# The pseudo-slit is always considered along the up-down direction. The routine returns the total flux and the number of spaxels used for the calculation

# Note: It is necessary the image be large enough so after the rotation, the entire nebula or galaxy remains in the new image. Sometimes an elongated nebula 
# (or even a galaxy) is observed in a specific PA, so entire nebula fits the field of view of the instrument (rotation angle=0 corresponds to observed PA=0) 
# (the orientation of the images is north up and east to the left) otherwise the user has to take into account the offset between the north direction and 
# the position of the sources.
# (C) Stavros Akras

import numpy as np
from scipy.ndimage.interpolation import rotate

def slitfluxline(flux,flux_err,angle,haa,hbb,widthslit,lengthslit,max_width):
       
        flux_rot=rotate(flux,angle,order=5)
        flux_err_rot=rotate(flux_err,angle,order=5)

        pixels_conversion=1 

        ha_rot=rotate(haa,angle,order=5)
        hb_rot=rotate(hbb,angle,order=5)
        size_rotx=flux_rot.shape[0]
        size_roty=flux_rot.shape[1]    
        centerx=(size_rotx)//2
        centery=(size_roty)//2
        
        widthslit_new=int(widthslit/pixels_conversion)
        lengthslit_new=int(lengthslit/pixels_conversion)
       
        ###for slits 
        if widthslit<max_width and lengthslit<=max_width:
            k=0
            sum=sumerror=0.0
            for i in range(centerx,centerx+lengthslit_new): # i increases vertially from bottom 
                if (widthslit_new % 2)==0:
                    for j in range(centery-widthslit_new//2,centery+widthslit_new//2): # j increases horizontally from left
                    ### Fluxes and errors
                        if (flux_rot[i,j] ==0):
                            sum=sum  
                        if flux_rot[i,j] > 0 and ha_rot[i,j] >= hb_rot[i,j]*2.86  and hb_rot[i,j] > 0:
                            sum=sum+flux_rot[i,j]
                            sumerror=sumerror+flux_err_rot[i,j]**2
                            k=k+1    

                            
                else:
                    for j in range(centery-(widthslit_new-1)//2,centery+(widthslit_new-1)//2+1): 
                        if (flux_rot[i,j] ==0):
                            sum=sum  
                        if flux_rot[i,j] > 0 and ha_rot[i,j] >= hb_rot[i,j]*2.86  and hb_rot[i,j] > 0:
                            sum=sum+flux_rot[i,j]
                            sumerror=sumerror+flux_err_rot[i,j]**2
                            k=k+1 
                            
        file10 = open('general_output_file.txt','a')
        if widthslit>max_width or lengthslit>max_width: print("Sorry, your slit width or/and length are larger that the true size of the image", file=file10)
        
         ###for the entire nebula
        if widthslit==max_width and lengthslit==max_width:
            k=0
            sum=0.0
            for i in range(centerx-lengthslit//2,centerx+lengthslit//2): # i increases vertically from bottom 
                for j in range(centery-widthslit//2,centery+widthslit//2): # j increases horizontally from left
                    if (flux_rot[i,j] ==0):
                        sum=sum  
                    if flux_rot[i,j] > 0 and ha_rot[i,j] >= hb_rot[i,j]*2.86  and hb_rot[i,j] > 0:
                        sum=sum+flux_rot[i,j]
                        sumerror=sumerror+flux_err_rot[i,j]**2
                        k=k+1    

        return sum,np.sqrt(sumerror),k

