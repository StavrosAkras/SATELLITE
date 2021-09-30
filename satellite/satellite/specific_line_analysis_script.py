# specific_line_analysis_script.py:
# Calculates the fluxes for the 10 specific pseudo-slits. The "rotate" routine first rotates the images clockwise.
# Note: The images has its 0,0 pixel at the upper-left corner and increase from left to right 
# and from top to bottom. This means that it is flipped compared to the ds9 representation. 
# The fits images has its 0,0 pixel at the lower-left corner and it increases from left 
# to right and from bottom to top.
# (C) Stavros Akras

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
import matplotlib.patches as patches
import math
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

    
def analysis_specific_line(flux,flux_err,angle,haa,hbb,widthslit,lengthslit,x_CS_old,y_CS_old,x_slit_old,y_slit_old,param_estimated,param_requered,num_slit):
       
        flux_rot=rotate(flux,angle,order=5)
        flux_err_rot=rotate(flux_err,angle,order=5)
        
        maxvalue=np.amax(flux)
        minvalue=maxvalue/1000
        pixels_conversion=1 

        ha_rot=rotate(haa,angle,order=5)
        hb_rot=rotate(hbb,angle,order=5)
        size_rotx=flux_rot.shape[0]
        size_roty=flux_rot.shape[1]    

        xoffset=x_slit_old-x_CS_old
        yoffset=y_slit_old-y_CS_old

        centerx=(size_rotx)/2+(xoffset*math.cos((-1)*math.pi*angle/180)-yoffset*math.sin((-1)*math.pi*angle/180))
        centery=(size_roty)/2+(xoffset*math.sin((-1)*math.pi*angle/180)+yoffset*math.cos((-1)*math.pi*angle/180))

################################################################################################
### This part is for the slit position testing module which is used to test the position of the 
### slit on the rotated image.
################################################################################################

        index=param_estimated.index('slit_position_testing')
        if param_requered[index]=="yes":
            
                fig= plt.figure(figsize=(10,10))
                ax=fig.add_subplot(111)
                plt.imshow(flux_rot,norm=LogNorm(minvalue,maxvalue))
                rect1 = patches.Rectangle((centerx-1-(widthslit/2),centery-(lengthslit/2)),(widthslit/pixels_conversion)-1,lengthslit/pixels_conversion,linewidth=1,edgecolor='k',facecolor='none')
                ax.add_patch(rect1)
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.45)
                plt.colorbar(cax=cax)
                plt.savefig("output_images/fig_slit"+str(num_slit)+".png")
                plt.savefig("output_images/fig_slit"+str(num_slit)+".pdf")
#               plt.show()

        centerx=int(centerx)
        centery=int(centery)
        widthslit_new=int(widthslit/pixels_conversion)
        lengthslit_new=int(lengthslit/pixels_conversion)
        
        #x,y must be the opposite!
        centery_new=centerx
        centerx_new=int(centery-(lengthslit/2))

        
        k=0
        sum=sumerror=0.0  
        for i in range(centerx_new,centerx_new+lengthslit_new):
            if (widthslit_new % 2)==0:
                for j in range(centery_new-widthslit_new//2,centery_new+widthslit_new//2): 
 
                    if (flux_rot[i,j] ==0):
                        sum=sum  
                    if flux_rot[i,j] > 0 and ha_rot[i,j] >= hb_rot[i,j]*2.86  and hb_rot[i,j] > 0:
                        sum=sum+flux_rot[i,j]
                        sumerror=sumerror+flux_err_rot[i,j]**2
                        k=k+1       
            else:
                for j in range(centery_new-(widthslit_new-1)//2,centery_new+(widthslit_new-1)//2+1):
                    if (flux_rot[i,j] ==0):
                        sum=sum  
                    if flux_rot[i,j] > 0 and ha_rot[i,j] >= hb_rot[i,j]*2.86  and hb_rot[i,j] > 0:
                        sum=sum+flux_rot[i,j]
                        sumerror=sumerror+flux_err_rot[i,j]**2
                        k=k+1  
                
        return sum,np.sqrt(sumerror),k      
                
    
