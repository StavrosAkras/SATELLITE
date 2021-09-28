# line_ratio_script.py:
# Compute the line ratio or logarithmic line ratio between a and b from the 2D maps. 
# The idex parameters defines which line is the nominator and which the denominator
# idex=0: the first line is the nominator and the second the denominator "the logarithmic ratio"
# idex=1: the first line is the denominar and the second the nominator "the logarithmic ratio" (inverse of idex=1)
# idex=2: the first line is the nominator and the second the denominator NOT logarithmic, simple ratio
# idex=3: the first line is the denominar and the second the nominator NOT logarithmic, simple ratio
# in case of negative spaxel values, the value -100000  (for the logarithmic ratio) and -0.000000001 (ratio) is given
# which are unrealistic and will be excluded.
# (C) Stavros Akras

import numpy as np



def lineratio(nom1,dom1,sx,sy,idex):
        a3=np.zeros((sy, sx))
       
        for j in range(0,sx):
            for i in range(0,sy):
                a3[i,j] =0 
                
                if idex==0:
                    if nom1[i,j] > 0 and dom1[i,j] >0:
                        a3[i,j]=np.log10(nom1[i,j]/dom1[i,j])
            
                    if nom1[i,j] <= 0 or dom1[i,j] <= 0:
                        a3[i,j]=-100000.0
                
                if idex==1:
                    if nom1[i,j] > 0 and dom1[i,j] >0:
                        a3[i,j]=np.log10(dom1[i,j]/nom1[i,j])
            
                    if nom1[i,j] <= 0 or dom1[i,j] <= 0:
                        a3[i,j]=-100000.0
 
                if idex==2:
                    if nom1[i,j] > 0 and dom1[i,j] >0:
                        a3[i,j]=nom1[i,j]/dom1[i,j]
            
                    if nom1[i,j] <= 0 or dom1[i,j] <= 0:
                        a3[i,j]=-0.000000001
                        
                if idex==3:
                    if nom1[i,j] > 0 and dom1[i,j] >0:
                        a3[i,j]=dom1[i,j]/nom1[i,j]
            
                    if nom1[i,j] <= 0 or dom1[i,j] <= 0:
                        a3[i,j]=-0.000000001


        return a3
