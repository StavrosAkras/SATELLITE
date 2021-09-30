# norm_flux_error_script.py:
# Computes the percentange error for a line ratio a1/a2 (emission line intensity with Hb=100)
# (C) Stavros Akras

import numpy as np

def flux_Error(a1,sa1,a2,sa2):
    
    if a1!=0 and a2!=0:
        par1=sa1/a1
        par2=sa2/a2
        
        error=(a1*100/a2)*np.sqrt(par1**2+par2**2)

    else:
        par1=par2=0
        error=0 

    return error


def ratio_Error(a1,sa1,a2,sa2):
    
    if a1!=0 and a2!=0:
        par1=sa1/a1
        par2=sa2/a2
    else:
        par1=par2=0
        
    error=np.sqrt(par1**2+par2**2)

    return error
