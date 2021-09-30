# calculations_excluding_outliers_script.py: 
# Finds the outliers from a sample of points (pixels)
# and computes the new mean and standard deviation.
# (C) Stavros Akras

import numpy as np


##################################################################################
### finds the outliers from a 2D (a2) array. 
### array a1 has only the positive values of the array a2
##################################################################################

def estimates_without_outliers(a1,a2,sx,sy):
    
##################################################################
### The parameters minn and maxx correspond to the percentile in %    
##################################################################    
    minn=5
    maxx=95
    
    per25chb=np.percentile(a1,minn)
    per75chb=np.percentile(a1,maxx)
    iqr=per75chb-per25chb

    if per25chb-1.5*iqr<0: per25chb=1.5*iqr
    
    sumchb=0.0
    k=0
    for j in range(0,sx):
        for i in range(0,sy):
            if (a2[i,j] > per25chb-1.5*iqr and a2[i,j] < per75chb+1.5*iqr):           
                sumchb= sumchb + a2[i,j]
                k=k+1

    meanchb=sumchb/k            
    sumchb=0.0
    k=0
    for j in range(0,sx):
        for i in range(0,sy):
            if (a2[i,j] > per25chb-1.5*iqr and a2[i,j] < per75chb+1.5*iqr):
                sumchb=sumchb+(a2[i,j]-meanchb)*(a2[i,j]-meanchb)
                k=k+1
            
    stdchb=np.sqrt(sumchb/(k-1))
    
    return meanchb,stdchb

##################################################################################
### finds the outliers from an 1D (a1) array
##################################################################################

def estimates_without_outliers1D(a1,number):

##################################################################
### The parameters minn and maxx correspond to the percentile in %    
##################################################################
    minn=5
    maxx=95
    
    a2=a1
    a1=np.absolute(a1)
    per25chb=np.percentile(a1,minn)
    per75chb=np.percentile(a1,maxx)
    iqr=per75chb-per25chb

    if per25chb-1.5*iqr<0: per25chb=1.5*iqr
    
    sumchb=0.0
    k=0
    for j in range(0,number):
        if (abs(a1[j]) > per25chb-1.5*iqr and abs(a1[j]) < per75chb+1.5*iqr):
            if a1[j]!=0:
                sumchb= sumchb + a1[j]
                k=k+1
                    
    if k!=0: 
        meanchb=sumchb/k            
    else:
        meanchb=0.0
        
    sumchb=0.0
    k=0
    for j in range(0,number):
        if (abs(a1[j]) > per25chb-1.5*iqr and abs(a1[j]) < per75chb+1.5*iqr):
            if a1[j]!=0:
                sumchb=sumchb+(a1[j]-meanchb)*(a1[j]-meanchb)
                k=k+1
    
    if k!=0:
        stdchb=np.sqrt(sumchb/(k-1))
    else:
        stdchb=0.0
    
    sumchb=0.0
    for j in range(0,number):
        sumchb= sumchb + a2[j]
        
    if sumchb<0:
        meanchb=(-1)*meanchb
        
    return meanchb,stdchb
