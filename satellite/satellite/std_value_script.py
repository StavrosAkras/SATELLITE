# std_value_script.py:
# Computes and returns the standard deviation of a 1D array (e.g. emission line ratio table).
# The "value" parameter defines the minimum value of the table than will be considered in the 
# estimation of the standard deviation, excluding the pixels wth a value of -100000. 
# The mean has been calculated in the mean_value_script.
# (C) Stavros Akras

import numpy as np


def stdvalue(a1,mean,value):
        sum=0.0
        l=0
        for j in range(0,len(a1)):
                if a1[j] > value:        
                    sum=sum+(a1[j]-mean)*(a1[j]-mean)
                    l=l+1
                    
        return np.sqrt(sum/(l-1))
