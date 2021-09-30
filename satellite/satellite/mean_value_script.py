# mean_value_script.py:
# Computes the mean value from a 2D table (e.g. emission line ratio table).
# The value parameter defines the minimum value of the table that should be considered for the 
# estimation of the mean value, excluding the spaxles for which the unrealistic value of -100000 
# had been given in the line_ratio_script.py
# (C) Stavros Akras




def meanvalue(a1,value):
        l1=0
        sum=0.0
        for j in range(0,len(a1)):
            if a1[j]>value:
                sum=sum+a1[j]
                l1=l1+1
                
        if l1==0:
            sum=0.0
            l1=1
        
        return sum/l1

def meanvalue2(a1,value1,value2):
        l1=0
        sum=0.0
        for j in range(0,len(a1)):
            if a1[j]>value1 and a1[j]<value2:
                sum=sum+a1[j]
                l1=l1+1
                
        if l1==0:
            sum=0.0
            l1=1
        
        return sum/l1
