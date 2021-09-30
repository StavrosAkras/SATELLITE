# positive_values_script.py:
# Finds if the nominator and the denomitor of a ratio 
# (e.g. line ratio, intensities) from two 2D maps are both positive 
# and save the values in an 1D array, otherwise put a value of -1000000
# (C) Stavros Akras



def positivevalues(a1,nom1,dom1,sx,sy,idex):
        l1=0
        datapos=[]
        for j in range(0,sx):
            for i in range(0,sy):
                datapos.append(idex)
  
                if nom1[i,j] > 0 and dom1[i,j] >0:
                    datapos[l1]=a1[i,j]
                    l1=l1+1
                    
                if nom1[i,j] <= 0 or dom1[i,j] <= 0:
                    datapos[l1]=-100000
                    l1=l1+1
                     
                 
        return datapos
