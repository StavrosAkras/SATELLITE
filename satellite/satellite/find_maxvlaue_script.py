# find_maxvlaue_script.py:
# Finds and returns the maximum value from 1D table.
# The maximum flux value can be find within the whole slit (1D array) or 
# in a specific region based on the rlimit parameter 
# (e.g. for a specific distance from the centrar stars larger than rlimit)
# (C) Stavros Akras



def maxvalue(flux,dis,rlimnname,rlimitvalue,pixsc):

        maxflux=0
        maxrad=0
        index_rlimit=rlimnname.index("limit_radial_in_arcsec")
        
        rlimit=float(rlimitvalue[index_rlimit])

        for i in range(0,len(flux)):
            if dis[i] > rlimit and flux[i]>maxflux:
                maxflux=flux[i]
                maxrad=dis[i]
                
        return maxflux,maxrad
