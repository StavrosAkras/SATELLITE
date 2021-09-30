# testing_unequal_arrays_script.py:
# Determines the lenght of 2 arrays 
# (C) Stavros Akras


def testing_unequal_arrays(test,test_error,test_radial,pixscale,slitlength):
    
    for j in range(0,slitlength):
        for i in range(0,len(test_radial)):
            if test_radial[i]!=i*pixscale:
                test_radial[i:i] = [i*pixscale] 
                test[i:i] = [0] 
                test_error[i:i] = [0]            

    
    return test,test_error,test_radial





