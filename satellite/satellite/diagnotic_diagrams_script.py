# diagnotic_diagrams_script.py:
# 
# (C) Stavros Akras

import numpy as np

from satellite import line_ratio_script as lrs
from satellite import positive_values_script as pvs
from satellite import mean_value_script as mvs
from satellite import std_value_script as svs

### NOT YET finished########################################3
        
def NIISIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    
    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,1)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,1)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log(Ha/[N II] 6548+6584)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log(Ha/[S II] 6716+6731)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos


def OIIINIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([N II] 6584/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos


def OIIISIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([S II] 6716+6731/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos
        

def OIIIHeII4686plot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log(He II 4686/Hb)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos
    

def OIIIHeII5412plot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log(He II 5412/Hb)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos



def OIIINIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([N I] 5199/Hb)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos

    
def OIIIOIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([O I] 6300/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos


    
    
def OIIINII5755plot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([N II] 5755/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos


    
    
def SII67166731HaSIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,2)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("[S II] 6716/6731","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([Ha/[S II] 6716+6731)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos


      
    
def SII67166731HaNIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,2)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("[S II] 6716/6731","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([Ha/[N II] 6548+6584)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos



    
    
def NIISIISIIOIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([N II] 6548+6584/[S II] 6716+6731)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([S II] 6716+6731/[O I] 6300+6363)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos


    
    
def NIISIINIIOIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([N II] 6548+6584/[S II] 6716+6731)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([N II] 658+6584/[O I] 6300+6363)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos

    
    
def OIIINII6584OIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([N II] 6584/[O I] 6300)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos



def OIIIOII732030OIIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([O II] 7320+7330/[O III] 5007)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos


    

def OIIIOII372729OIIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([O II] 3727+3729/[O III] 5007)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OIIIOII372729Hbplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([O II] 3727+3729/Hb)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos



def OIIIOII732030Haplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([O II] 7320+7330/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos





def OIIIHeIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log(He I 5876/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OIIIArIIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([Ar III] 7136/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OIIIOIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([O III] 5007/[O I] 6300)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OIIININIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([N I] 5199/[N II] 6584)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OIIINeIIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([Ne III] 3869/Hb)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OIIIOIIIHgplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    
    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([O III] 4363/Hg)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos



def OIIIArIVplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,2)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"[Ar IV] 4712/4740","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OIIICIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([C I] 8727/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OIIICIIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O III] 5007/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([C II] 6461/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos




def OICIplot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([O I] 6300/Ha)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log([C I] 8727/Ha)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos





def ArIVHeII4686plot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([Ar IV] 4712+4740/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log(He II 4686/Hb)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2, file=file10)

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos



def ArIVHeII5412plot(numerator1,denominator1,numerator2,denominator2,sizex,sizey):
    file10 = open('general_output_file.txt','a')

    data_pos_ratio1=[]

    data_ratio1=lrs.lineratio(numerator1,denominator1,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaNII.fits', dataNIIHa, hdr)
    data_pos_ratio1=pvs.positivevalues(data_ratio1,numerator1,denominator1,sizex,sizey,0)
    mean_ratio1=mvs.meanvalue(data_pos_ratio1,-10) 
    std_ratio1=svs.stdvalue(data_pos_ratio1,mean_ratio1,-10)

    data_pos_ratio2=[]
    #fits.writeto('results/test.fits', datasiiab/maxHb, hdr)
    data_ratio2= lrs.lineratio(numerator2,denominator2,sizex,sizey,0)
    #fits.writeto('results_models/abell14_logHaSII.fits', dataSIIHa, hdr)
    data_pos_ratio2=pvs.positivevalues(data_ratio2,numerator2,denominator2,sizex,sizey,0)
    mean_ratio2=mvs.meanvalue(data_pos_ratio2,-10) 
    std_ratio2=svs.stdvalue(data_pos_ratio2,mean_ratio2,-10) 

    print("log([Ar IV] 4712+4740/Hb)","%.3f" % mean_ratio1,"sigma","%.3f" % std_ratio1,"log(He II 5412/Hb)","%.3f" % mean_ratio2,"sigma","%.3f" % std_ratio2 >>file10, end=' ', file=file10) 

    maxper=75
    minper=25
    data_pos_ratio1 = [x for x in data_pos_ratio1 if x != -100000]
    data_pos_ratio2 = [x for x in data_pos_ratio2 if x != -100000]

#    per25n=np.percentile(data_pos_ratio1,minper)
#    per75n=np.percentile(data_pos_ratio1,maxper)
#    per25s=np.percentile(data_pos_ratio2,minper)
#    per75s=np.percentile(data_pos_ratio2,maxper)

    datax=[]
    datay=[]
    k=0
    for j in range(0,sizex):
        for i in range(0,sizey):
        
            datax.append(-99)
            datay.append(-99)

            if (numerator1[i,j] ==0 or denominator1[i,j] == 0):     
                datay[k]=-99  
                
            if denominator1[i,j] < 0: 
                datay[k]=-99        
        
            if numerator1[i,j] > 0 and denominator1[i,j] >0:             
                datay[k]=data_ratio1[i,j]

            if (numerator2[i,j] == 0  or denominator2[i,j] == 0): 
                datax[k]=-99
            
            if denominator2[i,j] <0: 
                datax[k]=-99
            
            if numerator2[i,j] > 0 and denominator2[i,j]  > 0:
                datax[k]=data_ratio2[i,j]

            k=k+1

    xpos=[]
    ypos=[]
    l=0 
    for j in range(0,len(datax)):
        if (datay[j] > -98 and datax[j] > -98): 
            xpos.append(-99)
            ypos.append(-99)
            xpos[l]=datax[j]
            ypos[l]=datay[j]
            l=l+1
#    print(per25n,per75n,per25s,per75s)
#    datax,datay,dataz,datad=xyzd(numerator1,denominator1,numerator2,denominator1,data_ratio1,data_ratio2,sizex,sizey,per25n,per75n,per25s,per75s)
#    xx1pos,yy1pos=diagnostic(datax,datay,dataz,datad,1,4) #idex=4 no BPT
    
    
    return xpos,ypos
