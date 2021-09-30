# generate_2D_lineratio_maps_script.py:
# Defines and computes the emission line ratios which are saves in 2D arrays,
# and then returns fits images for each line ratio.
# The routine also computes and saves in the  statisticnumbers_2.txt ascii file the
# Parameter Name, Nspaxels, 5% value, Q1 value, median value, Q3 value, 95% value, mean value and standard deviation.
# (C) Stavros Akras

import numpy as np
from satellite import save_FITSimages_script as sfs 
from satellite import statistics_calculations_script as scs


        
def lineratios(flux2D,line_names,lines_available,lines_radial,param_estimated,param_requered,param_mod_name,param_model_values,hdr):
    index_size=param_mod_name.index("total_num_pixels_horiz")
    sizex=param_model_values[index_size]
    sizey=param_model_values[index_size]
    
    
    class lineratios:
        logHeI5876s_HI6563s=np.zeros((sizex, sizey))
        logHeI6678s_HI6563s=np.zeros((sizex, sizey))
        logHeII5412s_HI4861s=np.zeros((sizex, sizey))
        logHeII4686s_HI4861s=np.zeros((sizex, sizey))
        logHeI5876s_HeII4686s=np.zeros((sizex, sizey))
        logHeI5876s_HeII5412s=np.zeros((sizex, sizey))
        logNII6583s_HI6563s=np.zeros((sizex, sizey)) 
        logNII6548sNII6583s_HI6563s=np.zeros((sizex, sizey))
        logNII6548sNII6583s_NII5755s=np.zeros((sizex, sizey))
        logNII6548sNII6583s_OIII4959sOIII5007s=np.zeros((sizex, sizey))
        logNI5199s_HI4861s=np.zeros((sizex, sizey))
        logSII6716sSII6731s_HI6563s=np.zeros((sizex, sizey))
        logSII6716s_S26731s=np.zeros((sizex, sizey))
        logSII6716sSII6731s_SIII6312sSIII9069s=np.zeros((sizex, sizey))     
        logO16300sO16363s_HI6563s=np.zeros((sizex, sizey))
        logO16300sO16363s_O15577s=np.zeros((sizex, sizey))               
        logO16300sO16363s_O34959sO35007s=np.zeros((sizex, sizey))
        logO16300sO16363s_O27320sO27330s=np.zeros((sizex, sizey))     
        logOI6300sOI6363s_OII3727sOII3729s=np.zeros((sizex, sizey))     
        logOI6300s_HI6563s=np.zeros((sizex, sizey))                           
        logOIII5007s_HI4861s=np.zeros((sizex, sizey))                           
        logOIII4959sOIII5007s_HI4861s=np.zeros((sizex, sizey))                
        logOIII4959sOIII5007s_OIII4363s=np.zeros((sizex, sizey))                
        logOII3727sOII3729s_HI4861s=np.zeros((sizex, sizey))                
        logOII3727sOII3729s_OIII4959sOIII5007s=np.zeros((sizex, sizey))     
        logOII7320sOII7330s_OIII4959sOIII5007s=np.zeros((sizex, sizey))      
        logArVI4712s_ArVI4740s=np.zeros((sizex, sizey))                           
        logArVI4712sArVI4740s_HI4861s=np.zeros((sizex, sizey))              
        logNeIII3868sNeIII3967s_HI4861s=np.zeros((sizex, sizey))              
        logClIII5517s_ClIII5538s=np.zeros((sizex, sizey))                         
        logClIII5517sClIII5538s_HI4861s=np.zeros((sizex, sizey))              
        logCI8727s_HI6563s=np.zeros((sizex, sizey))                           
        logCII6461s_HI6563s=np.zeros((sizex, sizey))                          
        
        
    
    class lineratios_error:
        logHeI5876s_HI6563s=np.zeros((sizex, sizey))
        logHeI6678s_HI6563s=np.zeros((sizex, sizey))
        logHeII5412s_HI4861s=np.zeros((sizex, sizey))
        logHeII4686s_HI4861s=np.zeros((sizex, sizey))
        logHeI5876s_HeII4686s=np.zeros((sizex, sizey))
        logHeI5876s_HeII5412s=np.zeros((sizex, sizey))
        logNII6583s_HI6563s=np.zeros((sizex, sizey)) 
        logNII6548sNII6583s_HI6563s=np.zeros((sizex, sizey))
        logNII6548sNII6583s_NII5755s=np.zeros((sizex, sizey))
        logNII6548sNII6583s_OIII4959sOIII5007s=np.zeros((sizex, sizey))
        logNI5199s_HI4861s=np.zeros((sizex, sizey))        
        logSII6716sSII6731s_HI6563s=np.zeros((sizex, sizey))
        logSII6716s_S26731s=np.zeros((sizex, sizey))
        logSII6716sSII6731s_SIII6312sSIII9069s=np.zeros((sizex, sizey)) 
        logO16300sO16363s_HI6563s=np.zeros((sizex, sizey))
        logO16300sO16363s_O15577s=np.zeros((sizex, sizey))  
        logO16300sO16363s_O34959sO35007s=np.zeros((sizex, sizey))
        logO16300sO16363s_O27320sO27330s=np.zeros((sizex, sizey))     
        logOI6300sOI6363s_OII3727sOII3729s=np.zeros((sizex, sizey))  
        logOI6300s_HI6563s=np.zeros((sizex, sizey))                           
        logOIII5007s_HI4861s=np.zeros((sizex, sizey))  
        logOIII4959sOIII5007s_HI4861s=np.zeros((sizex, sizey))                
        logOIII4959sOIII5007s_OIII4363s=np.zeros((sizex, sizey))
        logOII3727sOII3729s_HI4861s=np.zeros((sizex, sizey))                
        logOII3727sOII3729s_OIII4959sOIII5007s=np.zeros((sizex, sizey))     
        logOII7320sOII7330s_OIII4959sOIII5007s=np.zeros((sizex, sizey)) 
        logArVI4712s_ArVI4740s=np.zeros((sizex, sizey))                           
        logArVI4712sArVI4740s_HI4861s=np.zeros((sizex, sizey))  
        logNeIII3868sNeIII3967s_HI4861s=np.zeros((sizex, sizey))              
        logClIII5517s_ClIII5538s=np.zeros((sizex, sizey))                         
        logClIII5517sClIII5538s_HI4861s=np.zeros((sizex, sizey))
        logCI8727s_HI6563s=np.zeros((sizex, sizey))                           
        logCII6461s_HI6563s=np.zeros((sizex, sizey))

###################################################################
### computes the line ratios and save the values in 2D arrays
###################################################################

    for j in range(0,sizex):
        for i in range(0,sizey):
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.HeIa_5876[i,j]>0):
                if(flux2D.HeIa_5876[i,j]/flux2D.Ha_6563[i,j]<100000):
                    lineratios.logHeI5876s_HI6563s[i,j]=np.log10(flux2D.HeIa_5876[i,j]/flux2D.Ha_6563[i,j])
    
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.HeIb_6678[i,j]>0):
                if(flux2D.HeIb_6678[i,j]/flux2D.Ha_6563[i,j]<100000):
                    lineratios.logHeI6678s_HI6563s[i,j]=np.log10(flux2D.HeIb_6678[i,j]/flux2D.Ha_6563[i,j])
                                                                                         
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.HeIIb_5412[i,j]>0):
                if(flux2D.HeIIb_5412[i,j]/flux2D.Hb_4861[i,j]<100000):
                    lineratios.logHeII5412s_HI4861s[i,j]=np.log10(flux2D.HeIIb_5412[i,j]/flux2D.Hb_4861[i,j])
                                                                                        
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.HeIIa_4686[i,j]>0):
                if(flux2D.HeIIa_4686[i,j]/flux2D.Hb_4861[i,j]<100000):
                    lineratios.logHeII4686s_HI4861s[i,j]=np.log10(flux2D.HeIIa_4686[i,j]/flux2D.Hb_4861[i,j])
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.HeIIa_4686[i,j]>0 and flux2D.HeIa_5876[i,j]>0):
                if(flux2D.HeIa_5876[i,j]/flux2D.HeIIa_4686[i,j]<100000):
                    lineratios.logHeI5876s_HeII4686s[i,j]=np.log10(flux2D.HeIa_5876[i,j]/flux2D.HeIIa_4686[i,j])
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.HeIIb_5412[i,j]>0 and flux2D.HeIa_5876[i,j]>0):
                if(flux2D.HeIa_5876[i,j]/flux2D.HeIIb_5412[i,j]<100000):
                    lineratios.logHeI5876s_HeII5412s[i,j]=np.log10(flux2D.HeIa_5876[i,j]/flux2D.HeIIb_5412[i,j])    
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.NIIc_6584[i,j]>0):
                if(flux2D.NIIc_6584[i,j]/flux2D.Ha_6563[i,j]<100000):
                    lineratios.logNII6583s_HI6563s[i,j]=np.log10(flux2D.NIIc_6584[i,j]/flux2D.Ha_6563[i,j])    
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.NIIb_6548[i,j]>0 and flux2D.NIIc_6584[i,j]>0):
                if((flux2D.NIIb_6548[i,j]+flux2D.NIIc_6584[i,j])/flux2D.Ha_6563[i,j]<100000):
                    lineratios.logNII6548sNII6583s_HI6563s[i,j]=np.log10((flux2D.NIIb_6548[i,j]+flux2D.NIIc_6584[i,j])/flux2D.Ha_6563[i,j]) 
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.NIIa_5755[i,j]>0 and flux2D.NIIb_6548[i,j]>0 and flux2D.NIIc_6584[i,j]>0):
                if((flux2D.NIIb_6548[i,j]+flux2D.NIIc_6584[i,j])/flux2D.NIIa_5755[i,j]<100000):
                    lineratios.logNII6548sNII6583s_NII5755s[i,j]=np.log10((flux2D.NIIb_6548[i,j]+flux2D.NIIc_6584[i,j])/flux2D.NIIa_5755[i,j])  
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.NIIb_6548[i,j]>0 and flux2D.NIIc_6584[i,j]>0 and flux2D.OIIIb_4959[i,j]>0 and flux2D.OIIIc_5007[i,j]>0):                                
                if((flux2D.NIIb_6548[i,j]+flux2D.NIIc_6584[i,j])/(flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])<100000):
                    lineratios.logNII6548sNII6583s_OIII4959sOIII5007s[i,j]=np.log10((flux2D.NIIb_6548[i,j]+flux2D.NIIc_6584[i,j])/(flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j]))    
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.NI_5199[i,j]>0):                                
                if(flux2D.NI_5199[i,j]/flux2D.Hb_4861[i,j]<100000):
                    lineratios.logNI5199s_HI4861s[i,j]=np.log10(flux2D.NI_5199[i,j]/flux2D.Hb_4861[i,j])
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.SIIa_6716[i,j]>0 and flux2D.SIIb_6731[i,j]>0):                                
                if((flux2D.SIIa_6716[i,j]+flux2D.SIIb_6731[i,j])/flux2D.Ha_6563[i,j]<100000):
                    lineratios.logSII6716sSII6731s_HI6563s[i,j]=np.log10((flux2D.SIIa_6716[i,j]+flux2D.SIIb_6731[i,j])/flux2D.Ha_6563[i,j])
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.SIIa_6716[i,j]>0 and flux2D.SIIb_6731[i,j]>0):                                
                if(flux2D.SIIa_6716[i,j]/flux2D.SIIb_6731[i,j]<100000):
                    lineratios.logSII6716s_S26731s[i,j]=np.log10(flux2D.SIIa_6716[i,j]/flux2D.SIIb_6731[i,j])
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.SIIa_6716[i,j]>0 and flux2D.SIIb_6731[i,j]>0 and flux2D.SIIIa_6312[i,j]>0 and flux2D.SIIIb_9069[i,j]>0):                                
                if((flux2D.SIIa_6716[i,j]+flux2D.SIIb_6731[i,j])/(flux2D.SIIIa_6312[i,j]+flux2D.SIIIb_9069[i,j])<100000):
                    lineratios.logSII6716sSII6731s_SIII6312sSIII9069s[i,j]=np.log10((flux2D.SIIa_6716[i,j]+flux2D.SIIb_6731[i,j])/(flux2D.SIIIa_6312[i,j]+flux2D.SIIIb_9069[i,j]))
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIb_6300[i,j]>0 and flux2D.OIc_6363[i,j]>0):                                
                if((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/flux2D.Ha_6563[i,j]<100000):
                    lineratios.logO16300sO16363s_HI6563s[i,j]=np.log10((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/flux2D.Ha_6563[i,j])
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIb_6300[i,j]>0 and flux2D.OIc_6363[i,j]>0 and flux2D.OIa_5577[i,j]>0):                                
                if((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/flux2D.OIa_5577[i,j] <100000):
                    lineratios.logO16300sO16363s_O15577s[i,j]=np.log10((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/flux2D.OIa_5577[i,j])    
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIb_6300[i,j]>0 and flux2D.OIc_6363[i,j]>0 and flux2D.OIIIb_4959[i,j]>0 and flux2D.OIIIc_5007[i,j]>0):                               
                if((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/(flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])<100000):
                    lineratios.logO16300sO16363s_O34959sO35007s[i,j]=np.log10((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/(flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])) 
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIb_6300[i,j]>0 and flux2D.OIc_6363[i,j]>0 and flux2D.OIIc_7320[i,j]>0 and flux2D.OIId_7330[i,j]>0):                               
                if((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/(flux2D.OIIc_7320[i,j]+flux2D.OIId_7330[i,j])<100000):
                    lineratios.logO16300sO16363s_O27320sO27330s[i,j]=np.log10((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/(flux2D.OIIc_7320[i,j]+flux2D.OIId_7330[i,j]))    
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIb_6300[i,j]>0 and flux2D.OIc_6363[i,j]>0 and flux2D.OIIa_3727[i,j]>0 and flux2D.OIIb_3729[i,j]>0):                               
                if((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/(flux2D.OIIa_3727[i,j]+flux2D.OIIb_3729[i,j])<100000):
                    lineratios.logOI6300sOI6363s_OII3727sOII3729s[i,j]=np.log10((flux2D.OIb_6300[i,j]+flux2D.OIc_6363[i,j])/(flux2D.OIIa_3727[i,j]+flux2D.OIIb_3729[i,j]))    
             
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIb_6300[i,j]>0):                            
                if(flux2D.OIb_6300[i,j]/flux2D.Ha_6563[i,j] <100000):
                    lineratios.logOI6300s_HI6563s[i,j]=np.log10(flux2D.OIb_6300[i,j]/flux2D.Ha_6563[i,j]) 
           
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIIIc_5007[i,j]>0):                            
                if(flux2D.OIIIc_5007[i,j]/flux2D.Hb_4861[i,j]<100000):
                    lineratios.logOIII5007s_HI4861s[i,j]=np.log10(flux2D.OIIIc_5007[i,j]/flux2D.Hb_4861[i,j]) 
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIIIb_4959[i,j]>0 and flux2D.OIIIc_5007[i,j]>0):                            
                if((flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])/flux2D.Hb_4861[i,j]<100000):
                    lineratios.logOIII4959sOIII5007s_HI4861s[i,j]=np.log10((flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])/flux2D.Hb_4861[i,j])
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIIIb_4959[i,j]>0 and flux2D.OIIIc_5007[i,j]>0 and flux2D.OIIIa_4363[i,j]>0):                            
                if((flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])/flux2D.OIIIa_4363[i,j]<100000):
                    lineratios.logOIII4959sOIII5007s_OIII4363s[i,j]=np.log10((flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])/flux2D.OIIIa_4363[i,j])
      
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIIa_3727[i,j]>0 and flux2D.OIIb_3729[i,j]>0):                            
                if((flux2D.OIIa_3727[i,j]+flux2D.OIIb_3729[i,j])/flux2D.Hb_4861[i,j]<100000):
                    lineratios.logOII3727sOII3729s_HI4861s[i,j]=np.log10((flux2D.OIIa_3727[i,j]+flux2D.OIIb_3729[i,j])/flux2D.Hb_4861[i,j])
      
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIIa_3727[i,j]>0 and flux2D.OIIb_3729[i,j]>0 and flux2D.OIIIb_4959[i,j]>0 and flux2D.OIIIc_5007[i,j]>0):                            
                if((flux2D.OIIa_3727[i,j]+flux2D.OIIb_3729[i,j])/(flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])<100000):
                    lineratios.logOII3727sOII3729s_OIII4959sOIII5007s[i,j]=np.log10((flux2D.OIIa_3727[i,j]+flux2D.OIIb_3729[i,j])/(flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j]))
        
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.OIIc_7320[i,j]>0 and flux2D.OIId_7330[i,j]>0 and flux2D.OIIIb_4959[i,j]>0 and flux2D.OIIIc_5007[i,j]>0):                            
                if((flux2D.OIIc_7320[i,j]+flux2D.OIId_7330[i,j])/(flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j])<100000):
                    lineratios.logOII7320sOII7330s_OIII4959sOIII5007s[i,j]=np.log10((flux2D.OIIc_7320[i,j]+flux2D.OIId_7330[i,j])/(flux2D.OIIIb_4959[i,j]+flux2D.OIIIc_5007[i,j]))
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.ArIVa_4712[i,j]>0 and flux2D.ArIVb_4740[i,j]>0):                            
                if(flux2D.ArIVa_4712[i,j]/flux2D.ArIVb_4740[i,j] <100000):
                    lineratios.logArVI4712s_ArVI4740s[i,j]=np.log10(flux2D.ArIVa_4712[i,j]/flux2D.ArIVb_4740[i,j]) 
           
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.ArIVa_4712[i,j]>0 and flux2D.ArIVb_4740[i,j]>0):                            
                if((flux2D.ArIVa_4712[i,j]+flux2D.ArIVb_4740[i,j])/flux2D.Hb_4861[i,j] <100000):
                    lineratios.logArVI4712sArVI4740s_HI4861s[i,j]=np.log10((flux2D.ArIVa_4712[i,j]+flux2D.ArIVb_4740[i,j])/flux2D.Hb_4861[i,j]) 
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.NeIIIa_3868[i,j]>0 and flux2D.NeIIIb_3967[i,j]>0):                            
                if((flux2D.NeIIIa_3868[i,j]+flux2D.NeIIIb_3967[i,j])/flux2D.Hb_4861[i,j] <100000):
                    lineratios.logNeIII3868sNeIII3967s_HI4861s[i,j]=np.log10((flux2D.NeIIIa_3868[i,j]+flux2D.NeIIIb_3967[i,j])/flux2D.Hb_4861[i,j]) 
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.ClIIIa_5517[i,j]>0 and flux2D.ClIIIb_5538[i,j]>0):                            
                if((flux2D.ClIIIa_5517[i,j]+flux2D.ClIIIb_5538[i,j])/flux2D.Hb_4861[i,j]<100000):
                    lineratios.logClIII5517sClIII5538s_HI4861s[i,j]=np.log10((flux2D.ClIIIa_5517[i,j]+flux2D.ClIIIb_5538[i,j])/flux2D.Hb_4861[i,j]) 
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.ClIIIa_5517[i,j]>0 and flux2D.ClIIIb_5538[i,j]>0):                            
                if(flux2D.ClIIIa_5517[i,j]/flux2D.ClIIIb_5538[i,j] <100000):
                    lineratios.logClIII5517s_ClIII5538s[i,j]=np.log10(flux2D.ClIIIa_5517[i,j]/flux2D.ClIIIb_5538[i,j]) 
            
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.CI_8727[i,j]>0):                            
                if(flux2D.CI_8727[i,j]/flux2D.Ha_6563[i,j]<100000):
                    lineratios.logCI8727s_HI6563s[i,j]=np.log10(flux2D.CI_8727[i,j]/flux2D.Ha_6563[i,j])
                
            if(flux2D.Ha_6563[i,j]>0.0 and flux2D.Hb_4861[i,j]>0.0 and flux2D.CII_6461[i,j]>0):                            
                if(flux2D.CII_6461[i,j]/flux2D.Ha_6563[i,j]  <100000):
                    lineratios.logCII6461s_HI6563s[i,j]=np.log10(flux2D.CII_6461[i,j]/flux2D.Ha_6563[i,j])    


###################################################################
### creates and return FITs images for each emission line ratio
###################################################################
            
    sfs.save_image(lineratios.logHeI5876s_HI6563s,hdr,'log(HeI_5876s_HI_6563s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logHeI6678s_HI6563s,hdr,'log(HeI_6678s_HI_6563s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logHeII5412s_HI4861s,hdr,'log(HeII_5412s_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logHeII4686s_HI4861s,hdr,'log(HeII_4686s_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logHeI5876s_HeII4686s,hdr,'log(HeI_5876s_HeII_4686s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logHeI5876s_HeII5412s,hdr,'log(HeI_5876s_HeII_5412s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logNII6583s_HI6563s,hdr,'log(N2_6583s_HI_6563s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logNII6548sNII6583s_HI6563s,hdr,'log(N2_6548s+N2_6583s)_HI_6563s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logNII6548sNII6583s_NII5755s,hdr,'log((N2_6548s+N2_6583s)_N2_5755s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logNII6548sNII6583s_OIII4959sOIII5007s,hdr,'log((N2_6548s+N2_6583s)_(O3_4959s+O3_5007s))',param_estimated,param_requered)
    sfs.save_image(lineratios.logNI5199s_HI4861s,hdr,'log(N1_5199s_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logSII6716sSII6731s_HI6563s,hdr,'log((S2_6716s+S2_6731s)_HI_6563s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logSII6716s_S26731s,hdr,'log(S2_6716s_S2_6731s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logSII6716sSII6731s_SIII6312sSIII9069s,hdr,'log((S2_6716s+S2_6731s)_(S3_6312s+S3_9069s))',param_estimated,param_requered)
    sfs.save_image(lineratios.logO16300sO16363s_HI6563s,hdr,'log((O1_6300s+O1_6363s)_HI_6563s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logO16300sO16363s_O15577s,hdr,'log((O1_6300s+O1_6363s)_O1_5577s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logO16300sO16363s_O34959sO35007s,hdr,'log((O1_6300s+O1_6363s)_(O3_4959s+O3_5007s))',param_estimated,param_requered)
    sfs.save_image(lineratios.logO16300sO16363s_O27320sO27330s,hdr,'log((O1_6300s+O1_6363s)_(O2_7320s+O2_7330s))',param_estimated,param_requered)
    sfs.save_image(lineratios.logOI6300sOI6363s_OII3727sOII3729s,hdr,'log((O1_6300s+O1_6363s)_(O2_3727s+O2_3729s))',param_estimated,param_requered)
    sfs.save_image(lineratios.logOI6300s_HI6563s,hdr,'log(O1_6300s_HI_6563s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logOIII5007s_HI4861s,hdr,'log(O3_5007s_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logOIII4959sOIII5007s_HI4861s,hdr,'log((O3_4959s+O3_5007s)_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logOIII4959sOIII5007s_OIII4363s,hdr,'log((O3_4959s+O3_5007s)_O3_4363s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logOII3727sOII3729s_HI4861s,hdr,'log((O2_3727s+O2_3729s)_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logOII3727sOII3729s_OIII4959sOIII5007s,hdr,'log((O2_3727s+O2_3729s)_(O3_4959s+O3_5007s))',param_estimated,param_requered)
    sfs.save_image(lineratios.logOII7320sOII7330s_OIII4959sOIII5007s,hdr,'log((O2_7320s+O2_7330s)_(O3_4959s+O3_5007s))',param_estimated,param_requered)
    sfs.save_image(lineratios.logArVI4712s_ArVI4740s,hdr,'log(Ar4_4712s_Ar4_4740s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logArVI4712sArVI4740s_HI4861s,hdr,'log((Ar4_4712s+Ar4_4740s)_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logNeIII3868sNeIII3967s_HI4861s,hdr,'log((Ne3_3868s+Ne3_3967s)_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logClIII5517sClIII5538s_HI4861s,hdr,'log((Cl3_5517s+Cl3_5538s)_HI_4861s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logClIII5517s_ClIII5538s,hdr,'log(Cl3_5517s_Cl3_5538s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logCI8727s_HI6563s,hdr,'log(C1_8727s_HI_6563s)',param_estimated,param_requered)
    sfs.save_image(lineratios.logCII6461s_HI6563s,hdr,'log(C2_6461s_HI_6563s)',param_estimated,param_requered)


###################################################################
### computes Parameter Name, Nspaxels, 5% value, Q1 value, median value, 
### Q3 value, 95% value, mean value and standard deviation
### and save them in the statisticnumbers_2.txt ascii file
###################################################################

    file12 = open('statisticnumbers_2.txt','w')
    print('{:45}'.format('Parameter Name'), '{:10}'.format('Nspaxels'),'{:10}'.format('5% value'),'{:10}'.format('Q1 value'), '{:12}'.format('median value'), '{:10}'.format('Q3 value'), '{:10}'.format('95% value'), '{:10}'.format('mean value'), '{:10}'.format('sigma'), file=file12)
    print("#############################################################################################", file=file12)
    
    
    index=param_estimated.index('log(HeI_5876s_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logHeI5876s_HI6563s,sizex,sizey)
        print('{:45}'.format('log(HeI_5876s_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log(HeI_6678s_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logHeI6678s_HI6563s,sizex,sizey)
        print('{:45}'.format('log(HeI_6678s_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(HeII_5412s_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logHeII5412s_HI4861s,sizex,sizey)
        print('{:45}'.format('log(HeII_5412s_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(HeII_4686s_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logHeII4686s_HI4861s,sizex,sizey)
        print('{:45}'.format('log(HeII_4686s_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(HeI_5876s_HeII_4686s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logHeI5876s_HeII4686s,sizex,sizey)
        print('{:45}'.format('log(HeI_5876s_HeII_4686s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log(HeI_5876s_HeII_5412s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logHeI5876s_HeII5412s,sizex,sizey)
        print('{:45}'.format('log(HeI_5876s_HeII_5412s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(N2_6583s_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logNII6583s_HI6563s,sizex,sizey)
        print('{:45}'.format('log(N2_6583s_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(N2_6548s+N2_6583s)_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logNII6548sNII6583s_HI6563s,sizex,sizey)
        print('{:45}'.format('log(N2_6548s+N2_6583s)_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((N2_6548s+N2_6583s)_N2_5755s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logNII6548sNII6583s_NII5755s,sizex,sizey)
        print('{:45}'.format('log((N2_6548s+N2_6583s)_N2_5755s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((N2_6548s+N2_6583s)_(O3_4959s+O3_5007s))')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logNII6548sNII6583s_OIII4959sOIII5007s,sizex,sizey)
        print('{:45}'.format('log((N2_6548s+N2_6583s)_(O3_4959s+O3_5007s))'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log(N1_5199s_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logNI5199s_HI4861s,sizex,sizey)
        print('{:45}'.format('log(N1_5199s_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((S2_6716s+S2_6731s)_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logSII6716sSII6731s_HI6563s,sizex,sizey)
        print('{:45}'.format('log((S2_6716s+S2_6731s)_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(S2_6716s_S2_6731s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logSII6716s_S26731s,sizex,sizey)
        print('{:45}'.format('log(S2_6716s_S2_6731s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((S2_6716s+S2_6731s)_(S3_6312s+S3_9069s))')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logSII6716sSII6731s_SIII6312sSIII9069s,sizex,sizey)
        print('{:45}'.format('log((S2_6716s+S2_6731s)_(S3_6312s+S3_9069s))'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((O1_6300s+O1_6363s)_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logO16300sO16363s_HI6563s,sizex,sizey)
        print('{:45}'.format('log((O1_6300s+O1_6363s)_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log((O1_6300s+O1_6363s)_O1_5577s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logO16300sO16363s_O15577s,sizex,sizey)
        print('{:45}'.format('log((O1_6300s+O1_6363s)_O1_5577s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log((O1_6300s+O1_6363s)_(O3_4959s+O3_5007s))')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logO16300sO16363s_O34959sO35007s,sizex,sizey)
        print('{:45}'.format('log((O1_6300s+O1_6363s)_(O3_4959s+O3_5007s))'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((O1_6300s+O1_6363s)_(O2_7320s+O2_7330s))')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logO16300sO16363s_O27320sO27330s,sizex,sizey)
        print('{:45}'.format('log((O1_6300s+O1_6363s)_(O2_7320s+O2_7330s))'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((O1_6300s+O1_6363s)_(O2_3727s+O2_3729s))')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logOI6300sOI6363s_OII3727sOII3729s,sizex,sizey)
        print('{:45}'.format('log((O1_6300s+O1_6363s)_(O2_3727s+O2_3729s))'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(O1_6300s_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logOI6300s_HI6563s,sizex,sizey)
        print('{:45}'.format('log(O1_6300s_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(O3_5007s_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logOIII5007s_HI4861s,sizex,sizey)
        print('{:45}'.format('log(O3_5007s_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((O3_4959s+O3_5007s)_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logOIII4959sOIII5007s_HI4861s,sizex,sizey)
        print('{:45}'.format('log((O3_4959s+O3_5007s)_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((O3_4959s+O3_5007s)_O3_4363s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logOIII4959sOIII5007s_OIII4363s,sizex,sizey)
        print('{:45}'.format('log((O3_4959s+O3_5007s)_O3_4363s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log((O2_3727s+O2_3729s)_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logOII3727sOII3729s_HI4861s,sizex,sizey)
        print('{:45}'.format('log((O2_3727s+O2_3729s)_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log((O2_3727s+O2_3729s)_(O3_4959s+O3_5007s))')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logOII3727sOII3729s_OIII4959sOIII5007s,sizex,sizey)
        print('{:45}'.format('log((O2_3727s+O2_3729s)_(O3_4959s+O3_5007s))'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((O2_7320s+O2_7330s)_(O3_4959s+O3_5007s))')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logOII7320sOII7330s_OIII4959sOIII5007s,sizex,sizey)
        print('{:45}'.format('log((O2_7320s+O2_7330s)_(O3_4959s+O3_5007s))'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(Ar4_4712s_Ar4_4740s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logArVI4712s_ArVI4740s,sizex,sizey)
        print('{:45}'.format('log(Ar4_4712s_Ar4_4740s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log((Ar4_4712s+Ar4_4740s)_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logArVI4712sArVI4740s_HI4861s,sizex,sizey)
        print('{:45}'.format('log((Ar4_4712s+Ar4_4740s)_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log((Ne3_3868s+Ne3_3967s)_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logNeIII3868sNeIII3967s_HI4861s,sizex,sizey)
        print('{:45}'.format('log((Ne3_3868s+Ne3_3967s)_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log((Cl3_5517s+Cl3_5538s)_HI_4861s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logClIII5517sClIII5538s_HI4861s,sizex,sizey)
        print('{:45}'.format('log((Cl3_5517s+Cl3_5538s)_HI_4861s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(Cl3_5517s_Cl3_5538s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logClIII5517s_ClIII5538s,sizex,sizey)
        print('{:45}'.format('log(Cl3_5517s_Cl3_5538s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
    
    index=param_estimated.index('log(C1_8727s_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logCI8727s_HI6563s,sizex,sizey)
        print('{:45}'.format('log(C1_8727s_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)

    index=param_estimated.index('log(C2_6461s_HI_6563s)')
    if param_requered[index]=="yes":
        pxl_num,percent5,percentQ1,medianvalue,percentQ3,percent95,meanvalue,sigmavalue= scs.statistic_numbers(lineratios.logCII6461s_HI6563s,sizex,sizey)
        print('{:45}'.format('log(C2_6461s_HI_6563s)'), '{:10.0f}'.format(pxl_num),'{:8.2e}'.format(percent5), '{:8.2e}'.format(percentQ1),'{:10.2e}'.format(medianvalue), '{:8.2e}'.format(percent95),'{:8.2e}'.format(percentQ3), '{:8.2e}'.format(meanvalue),'{:8.2e}'.format(sigmavalue), file=file12)
   
    file12.close() 

    return      

        
