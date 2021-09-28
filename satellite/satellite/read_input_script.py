# read_input_script.py:
# Reads the names of the emission line provided (yes or no) in the input.txt file. 
# Reads the fits 2D flux and error images of emission lines and save them in 2D arrays. 
# These new 2D arrays are bigger that the first ones because the center of the nebula 
# is moved to the center of the new map and some rows/columns of pixels have been added.
# (C) Stavros Akras

import numpy as np
from astropy.io.fits import getdata




def read_input_images(line_name,requested,param,flux_or_error):
    
    global dataHanew,dataHbnew,dataHgnew,dataHdnew,dataHeIanew,dataHeIbnew,dataHeIIanew,dataHeIIbnew,dataNIIanew,dataNIIbnew,dataNIIcnew,dataNInew,dataOIIIanew,dataOIIIbnew,dataOIIIcnew,dataOIIanew,dataOIIbnew,dataOIIcnew,dataOIIdnew,dataOIanew,dataOIbnew,dataOIcnew,dataSIIanew,dataSIIbnew,dataSIIIanew,dataSIIIbnew,dataClIIIanew,dataClIIIbnew,dataArIIInew,dataArIVanew,dataArIVbnew,dataCInew,dataCIInew,dataNeIIIanew,dataNeIIIbnew,hdr

##########################################################################
### extra rows/columns of pixels are create and then added in the raw maps
### the number of rows and columns are given by the user.
##########################################################################
    z1=np.zeros((param[3],param[7]))  # above
    z2=np.zeros((param[4],param[7]))  # below
    z3=np.zeros((param[8],param[5]))  # left 
    z4=np.zeros((param[8],param[6]))  #right
    

    
    if flux_or_error=="fluxes":
        for i in range(0,len(line_name)):
            if line_name[i] == 'HI_6563s' and requested[i]=='yes':
                dataHa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHa=np.nan_to_num(dataHa)
                addupdown1=np.concatenate((z2,dataHa,z1), axis=0)
                dataHanew=np.concatenate((z3,addupdown1,z4), axis=1)
            sx=hdr['NAXIS1']
            sy=hdr['NAXIS2']
            if line_name[i] == 'HI_4861s' and requested[i]=='yes':
                dataHb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHb=np.nan_to_num(dataHb)
                addupdown1=np.concatenate((z2,dataHb,z1), axis=0)
                dataHbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HI_4861s' and requested[i]=='no':
                dataHb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHb,z1), axis=0)
                dataHbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HI_4340s' and requested[i]=='yes':
                dataHg, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHg=np.nan_to_num(dataHg)
                addupdown1=np.concatenate((z2,dataHg,z1), axis=0)
                dataHgnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HI_4340s' and requested[i]=='no':
                dataHg=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHg,z1), axis=0)
                dataHgnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HI_4101s' and requested[i]=='yes':
                dataHd, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHd=np.nan_to_num(dataHd)
                addupdown1=np.concatenate((z2,dataHd,z1), axis=0)
                dataHdnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HI_4101s' and requested[i]=='no':
                dataHd=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHd,z1), axis=0)
                dataHdnew=np.concatenate((z3,addupdown1,z4), axis=1)

            if line_name[i] == 'HeI_5876s' and requested[i]=='yes':
                dataHeIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHeIa=np.nan_to_num(dataHeIa)
                addupdown1=np.concatenate((z2,dataHeIa,z1), axis=0)
                dataHeIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HeI_5876s' and requested[i]=='no':
                dataHeIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIa,z1), axis=0)
                dataHeIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HeI_6678s' and requested[i]=='yes':
                dataHeIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHeIb=np.nan_to_num(dataHeIb)
                addupdown1=np.concatenate((z2,dataHeIb,z1), axis=0)
                dataHeIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HeI_6678s' and requested[i]=='no':
                dataHeIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIb,z1), axis=0)
                dataHeIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HeII_4686s' and requested[i]=='yes':    
                dataHeIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHeIIa=np.nan_to_num(dataHeIIa)
                addupdown1=np.concatenate((z2,dataHeIIa,z1), axis=0)
                dataHeIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HeII_4686s' and requested[i]=='no':
                dataHeIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIIa,z1), axis=0)
                dataHeIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HeII_5412s' and requested[i]=='yes':
                dataHeIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHeIIb=np.nan_to_num(dataHeIIb)
                addupdown1=np.concatenate((z2,dataHeIIb,z1), axis=0)
                dataHeIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HeII_5412s' and requested[i]=='no':
                dataHeIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIIb,z1), axis=0)
                dataHeIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            
            if line_name[i] == 'N2_5755s' and requested[i]=='yes':
                dataNIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNIIa=np.nan_to_num(dataNIIa)
                addupdown1=np.concatenate((z2,dataNIIa,z1), axis=0)
                dataNIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'N2_5755s' and requested[i]=='no':
                dataNIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIa,z1), axis=0)
                dataNIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'N2_6548s' and requested[i]=='yes':
                dataNIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNIIb=np.nan_to_num(dataNIIb)
                addupdown1=np.concatenate((z2,dataNIIb,z1), axis=0)
                dataNIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'N2_6548s' and requested[i]=='no':
                dataNIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIb,z1), axis=0)
                dataNIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'N2_6583s' and requested[i]=='yes':
                dataNIIc, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNIIc=np.nan_to_num(dataNIIc)
                addupdown1=np.concatenate((z2,dataNIIc,z1), axis=0)
                dataNIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'N2_6583s' and requested[i]=='no':
                dataNIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIc,z1), axis=0)
                dataNIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'N1_5199s' and requested[i]=='yes':
                dataNI, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNI=np.nan_to_num(dataNI)
                addupdown1=np.concatenate((z2,dataNI,z1), axis=0)
                dataNInew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'N1_5199s' and requested[i]=='no':
                dataNI=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNI,z1), axis=0)
                dataNInew=np.concatenate((z3,addupdown1,z4), axis=1)
            
            if line_name[i] == 'O3_4363s' and requested[i]=='yes':
                dataOIIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIIa=np.nan_to_num(dataOIIIa)
                addupdown1=np.concatenate((z2,dataOIIIa,z1), axis=0)
                dataOIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O3_4363s' and requested[i]=='no':
                dataOIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIa,z1), axis=0)
                dataOIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O3_4959s' and requested[i]=='yes':
                dataOIIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIIb=np.nan_to_num(dataOIIIb)
                addupdown1=np.concatenate((z2,dataOIIIb,z1), axis=0)
                dataOIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O3_4959s' and requested[i]=='no':
                dataOIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIb,z1), axis=0)
                dataOIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O3_5007s' and requested[i]=='yes':
                dataOIIIc, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIIc=np.nan_to_num(dataOIIIc)
                addupdown1=np.concatenate((z2,dataOIIIc,z1), axis=0)
                dataOIIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O3_5007s' and requested[i]=='no':
                dataOIIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIc,z1), axis=0)
                dataOIIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O2_3727s' and requested[i]=='yes':
                dataOIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIa=np.nan_to_num(dataOIIa)
                addupdown1=np.concatenate((z2,dataOIIa,z1), axis=0)
                dataOIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O2_3727s' and requested[i]=='no':
                dataOIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIa,z1), axis=0)
                dataOIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O2_3729s' and requested[i]=='yes':
                dataOIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIb=np.nan_to_num(dataOIIb)
                addupdown1=np.concatenate((z2,dataOIIb,z1), axis=0)
                dataOIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O2_3729s' and requested[i]=='no':
                dataOIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIb,z1), axis=0)
                dataOIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O2_7320s' and requested[i]=='yes':
                dataOIIc, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIc=np.nan_to_num(dataOIIc)
                addupdown1=np.concatenate((z2,dataOIIc,z1), axis=0)
                dataOIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O2_7320s' and requested[i]=='no':
                dataOIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIc,z1), axis=0)
                dataOIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O2_7330s' and requested[i]=='yes':
                dataOIId, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIId=np.nan_to_num(dataOIId)
                addupdown1=np.concatenate((z2,dataOIId,z1), axis=0)
                dataOIIdnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O2_7330s' and requested[i]=='no':
                dataOIId=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIId,z1), axis=0)
                dataOIIdnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O1_5577s' and requested[i]=='yes':
                dataOIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIa=np.nan_to_num(dataOIa)
                addupdown1=np.concatenate((z2,dataOIa,z1), axis=0)
                dataOIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O1_5577s' and requested[i]=='no':
                dataOIa=np.zeros((sy, sx))   
                addupdown1=np.concatenate((z2,dataOIa,z1), axis=0)
                dataOIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O1_6300s' and requested[i]=='yes':
                dataOIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIb=np.nan_to_num(dataOIb)
                addupdown1=np.concatenate((z2,dataOIb,z1), axis=0)
                dataOIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O1_6300s' and requested[i]=='no':
                dataOIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIb,z1), axis=0)
                dataOIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O1_6363s' and requested[i]=='yes':
                dataOIc, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIc=np.nan_to_num(dataOIc)
                addupdown1=np.concatenate((z2,dataOIc,z1), axis=0)
                dataOIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O1_6363s' and requested[i]=='no':
                dataOIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIc,z1), axis=0)
                dataOIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
                    
            if line_name[i] == 'S2_6716s' and requested[i]=='yes':
                dataSIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataSIIa=np.nan_to_num(dataSIIa)
                addupdown1=np.concatenate((z2,dataSIIa,z1), axis=0)
                dataSIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'S2_6716s' and requested[i]=='no':
                dataSIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIa,z1), axis=0)
                dataSIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'S2_6731s' and requested[i]=='yes':
                dataSIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataSIIb=np.nan_to_num(dataSIIb)
                addupdown1=np.concatenate((z2,dataSIIb,z1), axis=0)
                dataSIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'S2_6731s' and requested[i]=='no':
                dataSIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIb,z1), axis=0)
                dataSIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'S3_6312s' and requested[i]=='yes':
                dataSIIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataSIIIa=np.nan_to_num(dataSIIIa)
                addupdown1=np.concatenate((z2,dataSIIIa,z1), axis=0)
                dataSIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'S3_6312s' and requested[i]=='no':
                dataSIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIIa,z1), axis=0)
                dataSIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'S3_9069s' and requested[i]=='yes':
                dataSIIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataSIIIb=np.nan_to_num(dataSIIIb)
                addupdown1=np.concatenate((z2,dataSIIIb,z1), axis=0)
                dataSIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'S3_9069s' and requested[i]=='no':
                dataSIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIIb,z1), axis=0)
                dataSIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
                
            if line_name[i] == 'Cl3_5517s' and requested[i]=='yes':
                dataClIIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataClIIIa=np.nan_to_num(dataClIIIa)
                addupdown1=np.concatenate((z2,dataClIIIa,z1), axis=0)
                dataClIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Cl3_5517s' and requested[i]=='no':
                dataClIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataClIIIa,z1), axis=0)
                dataClIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Cl3_5538s' and requested[i]=='yes':
                dataClIIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataClIIIb=np.nan_to_num(dataClIIIb)
                addupdown1=np.concatenate((z2,dataClIIIb,z1), axis=0)
                dataClIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Cl3_5538s' and requested[i]=='no':
                dataClIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataClIIIb,z1), axis=0)
                dataClIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Ar3_7136s' and requested[i]=='yes':
                dataArIII, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataArIII=np.nan_to_num(dataArIII)
                addupdown1=np.concatenate((z2,dataArIII,z1), axis=0)
                dataArIIInew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ar3_7136s' and requested[i]=='no':
                dataArIII=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIII,z1), axis=0)
                dataArIIInew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Ar4_4712s' and requested[i]=='yes':
                dataArIVa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataArIVa=np.nan_to_num(dataArIVa)
                addupdown1=np.concatenate((z2,dataArIVa,z1), axis=0)
                dataArIVanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ar4_4712s' and requested[i]=='no':
                dataArIVa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIVa,z1), axis=0)
                dataArIVanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Ar4_4740s' and requested[i]=='yes':
                dataArIVb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataArIVb=np.nan_to_num(dataArIVb)
                addupdown1=np.concatenate((z2,dataArIVb,z1), axis=0)
                dataArIVbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ar4_4740s' and requested[i]=='no':
                dataArIVb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIVb,z1), axis=0)
                dataArIVbnew=np.concatenate((z3,addupdown1,z4), axis=1)

            if line_name[i] == 'C1_8727s' and requested[i]=='yes':
                dataCI, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataCI=np.nan_to_num(dataCI)
                addupdown1=np.concatenate((z2,dataCI,z1), axis=0)
                dataCInew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'C1_8727s' and requested[i]=='no':
                dataCI=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataCI,z1), axis=0)
                dataCInew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'C2_6461s' and requested[i]=='yes':
                dataCII, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataCII=np.nan_to_num(dataCII)
                addupdown1=np.concatenate((z2,dataCII,z1), axis=0)
                dataCIInew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'C2_6461s' and requested[i]=='no':
                dataCII=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataCII,z1), axis=0)
                dataCIInew=np.concatenate((z3,addupdown1,z4), axis=1)
    
            if line_name[i] == 'Ne3_3868s' and requested[i]=='yes':
                dataNeIIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNeIIIa=np.nan_to_num(dataNeIIIa)
                addupdown1=np.concatenate((z2,dataNeIIIa,z1), axis=0)
                dataNeIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ne3_3868s' and requested[i]=='no':
                dataNeIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNeIIIa,z1), axis=0)
                dataNeIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Ne3_3967s' and requested[i]=='yes':
                dataNeIIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNeIIIb=np.nan_to_num(dataNeIIIb)
                addupdown1=np.concatenate((z2,dataNeIIIb,z1), axis=0)
                dataNeIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ne3_3967s' and requested[i]=='no':
                dataNeIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNeIIIb,z1), axis=0)
                dataNeIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
    
    
    
    
        return dataHanew,dataHbnew,dataHgnew,dataHdnew,dataHeIanew,dataHeIbnew,dataHeIIanew,dataHeIIbnew,dataNIIanew,dataNIIbnew,dataNIIcnew,dataNInew,dataOIIIanew,dataOIIIbnew,dataOIIIcnew,dataOIIanew,dataOIIbnew,dataOIIcnew,dataOIIdnew,dataOIanew,dataOIbnew,dataOIcnew,dataSIIanew,dataSIIbnew,dataSIIIanew,dataSIIIbnew,dataClIIIanew,dataClIIIbnew,dataArIIInew,dataArIVanew,dataArIVbnew,dataCInew,dataCIInew,dataNeIIIanew,dataNeIIIbnew,hdr


    if flux_or_error=="errors":
        for i in range(0,len(line_name)):
            if line_name[i] == 'HI_6563e' and requested[i]=='yes':
                dataHa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHa=np.nan_to_num(dataHa)
                addupdown1=np.concatenate((z2,dataHa,z1), axis=0)
                dataHanew=np.concatenate((z3,addupdown1,z4), axis=1)
            sx=hdr['NAXIS1']
            sy=hdr['NAXIS2']
            if line_name[i] == 'HI_4861e' and requested[i]=='yes':
                dataHb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHb=np.nan_to_num(dataHb)
                addupdown1=np.concatenate((z2,dataHb,z1), axis=0)
                dataHbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HI_4861e' and requested[i]=='no':
                dataHb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHb,z1), axis=0)
                dataHbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HI_4340e' and requested[i]=='yes':
                dataHg, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHg=np.nan_to_num(dataHg)
                addupdown1=np.concatenate((z2,dataHg,z1), axis=0)
                dataHgnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HI_4340e' and requested[i]=='no':
                dataHg=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHg,z1), axis=0)
                dataHgnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HI_4101e' and requested[i]=='yes':
                dataHd, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHd=np.nan_to_num(dataHd)
                addupdown1=np.concatenate((z2,dataHd,z1), axis=0)
                dataHdnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HI_4101e' and requested[i]=='no':
                dataHd=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHd,z1), axis=0)
                dataHdnew=np.concatenate((z3,addupdown1,z4), axis=1)

            if line_name[i] == 'HeI_5876e' and requested[i]=='yes':
                dataHeIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHeIa=np.nan_to_num(dataHeIa)
                addupdown1=np.concatenate((z2,dataHeIa,z1), axis=0)
                dataHeIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HeI_5876e' and requested[i]=='no':
                dataHeIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIa,z1), axis=0)
                dataHeIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HeI_6678e' and requested[i]=='yes':
                dataHeIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHeIb=np.nan_to_num(dataHeIb)
                addupdown1=np.concatenate((z2,dataHeIb,z1), axis=0)
                dataHeIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HeI_6678e' and requested[i]=='no':
                dataHeIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIb,z1), axis=0)
                dataHeIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HeII_4686e' and requested[i]=='yes':    
                dataHeIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHeIIa=np.nan_to_num(dataHeIIa)
                addupdown1=np.concatenate((z2,dataHeIIa,z1), axis=0)
                dataHeIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HeII_4686e' and requested[i]=='no':
                dataHeIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIIa,z1), axis=0)
                dataHeIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'HeII_5412e' and requested[i]=='yes':
                dataHeIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataHeIIb=np.nan_to_num(dataHeIIb)
                addupdown1=np.concatenate((z2,dataHeIIb,z1), axis=0)
                dataHeIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'HeII_5412e' and requested[i]=='no':
                dataHeIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataHeIIb,z1), axis=0)
                dataHeIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            
            if line_name[i] == 'N2_5755e' and requested[i]=='yes':
                dataNIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNIIa=np.nan_to_num(dataNIIa)
                addupdown1=np.concatenate((z2,dataNIIa,z1), axis=0)
                dataNIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'N2_5755e' and requested[i]=='no':
                dataNIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIa,z1), axis=0)
                dataNIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'N2_6548e' and requested[i]=='yes':
                dataNIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNIIb=np.nan_to_num(dataNIIb)
                addupdown1=np.concatenate((z2,dataNIIb,z1), axis=0)
                dataNIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'N2_6548e' and requested[i]=='no':
                dataNIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIb,z1), axis=0)
                dataNIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'N2_6583e' and requested[i]=='yes':
                dataNIIc, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNIIc=np.nan_to_num(dataNIIc)
                addupdown1=np.concatenate((z2,dataNIIc,z1), axis=0)
                dataNIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'N2_6583e' and requested[i]=='no':
                dataNIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNIIc,z1), axis=0)
                dataNIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'N1_5199e' and requested[i]=='yes':
                dataNI, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNI=np.nan_to_num(dataNI)
                addupdown1=np.concatenate((z2,dataNI,z1), axis=0)
                dataNInew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'N1_5199e' and requested[i]=='no':
                dataNI=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNI,z1), axis=0)
                dataNInew=np.concatenate((z3,addupdown1,z4), axis=1)
            
            if line_name[i] == 'O3_4363e' and requested[i]=='yes':
                dataOIIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIIa=np.nan_to_num(dataOIIIa)
                addupdown1=np.concatenate((z2,dataOIIIa,z1), axis=0)
                dataOIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O3_4363e' and requested[i]=='no':
                dataOIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIa,z1), axis=0)
                dataOIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O3_4959e' and requested[i]=='yes':
                dataOIIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIIb=np.nan_to_num(dataOIIIb)
                addupdown1=np.concatenate((z2,dataOIIIb,z1), axis=0)
                dataOIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O3_4959e' and requested[i]=='no':
                dataOIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIb,z1), axis=0)
                dataOIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O3_5007e' and requested[i]=='yes':
                dataOIIIc, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIIc=np.nan_to_num(dataOIIIc)
                addupdown1=np.concatenate((z2,dataOIIIc,z1), axis=0)
                dataOIIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O3_5007e' and requested[i]=='no':
                dataOIIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIIc,z1), axis=0)
                dataOIIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O2_3727e' and requested[i]=='yes':
                dataOIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIa=np.nan_to_num(dataOIIa)
                addupdown1=np.concatenate((z2,dataOIIa,z1), axis=0)
                dataOIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O2_3727e' and requested[i]=='no':
                dataOIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIa,z1), axis=0)
                dataOIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O2_3729e' and requested[i]=='yes':
                dataOIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIb=np.nan_to_num(dataOIIb)
                addupdown1=np.concatenate((z2,dataOIIb,z1), axis=0)
                dataOIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O2_3729e' and requested[i]=='no':
                dataOIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIb,z1), axis=0)
                dataOIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O2_7320e' and requested[i]=='yes':
                dataOIIc, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIIc=np.nan_to_num(dataOIIc)
                addupdown1=np.concatenate((z2,dataOIIc,z1), axis=0)
                dataOIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O2_7320e' and requested[i]=='no':
                dataOIIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIIc,z1), axis=0)
                dataOIIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O2_7330e' and requested[i]=='yes':
                dataOIId, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIId=np.nan_to_num(dataOIId)
                addupdown1=np.concatenate((z2,dataOIId,z1), axis=0)
                dataOIIdnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O2_7330e' and requested[i]=='no':
                dataOIId=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIId,z1), axis=0)
                dataOIIdnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O1_5577e' and requested[i]=='yes':
                dataOIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIa=np.nan_to_num(dataOIa)
                addupdown1=np.concatenate((z2,dataOIa,z1), axis=0)
                dataOIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O1_5577e' and requested[i]=='no':
                dataOIa=np.zeros((sy, sx))   
                addupdown1=np.concatenate((z2,dataOIa,z1), axis=0)
                dataOIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O1_6300e' and requested[i]=='yes':
                dataOIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIb=np.nan_to_num(dataOIb)
                addupdown1=np.concatenate((z2,dataOIb,z1), axis=0)
                dataOIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O1_6300e' and requested[i]=='no':
                dataOIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIb,z1), axis=0)
                dataOIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'O1_6363e' and requested[i]=='yes':
                dataOIc, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataOIc=np.nan_to_num(dataOIc)
                addupdown1=np.concatenate((z2,dataOIc,z1), axis=0)
                dataOIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'O1_6363e' and requested[i]=='no':
                dataOIc=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataOIc,z1), axis=0)
                dataOIcnew=np.concatenate((z3,addupdown1,z4), axis=1)
                    
            if line_name[i] == 'S2_6716e' and requested[i]=='yes':
                dataSIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataSIIa=np.nan_to_num(dataSIIa)
                addupdown1=np.concatenate((z2,dataSIIa,z1), axis=0)
                dataSIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'S2_6716e' and requested[i]=='no':
                dataSIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIa,z1), axis=0)
                dataSIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'S2_6731e' and requested[i]=='yes':
                dataSIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataSIIb=np.nan_to_num(dataSIIb)
                addupdown1=np.concatenate((z2,dataSIIb,z1), axis=0)
                dataSIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'S2_6731e' and requested[i]=='no':
                dataSIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIb,z1), axis=0)
                dataSIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'S3_6312e' and requested[i]=='yes':
                dataSIIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataSIIIa=np.nan_to_num(dataSIIIa)
                addupdown1=np.concatenate((z2,dataSIIIa,z1), axis=0)
                dataSIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'S3_6312e' and requested[i]=='no':
                dataSIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIIa,z1), axis=0)
                dataSIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'S3_9069e' and requested[i]=='yes':
                dataSIIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataSIIIb=np.nan_to_num(dataSIIIb)
                addupdown1=np.concatenate((z2,dataSIIIb,z1), axis=0)
                dataSIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'S3_9069e' and requested[i]=='no':
                dataSIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataSIIIb,z1), axis=0)
                dataSIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
                
            if line_name[i] == 'Cl3_5517e' and requested[i]=='yes':
                dataClIIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataClIIIa=np.nan_to_num(dataClIIIa)
                addupdown1=np.concatenate((z2,dataClIIIa,z1), axis=0)
                dataClIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Cl3_5517e' and requested[i]=='no':
                dataClIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataClIIIa,z1), axis=0)
                dataClIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Cl3_5538e' and requested[i]=='yes':
                dataClIIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataClIIIb=np.nan_to_num(dataClIIIb)
                addupdown1=np.concatenate((z2,dataClIIIb,z1), axis=0)
                dataClIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Cl3_5538e' and requested[i]=='no':
                dataClIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataClIIIb,z1), axis=0)
                dataClIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Ar3_7136e' and requested[i]=='yes':
                dataArIII, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataArIII=np.nan_to_num(dataArIII)
                addupdown1=np.concatenate((z2,dataArIII,z1), axis=0)
                dataArIIInew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ar3_7136e' and requested[i]=='no':
                dataArIII=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIII,z1), axis=0)
                dataArIIInew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Ar4_4712e' and requested[i]=='yes':
                dataArIVa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataArIVa=np.nan_to_num(dataArIVa)
                addupdown1=np.concatenate((z2,dataArIVa,z1), axis=0)
                dataArIVanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ar4_4712e' and requested[i]=='no':
                dataArIVa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIVa,z1), axis=0)
                dataArIVanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Ar4_4740e' and requested[i]=='yes':
                dataArIVb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataArIVb=np.nan_to_num(dataArIVb)
                addupdown1=np.concatenate((z2,dataArIVb,z1), axis=0)
                dataArIVbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ar4_4740e' and requested[i]=='no':
                dataArIVb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataArIVb,z1), axis=0)
                dataArIVbnew=np.concatenate((z3,addupdown1,z4), axis=1)

            if line_name[i] == 'C1_8727e' and requested[i]=='yes':
                dataCI, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataCI=np.nan_to_num(dataCI)
                addupdown1=np.concatenate((z2,dataCI,z1), axis=0)
                dataCInew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'C1_8727e' and requested[i]=='no':
                dataCI=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataCI,z1), axis=0)
                dataCInew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'C2_6461e' and requested[i]=='yes':
                dataCII, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataCII=np.nan_to_num(dataCII)
                addupdown1=np.concatenate((z2,dataCII,z1), axis=0)
                dataCIInew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'C2_6461e' and requested[i]=='no':
                dataCII=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataCII,z1), axis=0)
                dataCIInew=np.concatenate((z3,addupdown1,z4), axis=1)
    
            if line_name[i] == 'Ne3_3868e' and requested[i]=='yes':
                dataNeIIIa, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNeIIIa=np.nan_to_num(dataNeIIIa)
                addupdown1=np.concatenate((z2,dataNeIIIa,z1), axis=0)
                dataNeIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ne3_3868e' and requested[i]=='no':
                dataNeIIIa=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNeIIIa,z1), axis=0)
                dataNeIIIanew=np.concatenate((z3,addupdown1,z4), axis=1)
            if line_name[i] == 'Ne3_3967e' and requested[i]=='yes':
                dataNeIIIb, hdr = getdata('image_data/'+line_name[i]+'.fit', 0, header=True)
                dataNeIIIb=np.nan_to_num(dataNeIIIb)
                addupdown1=np.concatenate((z2,dataNeIIIb,z1), axis=0)
                dataNeIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
            elif line_name[i] == 'Ne3_3967e' and requested[i]=='no':
                dataNeIIIb=np.zeros((sy, sx))
                addupdown1=np.concatenate((z2,dataNeIIIb,z1), axis=0)
                dataNeIIIbnew=np.concatenate((z3,addupdown1,z4), axis=1)
    
    
    
    
        return dataHanew,dataHbnew,dataHgnew,dataHdnew,dataHeIanew,dataHeIbnew,dataHeIIanew,dataHeIIbnew,dataNIIanew,dataNIIbnew,dataNIIcnew,dataNInew,dataOIIIanew,dataOIIIbnew,dataOIIIcnew,dataOIIanew,dataOIIbnew,dataOIIcnew,dataOIIdnew,dataOIanew,dataOIbnew,dataOIcnew,dataSIIanew,dataSIIbnew,dataSIIIanew,dataSIIIbnew,dataClIIIanew,dataClIIIbnew,dataArIIInew,dataArIVanew,dataArIVbnew,dataCInew,dataCIInew,dataNeIIIanew,dataNeIIIbnew,hdr





