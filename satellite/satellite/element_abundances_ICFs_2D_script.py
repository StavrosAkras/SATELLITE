# element_abundances_ICFs_2D_script.py:
# Computes the total elemental abundances and ICFs for all the pixels in 2D maps
# (C) Stavros Akras

import numpy as np
#import pyfits as pf
import pyneb as pn
from numpy import inf


from satellite import norm_flux_error_script as nfes


def total_abundances_ICFs(ion_abun,ion_abun_error,atom_abun_fake,i,j,line_names,lines_available,elem_abun_KB,elem_ICF_KB,elem_abun_DIMS,elem_ICF_DIMS,elem_abun_KB_error,elem_ICF_KB_error,elem_abun_DIMS_error,elem_ICF_DIMS_error,elem_abun_KB_ratio,elem_abun_KB_ratio_error,elem_abun_DIMS_ratio,elem_abun_DIMS_ratio_error,param_mod_name):

    index_atomicData=12
    name_atomicdata=str(param_mod_name[index_atomicData])
    
    if (name_atomicdata!="Chianti"):
        if (name_atomicdata!="IRAF_09"):
            pn.atomicData.includeFitsPath()
            pn.atomicData.includeDeprecatedPath()
            pn.atomicData.setDataFileDict(name_atomicdata)
        else:
            pn.atomicData.includeFitsPath()
            pn.atomicData.setDataFileDict(name_atomicdata)
    else:
        pn.atomicData.addAllChianti()

    ##########################################
    ### Helium
    ##########################################
    index1=index2=0
    index_line1=line_names.index("HeI_5876s")
    if lines_available[index_line1]=="yes": 
        He2_ab_1=ion_abun.HeIa[i,j]
        He2_ab_1_er=ion_abun_error.HeIa[i,j]
        index1=1                                            
                                                    
    index_line1=line_names.index("HeI_6678s")
    if lines_available[index_line1]=="yes": 
        He2_ab_2=ion_abun.HeIb[i,j]
        He2_ab_2_er=ion_abun_error.HeIb[i,j]
        index2=1                                                     

    He2_ab=0
    He2_ab_er=0
    if (index1==1 and index2==1):
         He2_ab=(He2_ab_1+He2_ab_2)/2.
         He2_ab_er=np.sqrt(He2_ab_1_er**2+He2_ab_2_er**2)/2.
    elif (index1==1 and index2==0):
         He2_ab=(He2_ab_1)
         He2_ab_er=(He2_ab_1_er)
    elif (index1==0 and index2==1):
         He2_ab=(He3_ab_2)
         He2_ab_er=(He2_ab_2_er)
    
    index1=index2=0
    index_line1=line_names.index("HeII_4686s")
    if lines_available[index_line1]=="yes": 
        He3_ab_1=ion_abun.HeIIa[i,j]
        He3_ab_1_er=ion_abun_error.HeIIa[i,j]
        index1=1                                            
                                                    
    index_line1=line_names.index("HeII_5412s")
    if lines_available[index_line1]=="yes": 
        He3_ab_2=ion_abun.HeIIb[i,j]
        He3_ab_2_er=ion_abun_error.HeIIb[i,j]
        index2=1  
    

    He3_ab=0
    He3_ab_er=0
    if (index1==1 and index2==1):
         He3_ab=(He3_ab_1+He3_ab_2)/2.
         He3_ab_er=np.sqrt(He3_ab_1_er**2+He3_ab_2_er**2)/2.
    elif (index1==1 and index2==0):
         He3_ab=(He3_ab_1)
         He3_ab_er=(He3_ab_1_er)
    elif (index1==0 and index2==1):
         He3_ab=He3_ab_2
         He3_ab_er=He3_ab_2_er
        

 
    
    ##########################################
    ### OXYGEN
    ##########################################
    index1=index2=index3=index4=0
    index_line1=line_names.index("O2_3727s")
    if lines_available[index_line1]=="yes": 
        O2_ab_1=ion_abun.OIIa[i,j]
        O2_ab_1_er=ion_abun_error.OIIa[i,j]
        index1=1 
    index_line1=line_names.index("O2_3729s")
    if lines_available[index_line1]=="yes": 
        O2_ab_2=ion_abun.OIIb[i,j]
        O2_ab_2_er=ion_abun_error.OIIb[i,j]
        index2=1
    index_line1=line_names.index("O2_7320s")
    if lines_available[index_line1]=="yes": 
        O2_ab_3=ion_abun.OIIc[i,j]
        O2_ab_3_er=ion_abun_error.OIIc[i,j]
        index3=1 
    index_line1=line_names.index("O2_7330s")
    if lines_available[index_line1]=="yes": 
        O2_ab_4=ion_abun.OIId[i,j]
        O2_ab_4_er=ion_abun_error.OIId[i,j]
        index4=1 
    
    O2_ab=0
    O2_ab_er=0
    if (index1==1 and index2==0 and index3==0 and index4==0):
        O2_ab=O2_ab_1
        O2_ab_er=O2_ab_1_er
    elif (index1==0 and index2==1 and index3==0 and index4==0):
        O2_ab=O2_ab_2
        O2_ab_er=O2_ab_2_er
    elif (index1==1 and index2==1 and index3==0 and index4==0):
        O2_ab=(O2_ab_1+O2_ab_2)/2.
        O2_ab_er=np.sqrt(O2_ab_1_er**2+O2_ab_2_er**2)/2.
    elif (index1==0 and index2==0 and index3==1 and index4==0):
        O2_ab=O2_ab_3
        O2_ab_er=O2_ab_3_er
    elif (index1==0 and index2==0 and index3==0 and index4==1):
        O2_ab=O2_ab_4
        O2_ab_er=O2_ab_4_er
    elif (index1==0 and index2==0 and index3==1 and index4==1):
        O2_ab=(O2_ab_3+O2_ab_4)/2.
        O2_ab_er=np.sqrt(O2_ab_3_er**2+O2_ab_4_er**2)/2.
    elif (index1==1 and index2==1 and index3==1 and index4==1):    
        O2_ab=(O2_ab_1+ O2_ab_2+ O2_ab_3+ O2_ab_4)/4.
        O2_ab_er=np.sqrt(O2_ab_1_er**2+O2_ab_2_er**2+O2_ab_3_er**2+O2_ab_4_er**2)/4.
    
        
    index1=index2=index3=index4=0
    index_line1=line_names.index("O3_4363s")
    if lines_available[index_line1]=="yes": 
        O3_ab_1=ion_abun.OIIIa[i,j]
        O3_ab_1_er=ion_abun_error.OIIIa[i,j]
        index1=1 
    index_line1=line_names.index("O3_4959s")
    if lines_available[index_line1]=="yes": 
        O3_ab_2=ion_abun.OIIIb[i,j]
        O3_ab_2_er=ion_abun_error.OIIIb[i,j]
        index2=1
    index_line1=line_names.index("O3_5007s")
    if lines_available[index_line1]=="yes": 
        O3_ab_3=ion_abun.OIIIc[i,j]
        O3_ab_3_er=ion_abun_error.OIIIc[i,j]
        index3=1 

    O3_ab=0 
    O3_ab_er=0
    if (index1==1 and index2==0 and index3==0):
        O3_ab=O3_ab_1
        O3_ab_er=O3_ab_1_er
    elif (index1==0 and index2==1 and index3==0 ):
        O3_ab=O3_ab_2
        O3_ab_er=O3_ab_2_er
    elif (index1==1 and index2==1 and index3==1 ):
        O3_ab=(O3_ab_1+O3_ab_2+O3_ab_3)/3.
        O3_ab_er=np.sqrt(O3_ab_1_er**2+O3_ab_2_er**2+O3_ab_3_er**2)/3.
    elif (index1==0 and index2==1 and index3==1 ):
        O3_ab=(O3_ab_2+O3_ab_3)/2.
        O3_ab_er=np.sqrt(O3_ab_2_er**2+O3_ab_3_er**2)/2.
  
    ##########################################
    ### Nitrogen
    ##########################################
    index1=index2=index3=0
    index_line1=line_names.index("N2_5755s")
    if lines_available[index_line1]=="yes": 
        N2_ab_1=ion_abun.NIIa[i,j]
        N2_ab_1_er=ion_abun_error.NIIa[i,j]
        index1=1 
    index_line1=line_names.index("N2_6548s")
    if lines_available[index_line1]=="yes": 
        N2_ab_2=ion_abun.NIIb[i,j]
        N2_ab_2_er=ion_abun_error.NIIb[i,j]
        index2=1
    index_line1=line_names.index("N2_6583s")
    if lines_available[index_line1]=="yes": 
        N2_ab_3=ion_abun.NIIc[i,j]
        N2_ab_3_er=ion_abun_error.NIIc[i,j]
        index3=1 

    N2_ab=0 
    N2_ab_er=0 
    if (index1==1 and index2==0 and index3==0):
        N2_ab=N2_ab_1
        N2_ab_er=N2_ab_1_er
    elif (index1==0 and index2==1 and index3==0 ):
        N2_ab=N2_ab_2
        N2_ab_er=N2_ab_2_er
    elif (index1==1 and index2==1 and index3==1 ):
        N2_ab=(N2_ab_1+N2_ab_2+N2_ab_3)/3.
        N2_ab_er=np.sqrt(N2_ab_1_er**2+N2_ab_2_er**2+N2_ab_3_er**2)/3.
    elif (index1==0 and index2==1 and index3==1 ):
        N2_ab=(N2_ab_2+N2_ab_3)/2.
        N2_ab_er=np.sqrt(N2_ab_2_er**2+N2_ab_3_er**2)/2.
  
   
    ##########################################
    ### Sulfurt
    ##########################################
    index1=index2=index3=index4=0
    index_line1=line_names.index("S2_6716s")
    if lines_available[index_line1]=="yes": 
        S2_ab_1=ion_abun.SIIa[i,j]
        S2_ab_1_er=ion_abun_error.SIIa[i,j]
        index1=1 
    index_line1=line_names.index("S2_6731s")
    if lines_available[index_line1]=="yes": 
        S2_ab_2=ion_abun.SIIb[i,j]
        S2_ab_2_er=ion_abun_error.SIIb[i,j]
        index2=1
    index_line1=line_names.index("S3_6312s")
    if lines_available[index_line1]=="yes": 
        S3_ab_1=ion_abun.SIIIa[i,j]
        S3_ab_1_er=ion_abun_error.SIIIa[i,j]
        index3=1 
    index_line1=line_names.index("S3_9069s")
    if lines_available[index_line1]=="yes": 
        S3_ab_2=ion_abun.SIIIb[i,j]
        S3_ab_2_er=ion_abun_error.SIIIb[i,j]
        index4=1 
        

    S2_ab=0 
    S2_ab_er=0 
    if (index1==1 and index2==0):
        S2_ab=S2_ab_1
        S2_ab_er=S2_ab_1_er
    elif (index1==0 and index2==1):
        S2_ab=S2_ab_2
        S2_ab_er=S2_ab_2_er
    elif (index1==1 and index2==1):
        S2_ab=(S2_ab_1+S2_ab_2)/2.
        S2_ab_er=np.sqrt(S2_ab_1_er**2+S2_ab_2_er**2)/2.
    
    S3_ab=0 
    S3_ab_er=0
    S3_ab_for_DIMS=-1
    if (index3==1 and index4==0):
        S3_ab=S3_ab_1
        S3_ab_er=S3_ab_1_er
    if (index3==0 and index4==1):
        S3_ab=S3_ab_2
        S3_ab_er=S3_ab_2_er
    elif (index3==1 and index4==1):
        S3_ab=(S3_ab_1+S3_ab_2)/2.
        S3_ab_er=np.sqrt(S3_ab_1_er**2+S3_ab_2_er**2)/2.
        
    if (S3_ab==0 and O2_ab!=0):   
        S3_ab=(4.677+(O3_ab/O2_ab)**0.433)*S2_ab
        a1=a2=a3=0
        a1=(S2_ab*(0.433*(O3_ab)**(-0.567)*(O2_ab)**(-0.433))*O3_ab_er)**2
        a2=(S2_ab*((-0.433)*O3_ab**(0.433)*(O2_ab)**(-1.433))*O2_ab_er)**2
        a3=((4.677+(O3_ab/O2_ab)**0.433)*S2_ab_er)**2
        S3_ab_er=np.sqrt(a1+a2+a3)
        S3_ab_for_DIMS=0
    
    ##########################################
    ### Neon
    ##########################################
    index1=index2=0
    index_line1=line_names.index("Ne3_3868s")
    if lines_available[index_line1]=="yes": 
        Ne3_ab_1=ion_abun.NeIIIa[i,j]
        Ne3_ab_1_er=ion_abun_error.NeIIIa[i,j]
        index1=1 
    index_line1=line_names.index("Ne3_3967s")
    if lines_available[index_line1]=="yes": 
        Ne3_ab_2=ion_abun.NeIIIb[i,j]
        Ne3_ab_2_er=ion_abun_error.NeIIIb[i,j]
        index2=1

    Ne3_ab=0 
    Ne3_ab_er=0 
    if (index1==1 and index2==0):
        Ne3_ab=Ne3_ab_1
        Ne3_ab_er=Ne3_ab_1_er
    elif (index1==0 and index2==1):
        Ne3_ab=Ne3_ab_2
        Ne3_ab_er=Ne3_ab_2_er
    elif (index1==1 and index2==1):
        Ne3_ab=(Ne3_ab_1+Ne3_ab_2)/2.
        Ne3_ab_er=np.sqrt(Ne3_ab_1_er**2+Ne3_ab_2_er**2)/2.
    
    ##########################################
    ### Argon
    ##########################################
    index1=index2=index3=0
    index_line1=line_names.index("Ar3_7136s")
    if lines_available[index_line1]=="yes": 
        Ar3_ab_1=ion_abun.ArIII[i,j]
        Ar3_ab_1_er=ion_abun_error.ArIII[i,j]
        index1=1 
    index_line1=line_names.index("Ar4_4712s")
    if lines_available[index_line1]=="yes": 
        Ar4_ab_1=ion_abun.ArIVa[i,j]
        Ar4_ab_1_er=ion_abun_error.ArIVa[i,j]
        index2=1
    index_line1=line_names.index("Ar4_4740s")
    if lines_available[index_line1]=="yes": 
        Ar4_ab_2=ion_abun.ArIVb[i,j]
        Ar4_ab_2_er=ion_abun_error.ArIVb[i,j]
        index3=1

    Ar3_ab=0 
    Ar3_ab_er=0 
    if (index1==1):
        Ar3_ab=Ar3_ab_1
        Ar3_ab_er=Ar3_ab_1_er
        
    Ar4_ab=0
    Ar4_ab_er=0
    if (index2==1 and index2==0):
        Ar4_ab=Ar4_ab_1
        Ar4_ab_er=Ar4_ab_1_er
    elif (index2==0 and index2==1):
        Ar4_ab=Ar4_ab_2
        Ar4_ab_er=Ar4_ab_2_er
    elif (index1==1 and index2==1):
        Ar4_ab=(Ar4_ab_1+Ar4_ab_2)/2.
        Ar4_ab_er=np.sqrt(Ar4_ab_1_er**2+Ar4_ab_2_er**2)/2.
    
    ##########################################
    ### Chlorine
    ##########################################
    index1=index2=0
    index_line1=line_names.index("Cl3_5517s")
    if lines_available[index_line1]=="yes": 
        Cl3_ab_1=ion_abun.ClIIIa[i,j]
        Cl3_ab_1_er=ion_abun_error.ClIIIa[i,j]
        index1=1 
    index_line1=line_names.index("Cl3_5538s")
    if lines_available[index_line1]=="yes": 
        Cl3_ab_2=ion_abun.ClIIIb[i,j]
        Cl3_ab_2_er=ion_abun_error.ClIIIb[i,j]
        index2=1
    
    Cl2_ab=0
    Cl2_ab_er=0 
    Cl4_ab=0
    Cl4_ab_er=0 
    
    
    Cl3_ab=0
    Cl3_ab_er=0 
    if (index1==1 and index2==0):
        Cl3_ab=Cl3_ab_1
        Cl3_ab_er=Cl3_ab_1
    elif (index2==0 and index2==1):
        Cl3_ab=Cl3_ab_2
        Cl3_ab_er=Cl3_ab_2
    elif (index1==1 and index2==1):
        Cl3_ab=(Cl3_ab_1+Cl3_ab_2)/2.
        Cl3_ab_er=np.sqrt(Cl3_ab_1_er**2+Cl3_ab_2_er**2)/2.
    
    

    if np.isnan(He2_ab)==True:
        He2_ab=0
    if np.isnan(He3_ab)==True:
        He3_ab=0
    if np.isnan(O2_ab)==True:
        O2_ab=0
    if np.isnan(O3_ab)==True:
        O3_ab=0
    if np.isnan(N2_ab)==True:
        N2_ab=0    
    if np.isnan(S2_ab)==True:
        S2_ab=0
    if np.isnan(S3_ab)==True:
        S3_ab=0
    if np.isnan(Ne3_ab)==True:
        Ne3_ab=0
    if np.isnan(Ar3_ab)==True:
        Ar3_ab=0        
    if np.isnan(Ar4_ab)==True:
        Ar4_ab=0    
    if np.isnan(Cl2_ab)==True:
        Cl2_ab=0
    if np.isnan(Cl3_ab)==True:
        Cl3_ab=0
    if np.isnan(Cl4_ab)==True:
        Cl4_ab=0
    
    if S3_ab_for_DIMS==-1:
        atom_abun = {'He2': He2_ab, 'He3':He3_ab,'O2': O2_ab, 'O3': O3_ab, 'N2':N2_ab,'S2': S2_ab, 'S3': S3_ab,'Ne3': Ne3_ab,'Ar3': Ar3_ab,'Ar4': Ar4_ab,'Cl3': Cl3_ab}

    if S3_ab_for_DIMS==0:
        atom_abun = {'He2': He2_ab, 'He3':He3_ab,'O2': O2_ab, 'O3': O3_ab, 'N2':N2_ab,'S2': S2_ab, 'S3': S3_ab,'Ne3': Ne3_ab,'Ar3': Ar3_ab,'Ar4': Ar4_ab,'Cl3': Cl3_ab}
        atom_abun_for_DIMS_S = {'He2': He2_ab, 'He3':He3_ab,'O2': O2_ab, 'O3': O3_ab, 'N2':N2_ab,'S2': S2_ab, 'S3': S3_ab_for_DIMS,'Ne3': Ne3_ab,'Ar3': Ar3_ab,'Ar4': Ar4_ab,'Cl3': Cl3_ab}


    icf = pn.ICF()

    ##########################################
    ### KB ICFS
    ##########################################


    
    elem_abun_KB.He[i,j]=(He2_ab+He3_ab) # Computes the He abundance
    elem_abun_KB_error.He[i,j]=(He2_ab_er+He3_ab_er) # Computes the He abundance

    if (He2_ab!=0):
        if (O2_ab!=0 or O3_ab!=0):
            elem_abun_KB.O[i,j] = icf.getElemAbundance(atom_abun, icf_list=['KB94_A10']).get("KB94_A10") # Computes the O abundance
            abunda=icf.getElemAbundance(atom_abun_fake, icf_list=['KB94_A10']).get("KB94_A10") # Computes the O abundance
            abunda[abunda == inf] =  0
            abunda=np.nan_to_num(abunda)
            abunda = [i2 for i2 in abunda if i2 != 0]
            abunda = [i2 for i2 in abunda if i2 < elem_abun_KB.O[i,j]*3 and i2>elem_abun_KB.O[i,j]*0.3]
            elem_abun_KB_error.O[i,j] = np.std(abunda)
            elem_ICF_KB.O[i,j]=elem_abun_KB.O[i,j]/(O2_ab+O3_ab)
            elem_ICF_KB_error.O[i,j]=elem_ICF_KB.O[i,j]*nfes.ratio_Error(elem_abun_KB.O[i,j],elem_abun_KB_error.O[i,j],O2_ab+O3_ab,O2_ab_er+O3_ab_er)
            
    if (elem_abun_KB.O[i,j]!=0 and O2_ab!=0 and N2_ab!=0):
        elem_abun_KB.N[i,j] = icf.getElemAbundance(atom_abun, icf_list=['KB94_A10','KB94_A1.10']).get('KB94_A1.10') # Computes the N abundance  
        abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['KB94_A10','KB94_A1.10']).get('KB94_A1.10') # Computes the N abundance 
        abunda[abunda == inf] =  0
        abunda=np.nan_to_num(abunda)
        abunda = [i2 for i2 in abunda if i2 != 0]
        abunda = [i2 for i2 in abunda if i2 < elem_abun_KB.N[i,j]*3 and i2>elem_abun_KB.N[i,j]*0.3]
        elem_abun_KB_error.N[i,j] = np.std(abunda)
        elem_ICF_KB.N[i,j]=elem_abun_KB.N[i,j]/N2_ab
        elem_ICF_KB_error.N[i,j]=elem_ICF_KB.N[i,j]*nfes.ratio_Error(elem_abun_KB.N[i,j],elem_abun_KB_error.N[i,j],N2_ab,N2_ab_er)

    
    if (elem_abun_KB.O[i,j]!=0 and O2_ab!=0 and O3_ab!=0 and S2_ab!=0 and S3_ab!=0):
        elem_abun_KB.S[i,j] = icf.getElemAbundance(atom_abun, icf_list=['KB94_A10','KB94_A36.10']).get("KB94_A36.10") # Computes the S abundance
        abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['KB94_A10','KB94_A36.10']).get("KB94_A36.10") # Computes the S abundance
        abunda[abunda == inf] =  0
        abunda=np.nan_to_num(abunda)
        abunda = [i2 for i2 in abunda if i2 != 0]
        abunda = [i2 for i2 in abunda if i2 < elem_abun_KB.S[i,j]*3 and i2>elem_abun_KB.S[i,j]*0.3]
        elem_abun_KB_error.S[i,j] = np.std(abunda)
        elem_ICF_KB.S[i,j]=elem_abun_KB.S[i,j]/(S2_ab+S3_ab)
        elem_ICF_KB_error.S[i,j]=elem_ICF_KB.S[i,j]*nfes.ratio_Error(elem_abun_KB.S[i,j],elem_abun_KB_error.S[i,j],S2_ab+S3_ab,S2_ab_er+S3_ab_er)

        
    if (elem_abun_KB.O[i,j]!=0 and O2_ab!=0 and O3_ab!=0 and Ne3_ab!=0):          
        elem_abun_KB.Ne[i,j] = icf.getElemAbundance(atom_abun, icf_list=['KB94_A10','KB94_A28.10']).get("KB94_A28.10") # Computes the Ne abundance
        abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['KB94_A10','KB94_A28.10']).get("KB94_A28.10") # Computes the Ne abundance
        abunda[abunda == inf] =  0
        abunda=np.nan_to_num(abunda)
        abunda = [i2 for i2 in abunda if i2 != 0]
        abunda = [i2 for i2 in abunda if i2 < elem_abun_KB.Ne[i,j]*3 and i2>elem_abun_KB.Ne[i,j]*0.3]
        elem_abun_KB_error.Ne[i,j] = np.std(abunda)
        elem_ICF_KB.Ne[i,j]=elem_abun_KB.Ne[i,j]/(Ne3_ab)
        elem_ICF_KB_error.Ne[i,j]=elem_ICF_KB.Ne[i,j]*nfes.ratio_Error(elem_abun_KB.Ne[i,j],elem_abun_KB_error.Ne[i,j],Ne3_ab,Ne3_ab_er)


    if (Ar3_ab==0 and Ar4_ab!=0):
        if (elem_abun_KB.Ne[i,j]!=0 and Ne3_ab!=0):
            elem_abun_KB.Ar[i,j] = icf.getElemAbundance(atom_abun, icf_list=['KB94_A10','KB94_A34']).get("KB94_A34") # Computes the Ar abundance 
            abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['KB94_A10','KB94_A34']).get("KB94_A34") # Computes the Ar abundance 
            abunda[abunda == inf] =  0
            abunda=np.nan_to_num(abunda)
            abunda = [i2 for i2 in abunda if i2 != 0]
            abunda = [i2 for i2 in abunda if i2 < elem_abun_KB.Ar[i,j]*3 and i2>elem_abun_KB.Ar[i,j]*0.3]
            elem_abun_KB_error.Ar[i,j] = np.std(abunda)
            elem_ICF_KB.Ar[i,j]=elem_abun_KB.Ar[i,j]/(Ar4_ab)
            elem_ICF_KB_error.Ar[i,j]=elem_ICF_KB.Ar[i,j]*nfes.ratio_Error(elem_abun_KB.Ar[i,j],elem_abun_KB_error.Ar[i,j],Ar4_ab,Ar4_ab_er)
    if (Ar3_ab!=0 and Ar4_ab==0):
            elem_abun_KB.Ar[i,j] = icf.getElemAbundance(atom_abun, icf_list=['KB94_A10','KB94_A32']).get("KB94_A32") # Computes the Ar abundance 
            abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['KB94_A10','KB94_A32']).get("KB94_A32") # Computes the Ar abundance
            abunda[abunda == inf] =  0
            abunda=np.nan_to_num(abunda)
            abunda = [i2 for i2 in abunda if i2 != 0]
            abunda = [i2 for i2 in abunda if i2 < elem_abun_KB.Ar[i,j]*3 and i2>elem_abun_KB.Ar[i,j]*0.3]
            elem_abun_KB_error.Ar[i,j] = np.std(abunda)
            elem_ICF_KB.Ar[i,j]=elem_abun_KB.Ar[i,j]/(Ar3_ab)
            elem_ICF_KB_error.Ar[i,j]=0.0
   
   
    if elem_abun_KB.N[i,j]!=0 and elem_abun_KB.O[i,j]!=0:
        elem_abun_KB_ratio.NO[i,j]=elem_abun_KB.N[i,j]/elem_abun_KB.O[i,j]
        elem_abun_KB_ratio_error.NO[i,j]=elem_abun_KB_ratio.NO[i,j]*np.sqrt((elem_abun_KB_error.N[i,j]/elem_abun_KB.N[i,j])**2 + (elem_abun_KB_error.O[i,j]/elem_abun_KB.O[i,j])**2)
    if elem_abun_KB.S[i,j]!=0 and elem_abun_KB.O[i,j]!=0:
        elem_abun_KB_ratio.SO[i,j]=elem_abun_KB.S[i,j]/elem_abun_KB.O[i,j]
        elem_abun_KB_ratio_error.SO[i,j]=elem_abun_KB_ratio.SO[i,j]*np.sqrt((elem_abun_KB_error.S[i,j]/elem_abun_KB.S[i,j])**2 + (elem_abun_KB_error.O[i,j]/elem_abun_KB.O[i,j])**2)
    if elem_abun_KB.Ne[i,j]!=0 and elem_abun_KB.O[i,j]!=0:
        elem_abun_KB_ratio.NeO[i,j]=elem_abun_KB.Ne[i,j]/elem_abun_KB.O[i,j]
        elem_abun_KB_ratio_error.NeO[i,j]=elem_abun_KB_ratio.NeO[i,j]*np.sqrt((elem_abun_KB_error.Ne[i,j]/elem_abun_KB.Ne[i,j])**2 + (elem_abun_KB_error.O[i,j]/elem_abun_KB.O[i,j])**2)
    if elem_abun_KB.Cl[i,j]!=0 and elem_abun_KB.O[i,j]!=0:
        elem_abun_KB_ratio.ClO[i,j]=elem_abun_KB.Cl[i,j]/elem_abun_KB.O[i,j]
        elem_abun_KB_ratio_error.ClO[i,j]=elem_abun_KB_ratio.ClO[i,j]*np.sqrt((elem_abun_KB_error.Cl[i,j]/elem_abun_KB.Cl[i,j])**2 + (elem_abun_KB_error.O[i,j]/elem_abun_KB.O[i,j])**2)
    if elem_abun_KB.Ar[i,j]!=0 and elem_abun_KB.O[i,j]!=0:
        elem_abun_KB_ratio.ArO[i,j]=elem_abun_KB.Ar[i,j]/elem_abun_KB.O[i,j]
        elem_abun_KB_ratio_error.ArO[i,j]=elem_abun_KB_ratio.ArO[i,j]*np.sqrt((elem_abun_KB_error.Ar[i,j]/elem_abun_KB.Ar[i,j])**2 + (elem_abun_KB_error.O[i,j]/elem_abun_KB.O[i,j])**2)

    ##########################################
    ### DIMS ICFS
    ##########################################
    
    
        
    elem_abun_DIMS.He[i,j]=(He2_ab+He3_ab) # Computes the He abundance
    elem_abun_DIMS_error.He[i,j]=(He2_ab_er+He3_ab_er) # Computes the He abundance
    
    if ((He2_ab+He3_ab)!=0 and (O2_ab+O3_ab)!=0):
        if (O2_ab!=0 or O3_ab!=0):
            if((He3_ab/(He2_ab+He3_ab))<0.95):
                elem_abun_DIMS.O[i,j]=icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12']).get("DIMS14_12") # Computes the O abundance
                abunda =icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12']).get("DIMS14_12") # Computes the O abundance
                abunda[abunda == inf] =  0
                abunda=np.nan_to_num(abunda)
                abunda = [i2 for i2 in abunda if i2 != 0]
                abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.O[i,j]*3 and i2>elem_abun_DIMS.O[i,j]*0.3]
                elem_abun_DIMS_error.O[i,j]= np.std(abunda)
                elem_ICF_DIMS.O[i,j]=elem_abun_DIMS.O[i,j]/(O2_ab+O3_ab)
                elem_ICF_DIMS_error.O[i,j]=elem_ICF_DIMS.O[i,j]*nfes.ratio_Error(elem_abun_DIMS.O[i,j],elem_abun_DIMS_error.O[i,j],O2_ab+O3_ab,O2_ab_er+O3_ab_er)

    
    if (elem_abun_DIMS.O[i,j]!=0 and N2_ab!=0):
        if ((He2_ab+He3_ab)!=0 and (O2_ab+O3_ab)!=0):
            if (O2_ab!=0 or O3_ab!=0): 
                if (He3_ab!=0 and (O3_ab/(O2_ab+O3_ab))<0.96):
                    elem_abun_DIMS.N[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_14']).get("DIMS14_14") # Computes the N abundance
                    abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_14']).get("DIMS14_14") # Computes the N abundance
                    abunda[abunda == inf] =  0
                    abunda=np.nan_to_num(abunda)
                    abunda = [i2 for i2 in abunda if i2 != 0]
                    abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.N[i,j]*3 and i2>elem_abun_DIMS.N[i,j]*0.3]
                    elem_abun_DIMS_error.N[i,j]= np.std(abunda)
                    elem_ICF_DIMS.N[i,j]=elem_abun_DIMS.N[i,j]/(N2_ab)
                    elem_ICF_DIMS_error.N[i,j]=elem_ICF_DIMS.N[i,j]*nfes.ratio_Error(elem_abun_DIMS.N[i,j],elem_abun_DIMS_error.N[i,j],N2_ab,N2_ab_er)
                            
                if (He3_ab==0):
                    elem_abun_DIMS.N[i,j] = icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_14b']).get("DIMS14_14b") # Computes the N abundance
                    abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_14b']).get("DIMS14_14b") # Computes the N abundance
                    abunda[abunda == inf] =  0
                    abunda=np.nan_to_num(abunda)
                    abunda = [i2 for i2 in abunda if i2 != 0]
                    abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.N[i,j]*3 and i2>elem_abun_DIMS.N[i,j]*0.3]
                    elem_abun_DIMS_error.N[i,j]= np.std(abunda)
                    elem_ICF_DIMS.N[i,j]=elem_abun_DIMS.N[i,j]/(N2_ab) 
                    elem_ICF_DIMS_error.N[i,j]=elem_ICF_DIMS.N[i,j]*nfes.ratio_Error(elem_abun_DIMS.N[i,j],elem_abun_DIMS_error.N[i,j],N2_ab,N2_ab_er)

     
    if (elem_abun_DIMS.O[i,j]!=0):
        if ((He2_ab+He3_ab)!=0 and (O2_ab+O3_ab)!=0):
            if (O2_ab!=0 or O3_ab!=0):
                if (S2_ab!=0 and S3_ab_for_DIMS==-1 and (O3_ab/(O2_ab+O3_ab))<0.95):
                    elem_abun_DIMS.S[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_26']).get("DIMS14_26") # Computes the S abundance
                    abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_26']).get("DIMS14_26") # Computes the S abundance
                    abunda[abunda == inf] =  0
                    abunda=np.nan_to_num(abunda)
                    abunda = [i2 for i2 in abunda if i2 != 0]
                    abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.S[i,j]*3 and i2>elem_abun_DIMS.S[i,j]*0.3]
                    elem_abun_DIMS_error.S[i,j]= np.std(abunda)
                    elem_ICF_DIMS.S[i,j]=elem_abun_DIMS.S[i,j]/(S2_ab+S3_ab)
                    elem_ICF_DIMS_error.S[i,j]=elem_ICF_DIMS.S[i,j]*nfes.ratio_Error(elem_abun_DIMS.S[i,j],elem_abun_DIMS_error.S[i,j],S2_ab+S3_ab,S2_ab_er+S3_ab_er)

                if (S2_ab!=0 and S3_ab_for_DIMS==0 and (He3_ab/(He2_ab+He3_ab))>0.02):
                    elem_abun_DIMS.S[i,j]= icf.getElemAbundance(atom_abun_for_DIMS_S, icf_list=['DIMS14_12','DIMS14_23']).get("DIMS14_23") # Computes the S abundance
                    abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_23']).get("DIMS14_23") # Computes the S abundance
                    abunda[abunda == inf] =  0
                    abunda=np.nan_to_num(abunda)
                    abunda = [i2 for i2 in abunda if i2 != 0]
                    abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.S[i,j]*3 and i2>elem_abun_DIMS.S[i,j]*0.3]
                    elem_abun_DIMS_error.S[i,j]= np.std(abunda)
                    elem_ICF_DIMS.S[i,j]=elem_abun_DIMS.S[i,j]/(S2_ab)
                    elem_ICF_DIMS_error.S[i,j]=elem_ICF_DIMS.S[i,j]*nfes.ratio_Error(elem_abun_DIMS.S[i,j],elem_abun_DIMS_error.S[i,j],S2_ab,S2_ab_er)


    if (elem_abun_DIMS.O[i,j]!=0):
        if (O2_ab!=0 or O3_ab!=0): 
            if (Ne3_ab!=0 and (O3_ab/(O2_ab+O3_ab))>0.1): 
                if (He3_ab==0 and (O2_ab+O3_ab)!=0):
                    elem_abun_DIMS.Ne[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_17b']).get("DIMS14_17b") # Computes the Ne abundance
                    abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_17b']).get("DIMS14_17b") # Computes the Ne abundance
                    abunda[abunda == inf] =  0
                    abunda=np.nan_to_num(abunda)
                    abunda = [i2 for i2 in abunda if i2 != 0]
                    abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.Ne[i,j]*3 and i2>elem_abun_DIMS.Ne[i,j]*0.3]
                    elem_abun_DIMS_error.Ne[i,j]= np.std(abunda)
                    elem_ICF_DIMS.Ne[i,j]=elem_abun_DIMS.Ne[i,j]/(Ne3_ab) 
                    elem_ICF_DIMS_error.Ne[i,j]=elem_ICF_DIMS.Ne[i,j]*nfes.ratio_Error(elem_abun_DIMS.Ne[i,j],elem_abun_DIMS_error.Ne[i,j],Ne3_ab,Ne3_ab_er)
                if (He3_ab!=0 and (O2_ab+O3_ab)!=0 and He3_ab/(He2_ab+He3_ab)<0.015):
                    elem_abun_DIMS.Ne[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_17c']).get("DIMS14_17c") # Computes the Ne abundance
                    abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_17c']).get("DIMS14_17c") # Computes the Ne abundance
                    abunda[abunda == inf] =  0
                    abunda=np.nan_to_num(abunda)
                    abunda = [i2 for i2 in abunda if i2 != 0]
                    abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.Ne[i,j]*3 and i2>elem_abun_DIMS.Ne[i,j]*0.3]
                    elem_abun_DIMS_error.Ne[i,j]= np.std(abunda)
                    elem_ICF_DIMS.Ne[i,j]=elem_abun_DIMS.Ne[i,j]/(Ne3_ab) 
                    elem_ICF_DIMS_error.Ne[i,j]=elem_ICF_DIMS.Ne[i,j]*nfes.ratio_Error(elem_abun_DIMS.Ne[i,j],elem_abun_DIMS_error.Ne[i,j],Ne3_ab,Ne3_ab_er)
                if (He3_ab!=0 and (O2_ab+O3_ab)!=0 and He3_ab/(He2_ab+He3_ab)>0.015):
                    elem_abun_DIMS.Ne[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_17a']).get("DIMS14_17a") # Computes the Ne abundance
                    abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_17a']).get("DIMS14_17a") # Computes the Ne abundance
                    abunda[abunda == inf] =  0
                    abunda=np.nan_to_num(abunda)
                    abunda = [i2 for i2 in abunda if i2 != 0]
                    abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.Ne[i,j]*3 and i2>elem_abun_DIMS.Ne[i,j]*0.3]
                    elem_abun_DIMS_error.Ne[i,j]= np.std(abunda)
                    elem_ICF_DIMS.Ne[i,j]=elem_abun_DIMS.Ne[i,j]/(Ne3_ab) 
                    elem_ICF_DIMS_error.Ne[i,j]=elem_ICF_DIMS.Ne[i,j]*nfes.ratio_Error(elem_abun_DIMS.Ne[i,j],elem_abun_DIMS_error.Ne[i,j],Ne3_ab,Ne3_ab_er)

    
    if (elem_abun_DIMS.O[i,j]!=0):
        if ((He2_ab+He3_ab)!=0 and (O2_ab+O3_ab)!=0):
            if (O2_ab!=0 or O3_ab!=0): 
                if (Ar3_ab!=0):
                    if ((O3_ab/(O2_ab+O3_ab))<=0.5 and (O3_ab/(O2_ab+O3_ab))>0.05):
                        elem_abun_DIMS.Ar[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_35']).get("DIMS14_35") # Computes the Ar abundance 
                        abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_35']).get("DIMS14_35") # Computes the Ar abundance 
                        abunda[abunda == inf] =  0
                        abunda=np.nan_to_num(abunda)
                        abunda = [i2 for i2 in abunda if i2 != 0]
                        abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.Ar[i,j]*3 and i2>elem_abun_DIMS.Ar[i,j]*0.3]
                        elem_abun_DIMS_error.Ar[i,j]= np.std(abunda)
                        elem_ICF_DIMS.Ar[i,j]=elem_abun_DIMS.Ar[i,j]/(Ar3_ab)
                        elem_ICF_DIMS_error.Ar[i,j]=elem_ICF_DIMS.Ar[i,j]*nfes.ratio_Error(elem_abun_DIMS.Ar[i,j],elem_abun_DIMS_error.Ar[i,j],Ar3_ab,Ar3_ab_er)

                    if ((O3_ab/(O2_ab+O3_ab))<0.95 and (O3_ab/(O2_ab+O3_ab))>0.5):
                        elem_abun_DIMS.Ar[i,j]=  icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_36']).get("DIMS14_36") # Computes the Ar abundance
                        abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_36']).get("DIMS14_36") # Computes the Ar abundance
                        abunda[abunda == inf] =  0
                        abunda=np.nan_to_num(abunda)
                        abunda = [i2 for i2 in abunda if i2 != 0]
                        abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.Ar[i,j]*3 and i2>elem_abun_DIMS.Ar[i,j]*0.3]
                        elem_abun_DIMS_error.Ar[i,j]= np.std(abunda)
                        elem_ICF_DIMS.Ar[i,j]=elem_abun_DIMS.Ar[i,j]/(Ar3_ab)
                        elem_ICF_DIMS_error.Ar[i,j]=elem_ICF_DIMS.Ar[i,j]*nfes.ratio_Error(elem_abun_DIMS.Ar[i,j],elem_abun_DIMS_error.Ar[i,j],Ar3_ab,Ar3_ab_er)



    if (elem_abun_DIMS.O[i,j]!=0):
        if ((He2_ab+He3_ab)!=0 and (O2_ab+O3_ab)!=0):
            if (O2_ab!=0 or O3_ab!=0): 
                if (Cl3_ab!=0 and Cl2_ab==0 and Cl4_ab==0 ):
                    if ((He3_ab/(He2_ab+He3_ab))>0.02 and (He3_ab/(He2_ab+He3_ab))<0.95):
                        elem_abun_DIMS.Cl[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_29']).get("DIMS14_29") # Computes the Ar abundance 
                        abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_29']).get("DIMS14_29") # Computes the Ar abundance 
                        abunda[abunda == inf] =  0
                        abunda=np.nan_to_num(abunda)
                        abunda = [i2 for i2 in abunda if i2 != 0]
                        abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.Cl[i,j]*3 and i2>elem_abun_DIMS.Cl[i,j]*0.3]
                        elem_abun_DIMS_error.Cl[i,j]= np.std(abunda)
                        elem_ICF_DIMS.Cl[i,j]=elem_abun_DIMS.Cl[i,j]/(Cl3_ab)
                        elem_ICF_DIMS_error.Cl[i,j]=elem_ICF_DIMS.Cl[i,j]*nfes.ratio_Error(elem_abun_DIMS.Cl[i,j],elem_abun_DIMS_error.Cl[i,j],Cl3_ab,Cl3_ab_er)
                if (Cl3_ab!=0 and Cl2_ab!=0 and Cl4_ab==0):
                    if ((He3_ab/(He2_ab+He3_ab))<=0.02):
                        elem_abun_DIMS.Cl[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_29b']).get("DIMS14_29b") # Computes the Ar abundance
                        abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_29b']).get("DIMS14_29b") # Computes the Ar abundance
                        abunda[abunda == inf] =  0
                        abunda=np.nan_to_num(abunda)
                        abunda = [i2 for i2 in abunda if i2 != 0]
                        abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.Cl[i,j]*3 and i2>elem_abun_DIMS.Cl[i,j]*0.3]
                        elem_abun_DIMS_error.Cl[i,j]= np.std(abunda)
                        elem_ICF_DIMS.Cl[i,j]=elem_abun_DIMS.Cl[i,j]/(Cl3_ab+Cl2_ab)
                        elem_ICF_DIMS_error.Cl[i,j]=elem_ICF_DIMS.Cl[i,j]*nfes.ratio_Error(elem_abun_DIMS.Cl[i,j],elem_abun_DIMS_error.Cl[i,j],Cl3_ab+Cl2_ab,Cl3_ab_er+Cl2_ab_er)
                if (Cl3_ab!=0 and Cl2_ab!=0 and Cl4_ab!=0):
                    if ((He3_ab/(He2_ab+He3_ab))>=0 and (He3_ab/(He2_ab+He3_ab))<=1):
                        elem_abun_DIMS.Cl[i,j]= icf.getElemAbundance(atom_abun, icf_list=['DIMS14_12','DIMS14_32']).get("DIMS14_32") # Computes the Ar abundance
                        abunda = icf.getElemAbundance(atom_abun_fake, icf_list=['DIMS14_12','DIMS14_32']).get("DIMS14_32") # Computes the Ar abundance
                        abunda[abunda == inf] =  0
                        abunda=np.nan_to_num(abunda)
                        abunda = [i2 for i2 in abunda if i2 != 0]
                        abunda = [i2 for i2 in abunda if i2 < elem_abun_DIMS.Cl[i,j]*3 and i2>elem_abun_DIMS.Cl[i,j]*0.3]
                        elem_abun_DIMS_error.Cl[i,j]= np.std(abunda)
                        elem_ICF_DIMS.Cl[i,j]=elem_abun_DIMS.Cl[i,j]/(Cl3_ab+Cl2_ab+Cl4_ab)
                        elem_ICF_DIMS_error.Cl[i,j]=elem_ICF_DIMS.Cl[i,j]*nfes.ratio_Error(elem_abun_DIMS.Cl[i,j],elem_abun_DIMS_error.Cl[i,j],Cl3_ab+Cl2_ab+Cl4_ab,Cl3_ab_er+Cl2_ab_er+Cl4_ab_er)
    
    if elem_abun_DIMS.N[i,j]!=0 and elem_abun_DIMS.O[i,j]!=0:
        elem_abun_DIMS_ratio.NO[i,j]=elem_abun_DIMS.N[i,j]/elem_abun_DIMS.O[i,j]
        elem_abun_DIMS_ratio_error.NO[i,j]=elem_abun_DIMS_ratio.NO[i,j]*np.sqrt((elem_abun_DIMS_error.N[i,j]/elem_abun_DIMS.N[i,j])**2 + (elem_abun_DIMS_error.O[i,j]/elem_abun_DIMS.O[i,j])**2)
    if elem_abun_DIMS.S[i,j]!=0 and elem_abun_DIMS.O[i,j]!=0:
        elem_abun_DIMS_ratio.SO[i,j]=elem_abun_DIMS.S[i,j]/elem_abun_DIMS.O[i,j]
        elem_abun_DIMS_ratio_error.SO[i,j]=elem_abun_DIMS_ratio.SO[i,j]*np.sqrt((elem_abun_DIMS_error.S[i,j]/elem_abun_DIMS.S[i,j])**2 + (elem_abun_DIMS_error.O[i,j]/elem_abun_DIMS.O[i,j])**2)
    if elem_abun_DIMS.Ne[i,j]!=0 and elem_abun_DIMS.O[i,j]!=0:
        elem_abun_DIMS_ratio.NeO[i,j]=elem_abun_DIMS.Ne[i,j]/elem_abun_DIMS.O[i,j]
        elem_abun_DIMS_ratio_error.NeO[i,j]=elem_abun_DIMS_ratio.NeO[i,j]*np.sqrt((elem_abun_DIMS_error.Ne[i,j]/elem_abun_DIMS.Ne[i,j])**2 + (elem_abun_DIMS_error.O[i,j]/elem_abun_DIMS.O[i,j])**2)
    if elem_abun_DIMS.Cl[i,j]!=0 and elem_abun_DIMS.O[i,j]!=0:
        elem_abun_DIMS_ratio.ClO[i,j]=elem_abun_DIMS.Cl[i,j]/elem_abun_DIMS.O[i,j]
        elem_abun_DIMS_ratio_error.ClO[i,j]=elem_abun_DIMS_ratio. ClO[i,j]*np.sqrt((elem_abun_DIMS_error.Cl[i,j]/elem_abun_DIMS.Cl[i,j])**2 + (elem_abun_DIMS_error.O[i,j]/elem_abun_DIMS.O[i,j])**2)
    if elem_abun_DIMS.Ar[i,j]!=0 and elem_abun_DIMS.O[i,j]!=0:
        elem_abun_DIMS_ratio.ArO[i,j]=elem_abun_DIMS.Ar[i,j]/elem_abun_DIMS.O[i,j]
        elem_abun_DIMS_ratio_error.ArO[i,j]=elem_abun_DIMS_ratio.ArO[i,j]*np.sqrt((elem_abun_DIMS_error.Ar[i,j]/elem_abun_DIMS.Ar[i,j])**2 + (elem_abun_DIMS_error.O[i,j]/elem_abun_DIMS.O[i,j])**2)
    
    return elem_abun_KB,elem_ICF_KB,elem_abun_DIMS,elem_ICF_DIMS,elem_abun_KB_error,elem_ICF_KB_error,elem_abun_DIMS_error,elem_ICF_DIMS_error,elem_abun_KB_ratio,elem_abun_KB_ratio_error,elem_abun_DIMS_ratio,elem_abun_DIMS_ratio_error


        
