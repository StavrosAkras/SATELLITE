# flux_angles_norm_definition_script.py:
# Defines the flux_angles_norm and flux_angles_norm_error arrays 
# in case the rotation analysis module is desactivated and returns two classes of arrays with zero values.
# (C) Stavros Akras




def flux_ang_norm_definition(param_mod_name,param_model_values):
        
           
    class flux_angles_norm:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]
    
    class flux_angles_norm_error:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]
    
    index_step_ang=param_mod_name.index("position_angle_step")
    index_star_ang=param_mod_name.index("minimum_position_angle")
    index_end_ang=param_mod_name.index("maximum_position_angle")
    start_angle=param_model_values[index_star_ang]
    end_angle=param_model_values[index_end_ang]
    step_angle=param_model_values[index_step_ang]
    
    ang=[]
    kkangle=0
    for i in range(start_angle,end_angle,step_angle):
        ang.append(0)
        ang[kkangle]=i
    
        
        flux_angles_norm.Ha_6563.append(-99)
        flux_angles_norm.Hb_4861.append(-99)
        flux_angles_norm.Hg_4340.append(-99)
        flux_angles_norm.Hd_4101.append(-99)
        flux_angles_norm.HeIa_5876.append(-99)
        flux_angles_norm.HeIb_6678.append(-99)
        flux_angles_norm.HeIIa_4686.append(-99)
        flux_angles_norm.HeIIb_5412.append(-99)
        flux_angles_norm.NIIa_5755.append(-99)
        flux_angles_norm.NIIb_6548.append(-99)
        flux_angles_norm.NIIc_6584.append(-99)
        flux_angles_norm.NI_5199.append(-99)
        flux_angles_norm.OIIIa_4363.append(-99)
        flux_angles_norm.OIIIb_4959.append(-99)
        flux_angles_norm.OIIIc_5007.append(-99)
        flux_angles_norm.OIIa_3727.append(-99)
        flux_angles_norm.OIIb_3729.append(-99)
        flux_angles_norm.OIIc_7320.append(-99)
        flux_angles_norm.OIId_7330.append(-99)
        flux_angles_norm.OIa_5577.append(-99)
        flux_angles_norm.OIb_6300.append(-99)
        flux_angles_norm.OIc_6363.append(-99)
        flux_angles_norm.SIIa_6716.append(-99)
        flux_angles_norm.SIIb_6731.append(-99)
        flux_angles_norm.SIIIa_6312.append(-99)
        flux_angles_norm.SIIIb_9069.append(-99)
        flux_angles_norm.ClIIIa_5517.append(-99)
        flux_angles_norm.ClIIIb_5538.append(-99)
        flux_angles_norm.ArIII_7136.append(-99)
        flux_angles_norm.ArIVa_4712.append(-99)
        flux_angles_norm.ArIVb_4740.append(-99)
        flux_angles_norm.CI_8727.append(-99)
        flux_angles_norm.CII_6461.append(-99)
        flux_angles_norm.NeIIIa_3868.append(-99)
        flux_angles_norm.NeIIIb_3967.append(-99)
           
        flux_angles_norm_error.Ha_6563.append(-99)
        flux_angles_norm_error.Hb_4861.append(-99)
        flux_angles_norm_error.Hg_4340.append(-99)
        flux_angles_norm_error.Hd_4101.append(-99)
        flux_angles_norm_error.HeIa_5876.append(-99)
        flux_angles_norm_error.HeIb_6678.append(-99)
        flux_angles_norm_error.HeIIa_4686.append(-99)
        flux_angles_norm_error.HeIIb_5412.append(-99)
        flux_angles_norm_error.NIIa_5755.append(-99)
        flux_angles_norm_error.NIIb_6548.append(-99)
        flux_angles_norm_error.NIIc_6584.append(-99)
        flux_angles_norm_error.NI_5199.append(-99)
        flux_angles_norm_error.OIIIa_4363.append(-99)
        flux_angles_norm_error.OIIIb_4959.append(-99)
        flux_angles_norm_error.OIIIc_5007.append(-99)
        flux_angles_norm_error.OIIa_3727.append(-99)
        flux_angles_norm_error.OIIb_3729.append(-99)
        flux_angles_norm_error.OIIc_7320.append(-99)
        flux_angles_norm_error.OIId_7330.append(-99)
        flux_angles_norm_error.OIa_5577.append(-99)
        flux_angles_norm_error.OIb_6300.append(-99)
        flux_angles_norm_error.OIc_6363.append(-99)
        flux_angles_norm_error.SIIa_6716.append(-99)
        flux_angles_norm_error.SIIb_6731.append(-99)
        flux_angles_norm_error.SIIIa_6312.append(-99)
        flux_angles_norm_error.SIIIb_9069.append(-99)
        flux_angles_norm_error.ClIIIa_5517.append(-99)
        flux_angles_norm_error.ClIIIb_5538.append(-99)
        flux_angles_norm_error.ArIII_7136.append(-99)
        flux_angles_norm_error.ArIVa_4712.append(-99)
        flux_angles_norm_error.ArIVb_4740.append(-99)
        flux_angles_norm_error.CI_8727.append(-99)
        flux_angles_norm_error.CII_6461.append(-99)
        flux_angles_norm_error.NeIIIa_3868.append(-99)
        flux_angles_norm_error.NeIIIb_3967.append(-99)

    return flux_angles_norm, flux_angles_norm_error,ang


def flux_ang_norm_definition_specificslits(param_mod_name,param_model_values):
        
           
    class flux_angles_norm:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]
    
    class flux_angles_norm_error:
        Ha_6563=[]
        Hb_4861=[]
        Hg_4340=[]
        Hd_4101=[]
        HeIa_5876=[]
        HeIb_6678=[]
        HeIIa_4686=[]
        HeIIb_5412=[]
        NIIa_5755=[]
        NIIb_6548=[]
        NIIc_6584=[]
        NI_5199=[]
        OIIIa_4363=[]
        OIIIb_4959=[]
        OIIIc_5007=[]
        OIIa_3727=[]
        OIIb_3729=[]
        OIIc_7320=[]
        OIId_7330=[]
        OIa_5577=[]
        OIb_6300=[]
        OIc_6363=[]
        SIIa_6716=[]
        SIIb_6731=[]
        SIIIa_6312=[]
        SIIIb_9069=[]
        ClIIIa_5517=[]
        ClIIIb_5538=[]
        ArIII_7136=[]
        ArIVa_4712=[]
        ArIVb_4740=[]
        CI_8727=[]
        CII_6461=[]
        NeIIIa_3868=[]
        NeIIIb_3967=[]
    
    for i in range(0,10):

        
        flux_angles_norm.Ha_6563.append(-99)
        flux_angles_norm.Hb_4861.append(-99)
        flux_angles_norm.Hg_4340.append(-99)
        flux_angles_norm.Hd_4101.append(-99)
        flux_angles_norm.HeIa_5876.append(-99)
        flux_angles_norm.HeIb_6678.append(-99)
        flux_angles_norm.HeIIa_4686.append(-99)
        flux_angles_norm.HeIIb_5412.append(-99)
        flux_angles_norm.NIIa_5755.append(-99)
        flux_angles_norm.NIIb_6548.append(-99)
        flux_angles_norm.NIIc_6584.append(-99)
        flux_angles_norm.NI_5199.append(-99)
        flux_angles_norm.OIIIa_4363.append(-99)
        flux_angles_norm.OIIIb_4959.append(-99)
        flux_angles_norm.OIIIc_5007.append(-99)
        flux_angles_norm.OIIa_3727.append(-99)
        flux_angles_norm.OIIb_3729.append(-99)
        flux_angles_norm.OIIc_7320.append(-99)
        flux_angles_norm.OIId_7330.append(-99)
        flux_angles_norm.OIa_5577.append(-99)
        flux_angles_norm.OIb_6300.append(-99)
        flux_angles_norm.OIc_6363.append(-99)
        flux_angles_norm.SIIa_6716.append(-99)
        flux_angles_norm.SIIb_6731.append(-99)
        flux_angles_norm.SIIIa_6312.append(-99)
        flux_angles_norm.SIIIb_9069.append(-99)
        flux_angles_norm.ClIIIa_5517.append(-99)
        flux_angles_norm.ClIIIb_5538.append(-99)
        flux_angles_norm.ArIII_7136.append(-99)
        flux_angles_norm.ArIVa_4712.append(-99)
        flux_angles_norm.ArIVb_4740.append(-99)
        flux_angles_norm.CI_8727.append(-99)
        flux_angles_norm.CII_6461.append(-99)
        flux_angles_norm.NeIIIa_3868.append(-99)
        flux_angles_norm.NeIIIb_3967.append(-99)
           
        flux_angles_norm_error.Ha_6563.append(-99)
        flux_angles_norm_error.Hb_4861.append(-99)
        flux_angles_norm_error.Hg_4340.append(-99)
        flux_angles_norm_error.Hd_4101.append(-99)
        flux_angles_norm_error.HeIa_5876.append(-99)
        flux_angles_norm_error.HeIb_6678.append(-99)
        flux_angles_norm_error.HeIIa_4686.append(-99)
        flux_angles_norm_error.HeIIb_5412.append(-99)
        flux_angles_norm_error.NIIa_5755.append(-99)
        flux_angles_norm_error.NIIb_6548.append(-99)
        flux_angles_norm_error.NIIc_6584.append(-99)
        flux_angles_norm_error.NI_5199.append(-99)
        flux_angles_norm_error.OIIIa_4363.append(-99)
        flux_angles_norm_error.OIIIb_4959.append(-99)
        flux_angles_norm_error.OIIIc_5007.append(-99)
        flux_angles_norm_error.OIIa_3727.append(-99)
        flux_angles_norm_error.OIIb_3729.append(-99)
        flux_angles_norm_error.OIIc_7320.append(-99)
        flux_angles_norm_error.OIId_7330.append(-99)
        flux_angles_norm_error.OIa_5577.append(-99)
        flux_angles_norm_error.OIb_6300.append(-99)
        flux_angles_norm_error.OIc_6363.append(-99)
        flux_angles_norm_error.SIIa_6716.append(-99)
        flux_angles_norm_error.SIIb_6731.append(-99)
        flux_angles_norm_error.SIIIa_6312.append(-99)
        flux_angles_norm_error.SIIIb_9069.append(-99)
        flux_angles_norm_error.ClIIIa_5517.append(-99)
        flux_angles_norm_error.ClIIIb_5538.append(-99)
        flux_angles_norm_error.ArIII_7136.append(-99)
        flux_angles_norm_error.ArIVa_4712.append(-99)
        flux_angles_norm_error.ArIVb_4740.append(-99)
        flux_angles_norm_error.CI_8727.append(-99)
        flux_angles_norm_error.CII_6461.append(-99)
        flux_angles_norm_error.NeIIIa_3868.append(-99)
        flux_angles_norm_error.NeIIIb_3967.append(-99)

    return flux_angles_norm, flux_angles_norm_error
