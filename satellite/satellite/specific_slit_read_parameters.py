# specific_slit_read_parameters.py:
# Reads the parameters of the 10 pseudo-slits (width, length, position angle, coordinates in x,y)
# (C) Stavros Akras




def read_input_slit_parameters(par_slit_name,par_slit_value,i):

    if i==1:
        index_line=par_slit_name.index("PA_for_specific_slit_1")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_1")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_1") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_1") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_1") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==2:
        index_line=par_slit_name.index("PA_for_specific_slit_2")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_2")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_2") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_2") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_2") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==3:
        index_line=par_slit_name.index("PA_for_specific_slit_3")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_3")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_3") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_3") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_3") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==4:
        index_line=par_slit_name.index("PA_for_specific_slit_4")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_4")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_4") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_4") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_4") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==5:
        index_line=par_slit_name.index("PA_for_specific_slit_5")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_5")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_5") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_5") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_5") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==6:
        index_line=par_slit_name.index("PA_for_specific_slit_6")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_6")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_6") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_6") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_6") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==7:
        index_line=par_slit_name.index("PA_for_specific_slit_7")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_7")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_7") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_7") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_7") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==8:
        index_line=par_slit_name.index("PA_for_specific_slit_8")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_8")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_8") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_8") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_8") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==9:
        index_line=par_slit_name.index("PA_for_specific_slit_9")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_9")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_9") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_9") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_9") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
    
    if i==10:
        index_line=par_slit_name.index("PA_for_specific_slit_10")        
        PA_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("width_for_specific_slit_10")
        wid_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("length_for_specific_slit_10") 
        len_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("x_coor_of_spec_slit_10") 
        xcoor_slit=par_slit_value[index_line]
        index_line=par_slit_name.index("y_coor_of_spec_slit_10") 
        ycoor_slit=par_slit_value[index_line]

        return PA_slit,wid_slit,len_slit,xcoor_slit,ycoor_slit
