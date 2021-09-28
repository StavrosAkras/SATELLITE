# read_input_lines_parameters_script.py:
# Reads all the input parameters that the code needs such as the names of the emission line provided by the user (yes or no), 
# the numerical parameters, the outputs that the code will return and the input for the diagnostic diagrams
# (C) Stavros Akras

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import array
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import rotate
from astropy.io import fits
from astropy.io.fits import getdata
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
import matplotlib.font_manager as font_manager
import random
import math




def read_input_lines_parameters():
    name=open("input.txt","r")
    datainput=name.read()
    datainput2 = datainput.split()
    
    
    line_names=[]
    line_provided=[]
    line_radial=[]
    line_ext_error=[]
    
    k1=k2=k3=k4=0
    for i in range(0,len(datainput2),4):
        line_names.append("none")
        line_names[k1]=datainput2[i]
        k1=k1+1
    for i in range(1,len(datainput2),4):
        line_provided.append("none")
        line_provided[k2]=datainput2[i]
        k2=k2+1
    for i in range(2,len(datainput2),4):
        line_radial.append("none")
        line_radial[k3]=datainput2[i]
        k3=k3+1    
    for i in range(3,len(datainput2),4):
        line_ext_error.append(0)
        line_ext_error[k4]=datainput2[i]
        k4=k4+1    
    
    name=open("outputs.txt","r")
    dataoutput=name.read()
    dataoutput2 = dataoutput.split()
    
    parameter_outputs=[]
    parameters_derived=[]
    
    k1=k2=0
    for i in range(0,len(dataoutput2),2):
        parameter_outputs.append("none")
        parameter_outputs[k1]=dataoutput2[i]
        k1=k1+1
    for i in range(1,len(dataoutput2),2):
        parameters_derived.append("none")
        parameters_derived[k2]=dataoutput2[i]
        k2=k2+1

    
    name=open("numerical_input.txt","r")    
    data_num_inp=name.read()
    data_num_inp2 = data_num_inp.split()

    parameter_model=[]
    parameter_values=[]
    
    k1=k2=0
    for i in range(0,len(data_num_inp2),2):
        parameter_model.append("none")
        parameter_model[k1]=data_num_inp2[i]
        k1=k1+1
    for i in range(1,len(data_num_inp2),2):
        parameter_values.append(-10e-10)
        parameter_values[k2]=int(data_num_inp2[i])
        k2=k2+1
        
        
    name=open("diagnostic_diagrams_input.txt","r")    
    data_input_DD=name.read()
    data_input_DD2 = data_input_DD.split()

    parameter_DD=[]
    parameter_values_DD0=[]
    parameter_values_DD1=[]
    parameter_values_DD2=[]
    parameter_values_DD3=[]
    parameter_values_DD4=[]
    
    k1=k2=k3=k4=k5=k6=0
    for i in range(0,len(data_input_DD2),6):
        parameter_DD.append("none")
        parameter_DD[k1]=data_input_DD2[i]
        k1=k1+1
    for i in range(1,len(data_input_DD2),6):
        parameter_values_DD0.append("none")
        parameter_values_DD0[k2]=data_input_DD2[i]
        k2=k2+1
    for i in range(2,len(data_input_DD2),6):
        parameter_values_DD1.append(-10e-10)
        parameter_values_DD1[k3]=data_input_DD2[i]
        k3=k3+1   
    for i in range(3,len(data_input_DD2),6):
        parameter_values_DD2.append(-10e-10)
        parameter_values_DD2[k4]=data_input_DD2[i]
        k4=k4+1   
    for i in range(4,len(data_input_DD2),6):
        parameter_values_DD3.append(-10e-10)
        parameter_values_DD3[k5]=data_input_DD2[i]
        k5=k5+1   
    for i in range(5,len(data_input_DD2),6):
        parameter_values_DD4.append(-10e-10)
        parameter_values_DD4[k6]=data_input_DD2[i]
        k6=k6+1   

    name=open("plots_parameters_input.txt","r")    
    datainputplots=name.read()
    datainputplots2 = datainputplots.split()

    parameter_plotname=[]
    parameter_plotymin=[]
    parameter_plotymax=[]

    
    k1=k2=k3=0
    for i in range(0,len(datainputplots2),3):
        parameter_plotname.append("none")
        parameter_plotname[k1]=datainputplots2[i]
        k1=k1+1
    for i in range(1,len(datainputplots2),3):
        parameter_plotymin.append(0)
        parameter_plotymin[k2]=datainputplots2[i]
        k2=k2+1
    for i in range(2,len(datainputplots2),3):
        parameter_plotymax.append(0)
        parameter_plotymax[k3]=datainputplots2[i]
        k3=k3+1   
    

    return line_names,line_ext_error,line_provided,line_radial,parameter_outputs,parameters_derived,parameter_model,parameter_values,parameter_DD,parameter_values_DD0,parameter_values_DD1,parameter_values_DD2,parameter_values_DD3,parameter_values_DD4,parameter_plotname,parameter_plotymin,parameter_plotymax
