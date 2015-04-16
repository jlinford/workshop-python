#---------------------- BEGIN small_parameters.py BEGIN ----------------------
# @file small_parameters.py                                                   
# @author jlinford                                                            
# @date 2015-04-15 21:17:47.427064                                            
# @brief Program parameters                                                   
#                                                                             
# Integration tolerances, program constants, and species indices.             
#                                                                             
# This file was generated by Kppa: http://www.paratools.com/Kppa              
#-----------------------------------------------------------------------------
import numpy as np
import sys


#-----------------------------------------------------------------------------
# Integration tolerances                                                      
#-----------------------------------------------------------------------------

# Absolute tolerance 
ATOLS           =np.float32(1.0)
# Relative tolerance 
RTOLS           =np.float32(0.001)


#-----------------------------------------------------------------------------
# Concentration constants                                                     
#-----------------------------------------------------------------------------

# Conversion factor 
CFACTOR         =np.float32(1.0)
# Default initialization value for variable species 
VAR_DEFAULT     =np.float32(0.0)
# Default initialization value for fixed species 
FIX_DEFAULT     =np.float32(0.0)


#-----------------------------------------------------------------------------
# Program constants                                                           
#-----------------------------------------------------------------------------

# Species count 
NSPEC           =int(7)
# Variable species count 
NVAR            =int(5)
# Fixed species count 
NFIX            =int(2)
# Active variable species count 
NVARACT         =int(5)
# Reaction (equation) count 
NREACT          =int(10)


#-----------------------------------------------------------------------------
# Numerical constants                                                         
#-----------------------------------------------------------------------------

ZERO            =np.float32(0.0)
HALF            =np.float32(0.5)
ONE             =np.float32(1.0)
TWO             =np.float32(2.0)


#-----------------------------------------------------------------------------
# Variable species indices                                                    
#-----------------------------------------------------------------------------

IND_O1D         =int(0)
IND_O           =int(1)
IND_O3          =int(2)
IND_NO          =int(3)
IND_NO2         =int(4)


#-----------------------------------------------------------------------------
# Fixed species indices                                                       
#-----------------------------------------------------------------------------

IND_M           =int(5)
IND_O2          =int(6)



#------------------------ END small_parameters.py END ------------------------