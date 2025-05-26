import numpy as np
import matplotlib.pyplot as plt
from MakeFunctions import MakeFunctions
import matplotlib
from CreateIterationProcess import MakeIterationFuntions
matplotlib.use('TkAgg')
# Simulation parameters
t = 1.0
U = -0.5
beta = 500
x_values = np.linspace(-5, 5, 100)
omega_values= np.linspace(-5, 5, 100)
# Class initialization
MkF = MakeFunctions(beta, x_values, t,omega_values ,U)
MkIF=MakeIterationFuntions(beta,x_values, t,omega_values,U)

D=((1/np.pi)*MkF.FermiFunction(x_values)*MkF.SumOfComplexRations(x_values))
# Calculate the Green function
G=MkF.G(omega_values)

#print()
# Calculate the function a = 1 - U * Y(0,0)
a = 1 - U * MkF.Y()[0]  # [0] is the result of integration, [1] is the error estimate
print("a = ", a)
# Plotting
MkF.PlotData('D function', x_values, D, 'x' ,'D', 'D function plotting' )
MkF.PlotData('Green Function Real Part',omega_values, D, 'ω', 'G(ω)' , 'Green function plotting' )



