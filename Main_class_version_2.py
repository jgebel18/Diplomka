import numpy as np
import matplotlib
import os

from SolveNonlinearEquation import SolveNonlinearEquation
from CreateIterationProcess import MakeIterationFuntions
from SolveQuarticEquation import SolveQuarticEquation
from MakeFunction_Gamma import MakeFunctions_Gamma
from SolveForGamma import MakeIterationProcessforGamma
matplotlib.use('TkAgg')


def ZerothIterationofSelfEnergy(U,t):
        Sigma= 1j *3* Gamma_0(U,t) / (2 *t)
        return Sigma

#Gamma_0 (Zeroth iteration of Gamma)
def Gamma_0(U,t):
    parameter = (2 * np.pi * t) / U
    gamma_0 = -(2 *t) ** 2 * np.exp(parameter)
    # print(gamma_0)
    return gamma_0

def NewGamma_0(U,t):
   return 1/np.sinh((2*np.pi*t)/(-U))*(-(2*t)**2)

Image_dir= 'Images'
File_dir= 'Files'
# Simulation parameters
t = 1.0  #Hopping parameter
U = -2.0 # Culoumb Interaction
beta = (30000*1/(np.abs(ZerothIterationofSelfEnergy(U,t))))#Inverse Temperature

ScaleFactor=1.5
while True:
    try:

        Rezolution= int(input("Suggest the"
                              " Rezolution:")) # Rezolution Sunils suggest Exeption for the future
        if Rezolution>=0:

            x_values= np.linspace(-(Gamma_0(U,t)/(2*t))*ScaleFactor, (Gamma_0(U,t)/(2*t))*ScaleFactor, Rezolution)
            #x_values = np.linspace(-0.5,0.5, Rezolution)
            omega_values= np.linspace(-10,10
                            , Rezolution)
            NumIteration= 100
            Tolerance= 1e-12
            # Class initialization
            #IP=MakeIterationFuntions(beta,x_values, t,omega_values, U, NumIteration,
             #                         Tolerance, Rezolution,
              #                         Gamma_0(U,t),ZerothIterationofSelfEnergy(U,t))
            IP= MakeIterationProcessforGamma(beta,x_values, t,omega_values, U, NumIteration,
                                      Tolerance, Rezolution,
                                      Gamma_0(U,t),ZerothIterationofSelfEnergy(U,t))
            IP.GiveFinalG()
            break
        else:
            print("Rezolution must be positive.")
    except ValueError:
        print("Wrong Value Try it again.")