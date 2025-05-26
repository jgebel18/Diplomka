import numpy as np
from MakeFunctions import MakeFunctions
from numpy.lib.scimath import sqrt as csqrt
import matplotlib.pyplot as plt
class MakeIterationFuntions:


  #Here is a construction of class variables
    def __init__(self, beta, x_values, t_value, omega_values, U):
        self.beta = beta            # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value #Hopping parameter
        self.omegavalues=omega_values            #omega for Greeen function                # Hopping parameter
        self.U = U   #ElectronIteraction
        self.Mkf= MakeFunctions(self.beta, self.x_values, self.t_value, self.omegavalues,self.U) # This calling of class is important for define another iterations


    #This is a First equation for a
    def a (self):
        a= 1-self.Mkf.Y()
        return a

    #Another Iterations of Gamma
    def Gamma (self):
        DotofItegrals= csqrt(2/(self.Mkf.D()*self.U))
        Gamma= (self.U/(self.beta*np.pi))*DotofItegrals*(1/np.sqrt(self.a))
        return Gamma

    # Here is a Another itertion of Sigma function
    def NewSigma(self, omega):
        Sigma=self.Gamma()*self.Mkf.G(omega)
        return Sigma







