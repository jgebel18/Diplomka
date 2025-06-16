import numpy as np
from MakeFunctions import MakeFunctions
from numpy.lib.scimath import sqrt as csqrt
import matplotlib.pyplot as plt
class MakeIterationFuntions:


  #Here is a construction of class variables
    def __init__(self, beta, x_values, t_value, omega_values, U, NumIterations, Tolerance):
        self.beta = beta            # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value #Hopping parameter
        self.omegavalues=omega_values            #omega for Greeen function                # Hopping parameter
        self.U = U   #ElectronIteraction
        self.Mkf= MakeFunctions(self.beta, self.x_values, self.t_value, self.omegavalues,self.U) # This calling of class is important for define another iterations
        self.NumIteration=NumIterations
        self. Tolerance=Tolerance


    #This is a First equation for a
    def a (self, Sigma):
        #print(self.Mkf.Y(Sigma))
        a= 1.0+self.Mkf.Y(Sigma)
        return a

    #Another Iterations of Gamma
    def Gamma (self, Sigma):
        DotofItegrals= csqrt(2/(self.Mkf.D(Sigma)*self.U))
        Gamma= (self.U/(self.beta*np.pi))*DotofItegrals*(1/np.sqrt(self.a(Sigma)))
        return Gamma

    # Here is a Another itertions of Sigma function
    def NewSigma(self, omega, Sigma):
        NewSigma=self.Gamma(Sigma)*self.Mkf.G(omega, Sigma)
        self.Mkf.PlotData("i-tá iteraci", self.omegavalues,self.Mkf.G(omega, Sigma), "ω", "G(ω)", "Greenova funkce " )
        return NewSigma

    # First iteration of self-energy (complex)
    def FirstIterationofSelfEnergy(self):
        return 1j * self.Gamma_0() / (2 * self.t_value)
        # Zeroth iteration of Gamma


    def Gamma_0(self):
        parameter = (2 * np.pi * self.t_value) / self.U
        gamma_0 = -(2 * self.t_value) ** 2 * np.exp(parameter)
        # print(gamma_0)
        return gamma_0

    #Main Iteration Process
    def Iterationprosess(self):
        Sigma=self.FirstIterationofSelfEnergy()
        for i in range(self.NumIteration):
            #I divided cases of NumIteration==0 and other case
                a_0=self.a(Sigma)
                Gamma_0=self.Gamma(Sigma)
                NewSigma=self.NewSigma(self.omegavalues, Sigma)
                # Tolerance check:
                if abs(NewSigma - Sigma) < self.Tolerance:
                    break
                Sigma = NewSigma
        return Sigma






















