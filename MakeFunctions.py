import numpy as np
import scipy.integrate as integrate
from numpy.lib.scimath import sqrt as csqrt
from PlotingFunctions import PlottigFunctions


class MakeFunctions:
    # Constructor with physical parameters
    def __init__(self, beta, x_values, t_value, omega_values, U):
        self.beta = beta            # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value     # Hopping parameter
        self.omega_vlaues= omega_values #omega values
        self.U = U      # Electron interaction
        self.PF= PlottigFunctions()




    # Imaginary part in integrand of Y
    def ComplexFrac_in_Y(self, x, Sigma):
        denom = csqrt(4 * self.t_value**2 - (-(x) + Sigma)**2)
        #print((1/csqrt(4 * self.t_value**2 - (-(-2) + Sigma)**2)).real)
        frac = 1 / (denom * (-x + Sigma))
        return frac.real

    # Fermi function
    def FermiFunction(self, x):
        return 1 / (np.exp(self.beta * x) + 1)

    # When we into D eqution in the outline for numerics we see that D is defined by Fermi function
    # and three rations here we defined all of them part by part and Here is a First part
    def FirstDRatio(self, x,Sigma):
        #Sigma=self.FirstIterationofSelfEnergy()
        First_ratio= 2*self.ComplexFrac_in_Y(x, Sigma)
        return First_ratio

    # Here is a second Part of
    def SecondDRation(self,x,Sigma):
        #Sigma=self.FirstIterationofSelfEnergy()
        denom=(csqrt(4 * self.t_value**2 - (-(x) + Sigma)**2))**3
        frac=(x-Sigma)/denom
        return -2*frac.real

    #Here is a Third Part
    def ThirdDRation(self,x,Sigma):
        FirstPartofDenom=(-x+Sigma)
        denom = (FirstPartofDenom*csqrt(4 * self.t_value ** 2 - (-(x) + Sigma) ** 2)) ** 3
        frac= (8*self.t_value**2-3*FirstPartofDenom**2)/denom
        return 8*(self.t_value**2)*frac.real

    #Here is a sumation of alll rations
    def SumOfComplexRations(self, x, Sigma):
        sum= ((1/16)*(self.FirstDRatio(x,Sigma)+self.SecondDRation(x,Sigma)+self.ThirdDRation(x,Sigma)))
        #self.PF.PlotIntegrands('Integrands of D', self.x_values, self.FirstDRatio(self.x_values,Sigma),'1st-Ratio' ,
                            #self.SecondDRation(self.x_values,Sigma),'2nd-Ratio', self.ThirdDRation(self.x_values,Sigma), '3rd-Ratio')
        return sum

    #Here is a function of Test integrand of left side
    def TestingIntegrandLeftSide(self, x, Sigma, beta):
        FirstPart = np.tanh((beta * x) / 2)
        SecondPart = (2 / (x - (1 + 1j) * Sigma) ** 3)
        return FirstPart * SecondPart.imag


    #Here is a testing integrand of left side
    def TestingIntegrandRightSide(self, x, Sigma ,beta):
        FirstPart = (((beta) / 2) ** 2 * (-2) * np.tanh((beta * x) / 2) * (
                    1 - np.tanh((beta * x) / 2) ** 2))
        SecondPart = (1 / (x - (1 + 1j) * Sigma))
        return FirstPart * SecondPart.imag

    # def InterpolateValuesofSigma(self, Sigma_values ):
    # Calculate Integral in D
    def D(self, Sigma):
        integrand = lambda x: (1 / np.pi) * self.FermiFunction(x) * self.SumOfComplexRations(x, Sigma(x))
        result, error = integrate.quad_vec(integrand, -np.inf,np.inf)
        print(result)
        #print(-4*(self.Gamma_0()/2))
        return result

    #Calculate Integral in Y(0,0)
    def Y(self,Sigma):
        integrand = lambda x: (1 / np.pi) * self.FermiFunction(x) * self.ComplexFrac_in_Y(x,Sigma(x))
        result,error = integrate.quad(integrand, -np.inf,  np.inf)
        return result


    #Here id a Green Function in the dependence of complex omega
    def G(self, omega ,Sigma):
        #Sigma=self.FirstIterationofSelfEnergy()
        reader= -1j
        denominator=csqrt(4 * self.t_value**2 - ((-1)*(omega) + Sigma(omega))**2)
        frac=reader/denominator
        return frac#, frac.imag