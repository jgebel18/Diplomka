import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
from numpy.lib.scimath import sqrt as csqrt


class MakeFunctions:
    # Constructor with physical parameters
    def __init__(self, beta, x_values, t_value, omega_values, U):
        self.beta = beta            # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value     # Hopping parameter
        self.omega_vlaues= omega_values
        self.U = U      # Electron interaction

    # Zeroth iteration of Gamma
    def Gamma_0(self):
        parameter = (2 * np.pi * self.t_value) / self.U
        gamma_0 = -(2 * self.t_value)**2 * np.exp(parameter)
        #print(gamma_0)
        return gamma_0


    #Possible 1st iteration of Self-energy
    def PossibeSigma(self,x):
        return 1j*self.Gamma_0()/(np.sqrt(4*self.t_value**2- x**2))


    # First iteration of self-energy (complex)
    def FirstIterationofSelfEnergy(self):
        return 1j * self.Gamma_0() / (2 * self.t_value)

    # Imaginary part in integrand of Y
    def ComplexFrac_in_Y(self, x):
        Sigma = self.FirstIterationofSelfEnergy()
        denom = csqrt(4 * self.t_value**2 - (-(x) + Sigma)**2)
        #print((1/csqrt(4 * self.t_value**2 - (-(-2) + Sigma)**2)).real)
        frac = 1 / (denom * (-x + Sigma))
        return frac.real

    # Fermi function
    def FermiFunction(self, x):
        return 1 / (np.exp(self.beta * x) + 1)


    # When we into D eqution in the outline for numerics we see that D is defined by Fermi function
    # and three rations here we defined all of them part by part and Here is a First part
    def FirstDRatio(self, x):
        Sigma=self.FirstIterationofSelfEnergy()
        First_ratio= 2*self.ComplexFrac_in_Y(x)
        return First_ratio


    # Here is a second Part of
    def SecondDRation(self,x):
        Sigma=self.FirstIterationofSelfEnergy()
        denom=(csqrt(4 * self.t_value**2 - (-(x) + Sigma)**2))**3
        frac=(x-Sigma)/denom
        return frac.real

    #Here is a Third Part
    def ThirdDRation(self,x):
        Sigma=self.FirstIterationofSelfEnergy()
        FirstPartofDenom=(-x+Sigma)
        denom = (FirstPartofDenom*csqrt(4 * self.t_value ** 2 - (-(x) + Sigma) ** 2)) ** 3
        frac= (8*self.t_value**2-3*FirstPartofDenom**2)/denom
        return 2*frac.imag

    #Here is a sumation of alll rations
    def SumOfComplexRations(self, x):
        sum= ((1/16)*(self.FirstDRatio(x)+self.SecondDRation(x)+self.ThirdDRation(x)))
        self.PlotIntegrands('Integrands of D', self.x_values, self.FirstDRatio(self.x_values),'1st-Ratio' ,
                            self.SecondDRation(self.x_values),'2nd-Ratio', self.ThirdDRation(self.x_values), '3rd-Ratio')
        return sum

    #Calculate Integral in D
    def D(self):
        integrand=lambda x:(1/np.pi)*self.FermiFunction(x)*self.SumOfComplexRations(x)
        result=integrate.quad(integrand,-np.inf, np.inf)
        return result

    #Calculate Integral in Y(0,0)
    def Y(self):
        integrand = lambda x: (1 / np.pi) * self.FermiFunction(x) * self.ComplexFrac_in_Y(x)
        result = integrate.quad(integrand, -np.inf,  np.inf)
        return result


    #Here id a Green Function in the dependence of complex omega
    def G(self, omega):
        Sigma=self.FirstIterationofSelfEnergy()
        reader= -1j
        denominator=csqrt(4 * self.t_value**2 - (-(omega) + Sigma)**2)
        frac=reader/denominator
        return frac.real#, frac.imag





    #I added here function to plot all integrnds in one picture and every part which is important to define
    # D function
    def PlotIntegrands(self, nazev, DataX, DataY, nazev1, DataY2, nazev2, DataY3, nazev3):
        plt.figure(figsize=(10, 10))
        plt.plot(DataX, DataY, label=nazev1)
        plt.plot(DataX, DataY2, label=nazev2)
        plt.plot(DataX, DataY3, label=nazev3)
        plt.title(nazev)
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('Integrands')
        plt.tight_layout()
        plt.savefig(nazev + '.png')  # uložit před zobrazením
        plt.show()

    # Here is another plotting function which is ploting one plot in one picture
    def PlotData(self, Nazev, dataX, dataY, Nazev_X, Nazev_Y, Title):
        plt.figure(figsize=(10, 10))
        plt.plot(dataX, dataY, label=Nazev)
        plt.xlabel(Nazev_X)
        plt.ylabel(Nazev_Y)
        plt.title(Title)
        plt.savefig(Title + 'pgn')
        plt.show()