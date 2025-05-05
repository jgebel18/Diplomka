import scipy.integrate as integrate
import numpy as np
from numpy.lib.scimath import sqrt as csqrt


class MakeFunctions:
    # Constructor with physical parameters
    def __init__(self, beta, x_values, t_value, U):
        self.beta = beta            # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value      # Hopping parameter
        self.U = U                  # Electron interaction

    # Zeroth iteration of Gamma
    def Gamma_0(self):
        parameter = (2 * np.pi * self.t_value) / self.U
        gamma_0 = -(2 * self.t_value)**2 * np.exp(parameter)
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
        # Debug výpis odmocniny a jejího argumentu
        #print(f"x = {x}")
        print(f"denom = {denom}")
        print(f"arg(denom) = {np.angle(denom):.4f} rad")  # fáze ve výchozím intervalu (-pi, pi]
        frac = 1 / (denom * (-x + Sigma))
        return frac.real  # or .imag if we  want the imaginary part

    # Fermi function
    def FermiFunction(self, x):
        return 1 / (np.exp(self.beta * x) + 1)

    # Calculates Y(0,0) for 1st iteration of self-energy
    def Y(self):
        integrand = lambda x: (1 / np.pi) * self.FermiFunction(x) * self.ComplexFrac_in_Y(x)
        result = integrate.quad(integrand, -np.inf,  np.inf)
        return result

