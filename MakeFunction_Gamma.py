import numpy as np
import scipy.integrate as integrate
from numpy.lib.scimath import sqrt as csqrt
from PlotingFunctions import PlottigFunctions
from scipy.interpolate import PchipInterpolator
from Nonlinear_Equation_Gamma_Solver import Nonlinear_Equation_Solver
class MakeFunctions_Gamma:

    def __init__(self, beta, x_values, t_value, omega_values, U):
        self.beta = beta
        self.x_values = x_values
        self.t_value = t_value
        self.omega_values = omega_values
        self.U = U
        self.PF = PlottigFunctions()
        self.D_Iterations = []
        self.Y_Iterations = []
        self.Nonlinear_Solver = Nonlinear_Equation_Solver(self.t_value, self.omega_values)


    # ================= BASIC FUNCTIONS =================
    def FermiFunction(self, x):
        return 1 / (np.exp(self.beta * x) + 1)

#This function call functions from scipy.Interpolate and interpolate
    # both parts of calculated values of self-energy
    def Sigma_Interpolation(self, Sigma_values):
        Sigma_komplex_inter= PchipInterpolator(self.omega_values,
                                               Sigma_values.imag,)
        Sigma_real_inter= PchipInterpolator(self.omega_values,
                                            Sigma_values.real)
        Sigma= lambda x: Sigma_real_inter(x)+1j*Sigma_komplex_inter(x)
        return Sigma

#This function call Nonlinear-equation-Solver and find values of Sigma
    def Sigma(self, Gamma):
        values=self.Nonlinear_Solver.GiveSolutionofNonlinear(Gamma)
        #print('omega, roots',(self.omegavalues,values))
        #self.PF.PlotRootOfEquation(-1j * values * 2 * self.t_value, self.omega_values)
        Sigma_values=-1j * values * 2 * self.t_value
        Sigma= self.Sigma_Interpolation(Sigma_values)
        return Sigma



    def ComplexFrac_in_Y(self, x, Sigma):
        if np.ndim(Sigma) != 0:
            raise ValueError("Sigma must be scalar")

        denom = csqrt(4 * self.t_value**2 - (-(x) + Sigma)**2)
        return (1 / (denom * (-x + Sigma))).real

    # ================= D RATIOS =================
    def FirstDRatio(self, x, Sigma):
        return 2 * self.ComplexFrac_in_Y(x, Sigma)

    def SecondDRation(self, x, Sigma):
        denom = csqrt(4 * self.t_value**2 - (-(x) + Sigma)**2)**3
        return -2 * ((x - Sigma) / denom).real

    def ThirdDRation(self, x,Sigma):

        fp = (-x + Sigma)
        denom = (fp * csqrt(4 * self.t_value**2 - (-(x) + Sigma)**2))**3
        frac = (8 * self.t_value**2 - 3 * fp**2) / denom
        return 8 * self.t_value**2 * frac.real

    def SumOfComplexRations(self, x, Sigma):
        return (1/16) * (
            self.FirstDRatio(x, Sigma)
            + self.SecondDRation(x, Sigma)
            + self.ThirdDRation(x, Sigma)
        )

    # ================= INTEGRALS (SCALAR ONLY) =================

# for saving the calculation time I modified Integration function into one
    def Integrator(self, Gamma , function):
        Sigma= self.Sigma(Gamma)
        def integrand(x):
            sig=Sigma(float(x)).item()
            integrand= (1/np.pi)*self.FermiFunction(x)*function(x,sig)
            return integrand
        result,_= integrate.quad_vec(integrand, -np.inf, np.inf, )
        return result

# These funct are working only as the caller the
    # current  function integral
    def D(self, Gamma):
        result = self.Integrator(Gamma,self.SumOfComplexRations)
        return result

    def Y(self, Gamma):
        result=self.Integrator(Gamma, self.ComplexFrac_in_Y)
        return result

    def Y_approx(self, Gamma):
         denom=-Gamma/(2*self.t_value)**2
         first_term= 1/(2*np.pi*self.t_value)
         numerator= np.sqrt(1+(denom**2)) +1
         return first_term *np.log(numerator/denom)

    def D_approx(self, Gamma):
        const_in_front = 1 / (16 * np.pi * self.t_value)
        denom = -Gamma / (2 * self.t_value) ** 2
        first_term = np.log((np.sqrt(1 + (Gamma / (2 * self.t_value) ** 2) ** 2) + 1) / (denom))
        second_term = np.sqrt(1 + (Gamma / (2 * self.t_value) ** 2) ** 2) / denom ** 2
        result = const_in_front * (first_term - second_term)
        return result * Gamma ** 2

    # ================= ITERATION DATA =================

    def D_Integrand(self, Gamma):
        Sigma=self.Sigma(Gamma)
        values = np.zeros_like(self.x_values, dtype=complex)
        for i, x in enumerate(self.x_values):
            sig = Sigma(float(x)).item()
            values[i] = (1 / np.pi) * self.FermiFunction(x) * self.SumOfComplexRations(x, sig)

        self.D_Iterations.append(values)

    def Y_Integrand(self, Gamma):
        Sigma = self.Sigma(Gamma)
        values = np.zeros_like(self.x_values, dtype=complex)
        for i, x in enumerate(self.x_values):
            sig = Sigma((x)).item()
            values[i] = (1 / np.pi) * self.FermiFunction(x) * self.ComplexFrac_in_Y(x, sig)

        self.Y_Iterations.append(values)

    # ================= GREEN FUNCTION =================

    def G(self, omega, Sigma):
        sig = Sigma(float(omega))
        denom = csqrt(4 * self.t_value**2 - ((-omega + sig)**2))
        return (-1j / denom)

    # ================= PLOTTING =================

    def CollectDataAndPlot(self, Sigma):
        self.Y_Integrand(Sigma)
        self.D_Integrand(Sigma)
        D_array = np.vstack(self.D_Iterations)
        Y_array = np.vstack(self.Y_Iterations)
        self.PF.PlotingItegrands(D_array, Y_array, self.x_values)
