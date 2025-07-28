import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
from numpy import number
from numpy.lib.scimath import sqrt as csqrt
from scipy.stats import trapezoid
from scipy.interpolate import interp1d
from scipy.integrate import tanhsinh
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



    def Gamma_0(self):
        parameter = (2 * np.pi * self.t_value) / self.U
        gamma_0 = -(2 * self.t_value) ** 2 * np.exp(parameter)
        # print(gamma_0)
        return gamma_0


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
        return frac.real

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

    def TestingIntegrandLeftSide(self, x, Sigma, beta):
        FirstPart = np.tanh((beta * x) / 2)
        SecondPart = (2 / (x - (1 + 1j) * Sigma) ** 3)
        # self.PF.PlotData('Left Side Integrand First Part', self.x_values,FirstPart, 'x',
        #            'Integrand', 'Integrand Left Side')
        # self.PF.PlotData('Left Side Integrand Second Part', self.x_values, SecondPart.imag, 'x',
        #             'Integrand', 'Integrand Left Side')
       # self.PF.PlotTwoSubplots('First Part', 'Second Part', self.x_values, FirstPart, SecondPart.imag, 'x', 'Integrands',
                             #'Both parts left side')
        return FirstPart * SecondPart.imag

    def TestingIntegrandRightSide(self, x, Sigma ,beta):
        FirstPart = (((beta) / 2) ** 2 * (-2) * np.tanh((beta * x) / 2) * (
                    1 - np.tanh((beta * x) / 2) ** 2))
        SecondPart = (1 / (x - (1 + 1j) * Sigma))
        #self.PlotData('Right Side Integrand First Part', self.x_values, FirstPart, 'x',
                     #'Integrand', 'Integrand Right Side')
        #self.PlotData('Right Side Integrand Second Part', self.x_values, SecondPart.imag , 'x',
                    #'Integrand', 'Integrand Right Side')
        #self.PF.PlotTwoSubplots('First Part', 'Second Part', self.x_values, FirstPart, SecondPart.imag,
         #                    'x', 'Integrands', 'Both parts left side')
        return FirstPart * SecondPart.imag

    # def InterpolateValuesofSigma(self, Sigma_values ):
    # Calculate Integral in D
    def D(self, Sigma):
        integrand = lambda x: (1 / np.pi) * self.FermiFunction(x) * self.SumOfComplexRations(x, Sigma)
        #values of results of integrals
        data_results_lhs=np.zeros((len(self.beta)))
        data_results_rhs=np.zeros((len(self.beta)))
        #values of Integrands for next plotting
        integrand_lhs_values = np.zeros((self.beta.size, self.x_values.size))
        integrand_rhs_values = np.zeros((self.beta.size, self.x_values.size))

        #Generation of arrays of  results and arrays of values of integrands for
        # different values of beta
        for i in range(self.beta.size):
            integrand_lhs= lambda x: self.TestingIntegrandLeftSide(x, Sigma, self.beta[i])
            integrand_rhs= lambda x: self.TestingIntegrandRightSide(x, Sigma, self.beta[i])
            result_lhs, error_lhs = integrate.quad_vec(integrand_lhs, -np.inf, np.inf)
            result_rhs, error_rhs = integrate.quad_vec(integrand_rhs, -np.inf, np.inf)
            integrand_lhs_values[i,:]=self.TestingIntegrandLeftSide(self.x_values, Sigma, self.beta[i])
            integrand_rhs_values[i,:]=self.TestingIntegrandRightSide(self.x_values, Sigma, self.beta[i])
            data_results_lhs[i]= result_lhs
            data_results_rhs[i]= result_rhs
        #Plotting of Integrals
        #self.PF.PlotIntegrandsfordifferentbetavalues(self.beta,self.x_values, integrand_lhs_values, integrand_rhs_values)
        result, error = integrate.quad_vec(integrand, -np.inf, np.inf)
        print((self.Gamma_0()/(2*self.t_value))**(-2))
        #Writing into txt file
        self.WritetoTxtFile(Sigma)
        #Ploting results of integrals
        self.PF.PlotREsultsofIntegrals(self.beta, data_results_lhs, data_results_rhs)
        print((-2 * self.t_value ** 3) / (3 * np.pi * self.Gamma_0() ** 2))
        return result

    def WritetoTxtFile(self, Sigma):
        lhs_all = []  # left side
        rhs_all = []  # right side
        for beta_i in self.beta:
            lhs = self.TestingIntegrandLeftSide(self.x_values, Sigma, beta_i)
            rhs = self.TestingIntegrandRightSide(self.x_values, Sigma, beta_i)
            lhs_all.append(lhs)
            rhs_all.append(rhs)

        # Převedeme na NumPy pole (každý řádek odpovídá jedné hodnotě beta)
        lhs_array = np.array(lhs_all).T  # shape (n_x, n_beta)
        rhs_array = np.array(rhs_all).T
        result = np.column_stack([self.x_values, lhs_array, rhs_array])

        # Sloupce pojmenujeme: x, LHS_beta0, LHS_beta1, ..., RHS_beta0, ...
        header_labels = ['x']
        header_labels += [f'LHS_beta={b:.3f}' for b in self.beta]
        header_labels += [f'RHS_beta={b:.3f}' for b in self.beta]
        header = ' '.join(header_labels)
        # Saving
        np.savetxt('AllIntegrands.txt', result, header=header, comments='', delimiter=' ')

    #Calculate Integral in Y(0,0)
    def Y(self,Sigma):
        integrand = lambda x: (1 / np.pi) * self.FermiFunction(x) * self.ComplexFrac_in_Y(x,Sigma)
        result,error = integrate.quad(integrand, -np.inf,  np.inf)
        return result


    #Here id a Green Function in the dependence of complex omega
    def G(self, omega ,Sigma):
        #Sigma=self.FirstIterationofSelfEnergy()
        reader= -1j
        denominator=csqrt(4 * self.t_value**2 - (-(omega) + Sigma)**2)
        frac=reader/denominator
        return frac#, frac.imag

