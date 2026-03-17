import numpy as np
import os
from scipy.optimize import root_scalar


#Hi Sunil I decided to create this class, because i neded to get our
# caluculations under control and for better orientation and clasification
# of other mistakes
class Nonlinear_Equation_Solver:
    def __init__(self, t, omegavalues,):
        self.t_value= t
        self.omega_values= omegavalues


#These functions are definig the key terms in nonlinear equation
    def gamma(self, Gamma):
        value = -Gamma / (2 * self.t_value) ** 2
        # print(value)
        return value

    def W(self, omega):
        return omega / (2 * self.t_value)
#This function gives the current value of nonlinear equation
    # for s eps , omega and gamma
    def NonlinearEquation(self, s, Gamma, eps, omega):
        W = self.W(omega) + eps * 1j
        numerator = self.gamma(Gamma)
        denominator = np.sqrt(1 + (s - (W) * 1j) ** 2)
        return s - numerator / denominator

#This function gives the fprime in Newtons method
    def NonlinearEquationDerivative(self, s, Gamma, eps, omega):
        W = (self.W(omega) + eps * 1j)
        numerator = -self.gamma(Gamma) * (s - 1j * W)
        denominator = (1 + (s - (W) * 1j) ** 2) ** (3 / 2)
        # denominator= (csqrt(1+(s-self.W(omega)*1j)**2)*csqrt(1+(s-self.W(omega)*1j)**2)
        #             *csqrt(1+(s-self.W(omega)*1j)**2))
        return (1 - numerator / denominator)


#This metod gives the solution of nonlinear equation
    #for every omega in omegavalues and unique Gamma
    def GiveSolutionofNonlinear(self, Gamma, ):
        complete_solution = np.empty(len(self.omega_values), dtype=complex)
        omega_0 = 1.0 + 1.0 * 1j
        eps=1e-8
        for i, omega in enumerate(self.omega_values):
            # solution = nonlin_equations.newton(f= self.QuarticEquation,
            #
            #                                  df=self.QuarticEquationDerivative, args=(Gamma,omega) , x0=0.95+0j)
            #This is the function from scipy and that function gives the solution of
            # nonlinear equation by Newton's method
            solution = root_scalar(f=self.NonlinearEquation, method='newton', args=(Gamma,eps, omega),
                                   fprime=self.NonlinearEquationDerivative, x0=omega_0, maxiter=1000)
            complete_solution[i] = solution.root
        return complete_solution

# Tis method tested the dependence of solution on different values of epsilon
    def GiveSolutionofNonlinearEps(self, Gamma, ):
        complete_solution = np.empty(( self.epsilon.size, self.omegavalues.size), dtype=complex)
        omega_0 = self.omegavalues[0]
        for i, epsilon  in enumerate(self.epsilon):
            for j, omega in enumerate(self.omegavalues):


                # solution = nonlin_equations.newton(f= self.QuarticEquation,
                #                                  df=self.QuarticEquationDerivative, args=(Gamma,omega) , x0=0.95+0j)
                solution = root_scalar(f=self.NonlinearEquation, method='newton', args=(Gamma,epsilon, omega),
                                        fprime=self.NonlinearEquationDerivative, x0=omega_0, maxiter=2000, xtol=1e-8)


                complete_solution[i][j] = solution.root
        self.PF.PlotRootOfEquationEps(-1j * complete_solution * 2 * self.t_value,
                                      self.omegavalues, self.epsilon)
        return complete_solution