import numpy as np
from numpy.lib.scimath import sqrt as csqrt


class Main_Terms_Equations:
    def __init__(self, beta, x_values, t_value, omega_values, U,
                 NumIterations, Tolerance, Rezolution):
        self.beta = beta  # Inverse temperature
        self.x_values = x_values  # Array of x values
        self.t_value = t_value  # Hopping parameter
        self.omegavalues = omega_values  # omega for Greeen function                # Hopping parameter
        self.U = U  # Electron Iteraction
        # This calling of class is important for define another iterations
        self.NumIteration = NumIterations  # Num Iterations
        self.Rezolution = Rezolution  # Rezolution
        self.Tolerance = Tolerance  # Tolerance



    # This is a First equation for a
    def a(self, Gamma, Mkf):
        # print(self.Mkf.Y(Sigma))
        a = 1.0 + self.U * Mkf.Y(Gamma)
        # print(a)
        return a

    # especially modified function for D in
    # dependence of one a and D
    def Gamma_function_a_D(self, beta, D, a):
        DotofItegrals = csqrt(2 / (D * self.U))
        Gamma_function = (self.U / (beta * np.pi)) * DotofItegrals * (1 / csqrt(a))
        # self.TestingClass.CheckDominantTerminD()
        # self.Mkf.CollectDataAndPlot()
        return Gamma_function

    def GammaNonlinearFunction(self, Gamma, beta, Mkf):
        Gamma_function = self.Gamma_function(Gamma, beta, Mkf)
        return Gamma - Gamma_function

    # Another Iterations of Gamma
    def Gamma_function(self, Gamma, beta, Mkf):
        DotofItegrals = csqrt(2 / (Mkf.D(Gamma) * self.U))
        Gamma_function = (self.U / (beta * np.pi)) * DotofItegrals * (1 / csqrt(self.a(Gamma, Mkf)))
        # self.TestingClass.CheckDominantTerminD()
        # self.Mkf.CollectDataAndPlot()
        return Gamma_function

    # Trhis function return zero omega cese which is more comfortable to solve analithicaly
    def ZeroOmegaCase(self, Gamma):
        values_real_1, values_real_2 = (+csqrt(-0.5 + csqrt((Gamma / (2 * self.t_value) ** 2) ** 2 + 0.25)),
                                        -csqrt(-0.5 + csqrt((Gamma / (2 * self.t_value) ** 2) ** 2 + 0.25)))
        values_imag_1, values_imag_2 = (+1j * csqrt(0.5 + csqrt((Gamma / (2 * self.t_value) ** 2) ** 2 + 0.25)),
                                        -1j * csqrt(0.5 + csqrt((Gamma / (2 * self.t_value) ** 2) ** 2 + 0.25)))
        print(np.array([values_real_1, values_real_2, values_imag_1, values_imag_2]))
        return np.array([values_real_1, values_real_2, values_imag_1, values_imag_2])

    # Here is a function to plotting every parameter of this Iteration Process
    def PrintEveryParameter(self, i, Sigma):
        print(f"""Values of parameters in {i}-Iteration
              β = {self.beta}
              U = {self.U}
              a = {self.a(Sigma)}
              Γ = {self.Gamma(Sigma )}
              D = {self.Mkf.D(Sigma)}
              Y = {self.Mkf.Y(Sigma)}	
              Σ_{0} = {self.Sigma_0}
              Γ_{0} = {self.Gamma_0}""")

    def NonlinearEquationa(self, Gamma, beta , a, D):

        DotsOfIntegrals= 2/(self.U*D)*(1/Gamma**2)
        SecondTerm= (self.U/(np.pi*beta))**2
        FinalEquation= SecondTerm*DotsOfIntegrals-a
        return FinalEquation

    def NonlinearEquationA(self, Gamma, beta , Mkf):
        a= self.a(Gamma, Mkf)
        DotsOfIntegrals= 2/(self.U*Mkf.D(Gamma))*(1/Gamma**2)
        SecondTerm= (self.U/(np.pi*beta))**2
        FinalEquation= SecondTerm*DotsOfIntegrals-a
        return FinalEquation


