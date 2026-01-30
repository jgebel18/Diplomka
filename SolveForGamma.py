import numpy as np
from MakeFunctions import MakeFunctions
from scipy.interpolate import UnivariateSpline
from scipy.optimize import newton
from numpy.lib.scimath import sqrt as csqrt
from scipy.interpolate import PchipInterpolator
from PlotingFunctions import PlottigFunctions
from MakeFunction_Gamma import MakeFunctions_Gamma
from Class_for_testing_calculations import TestingCalculations
class MakeIterationProcessforGamma:

  #Here is a construction of class variables
    def __init__(self, beta, x_values, t_value, omega_values, U,
                 NumIterations, Tolerance, Rezolution,Gamma_0, Sigma_0):
        self.beta = beta            # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value #Hopping parameter
        self.omegavalues=omega_values #omega for Greeen function                # Hopping parameter
        self.U = U   #Electron Iteraction
        self.Mkf= MakeFunctions_Gamma(self.beta, self.x_values, self.t_value, self.omegavalues,self.U) # This calling of class is important for define another iterations
        self.NumIteration=NumIterations #Num Iterations
        self.Rezolution=Rezolution #Rezolution
        self.Tolerance=Tolerance #Tolerance
        self.PF=PlottigFunctions() # Calling of class with plotiing methods
        self.Gamma_0= 0
        self.step= 2
        self.num_steps=5
        self.Sigma_0=Sigma_0
        self.Gamma_values_left=np.linspace(self.Gamma_0-self.num_steps*self.step, self.Gamma_0, self.num_steps+1)
        self.Gamma_values_right= np.linspace(self.Gamma_0, self.Gamma_0+0.05, self.num_steps+1)
        self.Gamma_values= np.unique(np.concatenate((self.Gamma_values_left,
                                                     self.Gamma_values_right)))
        self.TestingClass = TestingCalculations(self.beta, self.x_values, self.t_value, self.omegavalues, self.U,
                                                self.NumIteration, self.Tolerance,
                                                self.Rezolution, self.Gamma_0,
                                                self.Sigma_0)
        self.Gamma=self.Gamma_0
        self.a_values=None
        self.D_values= None
        self.D2_values_Sigma= None


        # This is a First equation for a
    def a(self, Sigma):
        # print(self.Mkf.Y(Sigma))
        a = 1.0 + self.U * self.Mkf.Y(Sigma)
        print(a)
        return a

    # Another Iterations of Gamma
    def Gamma(self, Sigma):
        DotofItegrals = csqrt(2 / (self.Mkf.D(Sigma) * self.U))
        Gamma = (self.U / (self.beta * np.pi)) * DotofItegrals * (1 / csqrt(self.a(Sigma)))
        # self.TestingClass.CheckDominantTerminD()
        # self.Mkf.CollectDataAndPlot()
        return Gamma

    # little Gamma
    def gamma(self, Gamma):
        value = -Gamma / (2 * self.t_value) ** 2
        # print(value)
        return value

    # W-coeficent in the dependence of omega
    def W(self, omega):
        return omega / (2 * self.t_value)

    # Quartic equation
    def QuarticEquation(self, s, Gamma, omega):
        equation = s ** 4 - 2j * self.W(omega) * s ** 3 + (1 - self.W(omega) ** 2) * s ** 2 - self.gamma(Gamma) ** 2
        return equation

    # Quartic equation derivative
    def QuarticEquationDerivative(self, s, Gamma, omega):
        equation = 4 * s ** 3 - 6j * self.W(omega) * s ** 2 + 2 * (1 - self.W(omega) ** 2) * s
        return equation

    # Phi function for simple iteration method
    def QuarticFunctionPhi(self, s, Gamma, omega):
        numerator = (self.gamma(Gamma) ** 2)
        denom = (s ** 3 - 2j * self.W(omega) * s ** 2 + (1 - self.W(omega) ** 2) * s)
        return numerator / denom

    # Nonlinear equation solver
    def GiveSolutionofNonlinear(self, Gamma, omega):
        # solution = nonlin_equations.newton(f= self.QuarticEquation,
        #                                  df=self.QuarticEquationDerivative, args=(Gamma,omega) , x0=0.95+0j)
        solution = newton(func=self.QuarticEquation, args=(Gamma, omega),
                          fprime=self.QuarticEquationDerivative, x0=2.0, maxiter=1000)
        root = solution
        return root

    # Quartic equation solver
    def GiveSolutionofQuartic(self, Gamma, omega):
        coeffs = [1, -2j * self.W(omega), (1 - self.W(omega) ** 2), 0, -(self.gamma(Gamma)) ** 2]
        # coeffs = [ -2j , -2*(self.W(omega)-1 ), 0, -(self.gamma(Gamma)) ** 2]
        return np.roots(coeffs)

    # Complete values of root for every values of omega
    def GiveValuesOfRoots(self, Gamma):
        ValuesOfRoots = np.zeros(len(self.omegavalues), dtype=complex)
        ValuesOfRoots_2 = np.zeros((len(self.omegavalues), 4), dtype=complex)
        for i, omega in enumerate(self.omegavalues):
            ValuesOfRoots[i] = self.GiveSolutionofNonlinear(Gamma, omega)
            ValuesOfRoots_2[i] = self.GiveSolutionofQuartic(Gamma, omega)
        self.ZeroOmegaCase(Gamma)
        #print(self.GiveSolutionofNonlinear(Gamma, 0.0))
        return ValuesOfRoots, ValuesOfRoots_2

    # Trhis function return zero omega cese which is more comfortable to solve analithicaly
    def ZeroOmegaCase(self, Gamma):
        values_real_1, values_real_2 = (+csqrt(-0.5 + csqrt((Gamma / (2 * self.t_value) ** 2) ** 2 + 0.25)),
                                        -csqrt(-0.5 + csqrt((Gamma / (2 * self.t_value) ** 2) ** 2 + 0.25)))
        values_imag_1, values_imag_2 = (+1j * csqrt(0.5 + csqrt((Gamma / (2 * self.t_value) ** 2) ** 2 + 0.25)),
                                        -1j * csqrt(-0.5 + csqrt((Gamma / (2 * self.t_value) ** 2) ** 2 + 0.25)))
        #print(np.array([values_real_1, values_real_2, values_imag_1, values_imag_2]))

    # This function generate the solution of both equation
    def GenerateSolution(self, Gamma):
        values, values_quartic = self.GiveValuesOfRoots(Gamma)
        values = values.reshape(-1, 1)
        complete_values = np.hstack([-1j * values * 2 * self.t_value, -1j * values_quartic * 2 * self.t_value])
        #self.PF.PlotRootsOfEquation(complete_values.T, self.omegavalues)
        return -1j * values * 2 * self.t_value

    # Here is a function to plotting every parameter of this Iteration Process
    def PrintEveryParameter(self, i, Sigma):
        print(f"""Values of parameters in {i}-Iteration
            β = {self.beta}
            U = {self.U}
            a = {self.a(Sigma)}
            Γ = {self.Gamma(Sigma)}
            D = {self.Mkf.D(Sigma)}
            Y = {self.Mkf.Y(Sigma)}	
            Σ_{0} = {self.Sigma_0}
            Γ_{0} = {self.Gamma_0}""")

    # Main Iteration Process I have upgraded and I added an array to collect values in every different iteration
    def Iterationprocess(self):
        Sigma = self.Sigma_0
        SigmaIterations = []
        SigmaValues = np.full(len(self.omegavalues), Sigma, dtype=complex)
        NewSigma = None
        # SigmaIterations.append(SigmaValues)
        for i in range(self.NumIteration):
            Sigma = self.InterpolateOfSigmaValues(SigmaValues, self.omegavalues)
            self.Mkf.CollectDataAndPlot(Sigma)
            self.PF.PlotSigmaFunction(self.omegavalues, Sigma, SigmaValues)
            self.PrintEveryParameter(i, Sigma)
            SigmaIterations.append(Sigma)
            Gamma = self.Gamma(Sigma)
            NewSigma = self.GenerateSolution(Gamma)  # self.NewSigma(self.omegavalues, Sigma, i) #Calculating Sigma

            # Tolerance check
            if self.CheckingforContinue(SigmaValues, NewSigma) == True:
                return NewSigma, np.array(SigmaIterations)
            else:
                SigmaValues = NewSigma  # aktualization of sigma vlaues
        return NewSigma, np.array(SigmaIterations)

    # Here is the Interpolation of SelfEnergy. I decided for start with this iterpolation function. In future
    # I will probably change it

    # Sunil Here is the key method which is in the middle of our interest
    def InterpolateOfSigmaValues(self, SigmaValues, OmegaValues, ):
        # Interpolace reálné a imaginární části zvlášť
        SplineOmegaImag = PchipInterpolator(OmegaValues,
                                            SigmaValues.imag)  # UnivariateSpline(OmegaValues, np.imag(SigmaValues), s=0)
        SplineOmegaReal = PchipInterpolator(OmegaValues,
                                            SigmaValues.real)  # UnivariateSpline(OmegaValues, np.real(SigmaValues), s=0)
        SplineOmega = lambda x: SplineOmegaReal(x) + 1j * SplineOmegaImag(x)
        return SplineOmega

    # Here is a method checking for continue of comparing diff betwen actual and past iteration of Sigma
    # We have a sorted caseses for imag part and real part
    def CheckingforContinue(self, SigmaValues, NewSigmaValues):
        diff_real = np.std((NewSigmaValues.real - SigmaValues.real))
        diff_imag = np.std((NewSigmaValues.imag - SigmaValues.imag))
        print(f"std_real={diff_real}")
        print(f"std_imag={diff_imag}")
        diff = np.std(np.abs(NewSigmaValues - SigmaValues))
        print(f"std={diff}")
        if diff_real < self.Tolerance and diff_imag < self.Tolerance:
            return True
        else:
            return False

    # This function Gives  Green function
    def GiveFinalG(self):
        self.D2_values_Sigma=np.array([self.GenerateSolution(Gamma) for Gamma in self.Gamma_values])
        Sigma_approxmations= np.array([self.InterpolateOfSigmaValues(Sigma, self.omegavalues) for Sigma in (self.D2_values_Sigma)])


        self.a_values= np.array([self.Mkf.Y(Sigma) for Sigma in Sigma_approxmations])
        self.D_values= np.array([self.Mkf.D(Sigma ) for Sigma in Sigma_approxmations])
        D_approx_values = np.array([self.Mkf.D_approx(Gamma) for Gamma in self.Gamma_values])
        self.PF.Plot_Values_of_a_and_D(self.a_values, self.D_values,D_approx_values,  self.Gamma_values)


        #SigmaValues, SigmaIterations = self.Iterationprocess()
        # Sigma= self.InterpolateOfSigmaValues( SigmaValues,self.omegavalues)
        # self.TestingClass.CheckLeftSideIntegrand()
        # self.TestingClass.CheckDominantTerminD()
        # self.PF.PlotSigmaFunction(self.omegavalues, SigmaIterations)
        # self.PF.PlotGreenFunction(self.omegavalues, self.Mkf.G,  SigmaIterations,SigmaIterations.size)