import numpy as np
from MakeFunctions import MakeFunctions
from scipy.interpolate import UnivariateSpline
from scipy.optimize import newton
from scipy.optimize import root
from scipy.optimize import root_scalar
from numpy.lib.scimath import sqrt as csqrt
from scipy.interpolate import PchipInterpolator
from PlotingFunctions import PlottigFunctions
from MakeFunction_Gamma import MakeFunctions_Gamma
from Class_for_testing_calculations import TestingCalculations
class MakeIterationProcessforGamma:
  #Here is a construction of class variables
    def __init__(self, beta, x_values, t_value, omega_values, U,
                 NumIterations, Tolerance, Rezolution,Gamma_0, Sigma_0):
        self.beta = beta           # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value #Hopping parameter
        self.omegavalues=omega_values #omega for Greeen function                # Hopping parameter
        self.U = U   #Electron Iteraction
         # This calling of class is important for define another iterations
        self.NumIteration=NumIterations #Num Iterations
        self.Rezolution=Rezolution #Rezolution
        self.Tolerance=Tolerance #Tolerance
        self.PF=PlottigFunctions() # Calling of class with plotiing methods
        self.Gamma_0=0
        self.step= 0.01
        self.num_steps=100
        self.Path_Files='Files'
        self.Path_Images='Images'
        self.Sigma_0=Sigma_0
        self.Gamma_values=np.linspace(self.Gamma_0-self.num_steps*self.step, self.Gamma_0, self.num_steps+1)
        #self.Gamma_values_right= np.linspace(self.Gamma_0, self.Gamma_0+self.num_steps*self.step, self.num_steps+1)
        #np.unique(np.concatenate((self.Gamma_values_left,                    #                         self.Gamma_values_right)))
        self.TestingClass = TestingCalculations(self.beta, self.x_values, self.t_value, self.omegavalues, self.U,
                                                self.NumIteration, self.Tolerance,
                                                self.Rezolution, self.Gamma_0,
                                                self.Sigma_0)


        self.D2_values_Sigma= None


        # This is a First equation for a
    def a(self, Gamma, Mkf):
        # print(self.Mkf.Y(Sigma))
        a = 1.0 + self.U * Mkf.Y(Gamma)
        #print(a)
        return a

    # especially modified function for D in
  # dependence of one a and D
    def Gamma_function_a_D(self, beta, D, a):

        DotofItegrals = csqrt(2 / (D * self.U))
        Gamma_function = (self.U / (beta * np.pi)) * DotofItegrals * (1 / csqrt(a))
        # self.TestingClass.CheckDominantTerminD()
        # self.Mkf.CollectDataAndPlot()
        return Gamma_function

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




# I modified this function to gives 3D numpy.array for every important termi in solved equations
  #
    def GiveFinalG(self):
            # Inicializace jednoho 3D pole: (4 metriky, počet_beta, počet_gamma)
        res = np.empty((4, self.beta.size, self.Gamma_values.size))
        for i, beta in enumerate(self.beta):
            mkf = MakeFunctions_Gamma(beta, self.x_values, self.t_value, self.omegavalues, self.U)
            for j, gamma in enumerate(self.Gamma_values):
                res[0, i, j] = mkf.Y(gamma)
                res[1, i, j] = mkf.D(gamma)*gamma**2
                # This is for saving the calculation times I decide to calculate Y
                # and D once I also created special functions in this class
                res[2, i, j] = 1+self.U*res[0, i, j]#self.Gamma_function(gamma, beta, mkf)
                res[3, i, j] = self.Gamma_function_a_D(beta, res[1, i, j]/gamma**2 ,res[2, i, j] )
        self.PF.Plot_Beta_Gamma_Dependence(res, self.beta, self.Gamma_values)
        #self.PF.Plot_Values_of_a_and_D(Y_values, D_values,Y_approx,D_approx,
        #                               self.Gamma_values, self.beta)


        #SigmaValues, SigmaIterations = self.Iterationprocess()
        # Sigma= self.InterpolateOfSigmaValues( SigmaValues,self.omegavalues)
        # self.TestingClass.CheckLeftSideIntegrand()
        # self.TestingClass.CheckDominantTerminD()
        # self.PF.PlotSigmaFunction(self.omegavalues, SigmaIterations)
        # self.PF.PlotGreenFunction(self.omegavalues, self.Mkf.G,  SigmaIterations,SigmaIterations.size)