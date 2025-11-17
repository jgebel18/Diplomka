import numpy as np
from fontTools.misc.bezierTools import epsilon
from scipy.interpolate import UnivariateSpline
from PlotingFunctions import PlottigFunctions
from MakeFunctions import MakeFunctions
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import quad, quad_vec

matplotlib.use('TkAgg')

from numpy.lib.scimath import sqrt as  csqrt

class TestingCalculations:
    def __init__(self, beta, x_values, t_value, omega_values, U, NumIterations, Tolerance, Rezolution, Gamma_0, Sigma_0):
        self.beta = beta  # Inverse temperature
        self.x_values = x_values  # Array of x values
        self.t_value = t_value  # Hopping parameter
        self.omegavalues = omega_values  # omega for Greeen function                # Hopping parameter
        self.U = U  # ElectronIteraction
        self.epsilon= None
        self.Mkf = MakeFunctions(self.beta, self.x_values, self.t_value, self.omegavalues,
                                 self.U)  # This calling of class is important for define another iterations
        self.Gamma_0=(Gamma_0)
        self.Sigma_0=(Sigma_0)
        self.NumIteration = NumIterations
        self.Rezolution = Rezolution
        self.Tolerance = Tolerance
        self.PF = PlottigFunctions()
        self.beta_test=np.linspace(0.1, 30, 100)*1/np.abs(self.Sigma_0)




    def SetValueEpsilon(self):
        Sigma = self.Sigma_0
        Helparray= np.linspace(0,10, 100)
        self.epsilon = Helparray*np.abs(Sigma)

    def TestingIntegrandLeftSide(self, x, Sigma, beta,epsilon):
        FirstPart = np.tanh((beta * x) / 2)
        SecondPart = (2 / (x - (1) * Sigma-epsilon) ** 3)
        return FirstPart * SecondPart.imag


    def Resultofanalythicalbetainfinity(self,):
        FirstPart= -4*(self.Gamma_0/(2*self.t_value))*self.epsilon
        SecondPart= ((self.Gamma_0/(2*self.t_value))**2+self.epsilon**2)
        return FirstPart/(SecondPart**2)

    def CheckLeftSideIntegrand(self ):
        self.SetValueEpsilon()
        Sigma= self.Sigma_0
        #print(self.epsilon)
        list_results_of_integral=np.zeros(len(self.epsilon))
        for i in range(len(self.epsilon)):
            integrand_left_side= lambda x:self.TestingIntegrandLeftSide(x, Sigma, self.beta, self.epsilon[i])
            result_left_side, error_left_side= quad_vec(integrand_left_side, -np.inf, np.inf)
            #print(f"Integral for beta={self.beta}:",result_left_side)
            list_results_of_integral[i]=result_left_side
        IdealCase=self.Resultofanalythicalbetainfinity()
        print("Difference:=", (list_results_of_integral-IdealCase)/IdealCase)
        self.PF.PlotResultsofIntegrals( self.epsilon,list_results_of_integral,IdealCase)
        #print(f"Case for beta goes to infinity",IdealCase )

    def Dominant_D_Result(self, Sign):
        alpha = self.Gamma_0 / (2 * self.t_value) ** 2
        z = Sign * 1j * ((-alpha) - csqrt(1 + alpha ** 2))
        FirstTerm = -(self.t_value ** 3 / (np.pi * self.Gamma_0 ** 2))
        FirstRatio = (((Sign * 5j * z) * (-alpha) - 3) / ((Sign * 1j * z) * (-alpha) - 1))
        SecondRatio = ((Sign * 1j * z) / (csqrt(1 + alpha ** 2)))
        Finalresult = FirstTerm * (Sign * 1j * z - alpha * (FirstRatio + SecondRatio + 1))
        return Finalresult


    def DominantDtermiinbetagoestoinfinity(self, Sign ):
        if Sign==(-1) or Sign==1:
            Dominant_D_result=self.Dominant_D_Result(Sign)
            return Dominant_D_result
        else:
            print("Only -1 or 1 values are free")

    def ThirdDRation(self, x, Sigma):
        FirstPartofDenom = (-x + Sigma)
        denom = (FirstPartofDenom * csqrt(4 * self.t_value ** 2 - (-(x) + Sigma) ** 2)) ** 3
        frac = (8 * self.t_value ** 2 - 3 * FirstPartofDenom ** 2) / denom
        return 8 * (self.t_value ** 2) * frac.real

    def CheckDominantTerminD(self):
        Sigma= self.Sigma_0
        #Helparray= np.linspace(-10,0, 100)*np.abs(Sigma)-2
        Integrand_Dominant_term_D=lambda x: (1/np.pi)*(1/16)*self.Mkf.FermiFunction(x)*self.ThirdDRation(x,Sigma)
        y = np.array([Integrand_Dominant_term_D(xi) for xi in Helparray])
        self.PF.PlotData(Helparray,y)
        Result_Dominant_D, Error_Dominant_D= quad(Integrand_Dominant_term_D, -np.inf, np.inf)
        #print(f"Main term integral for beta={self.beta}", Result_Dominant_term_D)
        IdealCaseDominantTerm= self.Dominant_D_Result(1)
        print(f"Difference", (Result_Dominant_D-IdealCaseDominantTerm)/IdealCaseDominantTerm)