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
    #Here is a Constructor with new special variables epsilon and betatest the values of these variables will be set in setters in this class
    def __init__(self, beta, x_values, t_value, omega_values, U, NumIterations, Tolerance, Rezolution):
        self.beta = beta  # Inverse temperature
        self.x_values = x_values  # Array of x values
        self.t_value = t_value  # Hopping parameter
        self.omegavalues = omega_values  # omega for Greeen function                # Hopping parameter
        self.U = U  # ElectronIteraction
        self.epsilon= None
        self.betatest=None
        self.Mkf = MakeFunctions(self.beta, self.x_values, self.t_value, self.omegavalues,
                                 self.U)  # This calling of class is important for define another iterations
        self.NumIteration = NumIterations
        self.Rezolution = Rezolution
        self.Tolerance = Tolerance
        self.PF = PlottigFunctions()

#Here is  a function for Zeroth iteration of Gamma
    def Gamma_0(self):
        parameter = (2 * np.pi * self.t_value) / self.U
        gamma_0 = -(2 * self.t_value) ** 2 * np.exp(parameter)
        #print(gamma_0)
        return gamma_0

# Here is a Zeroth iteration of Sigma
    def ZerothIterationofSelfEnergy(self):
        Sigma= 1j * self.Gamma_0() / (2 * self.t_value)
        print("Sigma = ", Sigma)
        return Sigma
#Here is a setter for test beta values
    def SetBetaTest(self):
        self.betatest=np.linspace(0.1, 10, 100)*1/np.abs(self.ZerothIterationofSelfEnergy())

#Here is a setterfor Epsilon
    def SetValueEpsilon(self):
        Sigma = self.ZerothIterationofSelfEnergy()
        self.epsilon = 1.2*np.abs(Sigma)


#Here is a function for testing integrands for D function for confim wich shape of equation is
    # most comfortable for numerical integration methods in Python
    def TestingIntegrandLeftSide(self, x, Sigma, beta,epsilon):
        FirstPart = np.tanh((beta * x) / 2)
        SecondPart = (2 / (x - (1) * Sigma-epsilon) ** 3)
        return FirstPart * SecondPart.imag

#Here is a function which inicializes idealized case for beta goes to infinity
    def Resultofanalythicalbetainfinity(self,):
        FirstPart= -4*(self.Gamma_0()/(2*self.t_value))*self.epsilon
        SecondPart= ((self.Gamma_0()/(2*self.t_value))**2+self.epsilon**2)
        return FirstPart/(SecondPart**2)

# Here is a function which comparsion idealized case and result of numerical integration
    def CheckLeftSideIntegrand(self ):
        #Calling setters of Epsilon/Beta
        self.SetValueEpsilon()
        self.SetBetaTest()
        Sigma= self.ZerothIterationofSelfEnergy()
        #print(self.epsilon)
        list_results_of_integral=np.zeros(len(self.betatest))
        # for cykle for calculating integrals for different beta/epsilon
        for i in range(len(self.betatest)):
            integrand_left_side= lambda x:self.TestingIntegrandLeftSide(x, Sigma, self.betatest[i], self.epsilon)
            result_left_side, error_left_side= quad_vec(integrand_left_side, -np.inf, np.inf)
            #print(f"Integral for beta={self.beta}:",result_left_side)
            list_results_of_integral[i]=result_left_side
            #Idealizet case for beta goes to infitinty
        IdealCase=self.Resultofanalythicalbetainfinity()
        #print("Difference:=", list_results_of_integral-IdealCase)
        self.PF.PlotResultsofIntegrals(self.betatest,list_results_of_integral)
        #print(f"Case for beta goes to infinity",IdealCase )


#Here is a idealized D term for beta goes to infinity
    def Dominant_D_Result(self, Sign):
        alpha= self.Gamma_0()/(2*self.t_value)**2
        z=Sign*1j*((-alpha)-csqrt(1+alpha**2))
        FirstTerm=-(self.t_value**3/(np.pi*self.Gamma_0()**2))
        FirstRatio= (((Sign*5j*z)*(-alpha)-3)/((Sign*1j*z)*(-alpha)-1))
        SecondRatio= ((Sign*1j*z)/(csqrt(1+alpha**2)))
        Finalresult= FirstTerm*(Sign*1j*z-alpha*(FirstRatio+SecondRatio+1))
        return Finalresult

# Here is a ideal case for specific choose of sign
    def DominantDtermiinbetagoestoinfinity(self, Sign ):
        if Sign==(-1) or Sign==1:
            Dominant_D_result=self.Dominant_D_Result(Sign)
            return Dominant_D_result
        else:
            print("Only -1 or 1 values are free")


# Here is third dominant term in equation for D
    def ThirdDRation(self, x, Sigma):
        FirstPartofDenom = (-x + Sigma)
        denom = (FirstPartofDenom * csqrt(4 * self.t_value ** 2 - (-(x) + Sigma) ** 2)) ** 3
        frac = (8 * self.t_value ** 2 - 3 * FirstPartofDenom ** 2) / denom
        return 8 * (self.t_value ** 2) * frac.real

#Here is function for comparsion and print Idealized D for beta goes to infinity a nd Integral of Third D for high beta value
    def CheckDominantTerminD(self):
        Sigma= self.ZerothIterationofSelfEnergy()
        Integrand_Dominant_term_D=lambda x: (1/np.pi)*(1/16)*self.Mkf.FermiFunction(x)*self.Mkf.ThirdDRation(x,Sigma)
        Result_Dominant_term_D, Error_Dominant_term_D= quad_vec(Integrand_Dominant_term_D, -np.inf, np.inf)
        print(f"Main term integral for beta={self.beta}", Result_Dominant_term_D)
        IdealCaseDominantTerm= self.Dominant_D_Result(-1)
        print(f"Case for beta=infinity", IdealCaseDominantTerm)


