import numpy as np
from MakeFunctions import MakeFunctions
from scipy.interpolate import UnivariateSpline
from numpy.lib.scimath import sqrt as csqrt
from PlotingFunctions import PlottigFunctions
from Class_for_testing_calculations import TestingCalculations
class MakeIterationFuntions:

  #Here is a construction of class variables
    def __init__(self, beta, x_values, t_value, omega_values, U,
                 NumIterations, Tolerance, Rezolution,Gamma_0,Sigma_0):
        self.beta = beta            # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value #Hopping parameter
        self.omegavalues=omega_values #omega for Greeen function                # Hopping parameter
        self.U = U   #Electron Iteraction
        self.Mkf= MakeFunctions(self.beta, self.x_values, self.t_value, self.omegavalues,self.U) # This calling of class is important for define another iterations
        self.NumIteration=NumIterations #Num Iterations
        self.Rezolution=Rezolution #Rezolution
        self.Tolerance=Tolerance #Tolerance

        self.PF=PlottigFunctions() # Calling of class with plotiing methods
        self.Gamma_0= Gamma_0
        self.Sigma_0=Sigma_0
        self.TestingClass = TestingCalculations(self.beta, self.x_values, self.t_value, self.omegavalues, self.U,
                                                self.NumIteration, self.Tolerance, self.Rezolution, self.Gamma_0, self.Sigma_0)

    #This is a First equation for a
    def a (self, Sigma):
        #print(self.Mkf.Y(Sigma))
        a= 1.0+self.U*self.Mkf.Y(Sigma)
        print(a)
        return a

    #Another Iterations of Gamma
    def Gamma (self, Sigma):
        DotofItegrals= csqrt(2/(self.Mkf.D(Sigma)*self.U))
        Gamma= (self.U/(self.beta*np.pi))*DotofItegrals*(1/csqrt( self.a(Sigma)))
        # self.TestingClass.CheckDominantTerminD()
        return Gamma

    # Here is a Another itertions of Sigma function
    def NewSigma(self, omega, Sigma, i):
        Gamma= self.Gamma(Sigma)
        if i == 0:
            NewSigma=(-1)*Gamma*self.Mkf.G(omega, Sigma)
            return NewSigma
       # self.PF.PlotData("i-tá iterace G(ω)", self.omegavalues,self.Mkf.G(omega, Sigma), "ω", "G(ω)", "Greenova funkce " )
        else:
            NewSigma= (-1)*Gamma*self.Mkf.G(omega, Sigma)
        return NewSigma

    #Here is a function to plotting every parameter of this Iteration Process
    def PrintEveryParameter(self,i,Sigma):
        print(f"""Values of parameters in {i}-Iteration
        β = {self.beta}
        U = {self.U}
        a = {self.a(Sigma)}
        Γ = {self.Gamma(Sigma)}
        D = {self.Mkf.D(Sigma)}
        Y = {self.Mkf.Y(Sigma)}	
        Σ_{0} = {self.Sigma_0}
        Γ_{0} = {self.Gamma_0}""")


    #Main Iteration Process I have upgraded and I added an array to collect values in every different iteration
    def Iterationprocess(self):
        Sigma = self.Sigma_0
        SigmaIterations = []
        SigmaValues = np.full(len(self.omegavalues), Sigma)
        NewSigma=None
        #SigmaIterations.append(SigmaValues)
        for i in range(self.NumIteration):
            Sigma = self.InterpolateOfSigmaValues(SigmaValues, self.omegavalues)
            self.PrintEveryParameter(i,Sigma)
            SigmaIterations.append(Sigma)
            NewSigma = self.NewSigma(self.omegavalues, Sigma, i) #Calculating Sigma

            self.PF.PlotGreenFunction(self.omegavalues, self.Mkf.G, SigmaIterations, np.array(SigmaIterations).size)
            self.PF.PlotSigmaFunction(self.omegavalues, np.array(SigmaIterations))#Plotting of all iterations of sigma
            # Tolerance check
            if self.CheckingforContinue(SigmaValues,NewSigma)==True:
                return NewSigma, np.array(SigmaIterations)
            else:
                SigmaValues = NewSigma  # aktualization of sigma vlaues
        return NewSigma, np.array(SigmaIterations)

#Here is the Interpolation of SelfEnergy. I decided for start with this iterpolation function. In future
#I will probably change it
    def InterpolateOfSigmaValues(self, SigmaValues, OmegaValues,):
        # Interpolace reálné a imaginární části zvlášť
        SplineOmegaImag = UnivariateSpline(OmegaValues, np.imag(SigmaValues), s=0)
        SplineOmegaReal = UnivariateSpline(OmegaValues, np.real(SigmaValues), s=0)
        SplineOmega = lambda x: SplineOmegaReal(x) + 1j * SplineOmegaImag(x)
        return SplineOmega

#Here is a method checking for continue of comparing diff betwen actual and past iteration of Sigma
#We have a sorted caseses for imag part and real part
    def CheckingforContinue(self, SigmaValues, NewSigmaValues):
        diff_real= np.std((NewSigmaValues.real - SigmaValues.real))
        diff_imag= np.std((NewSigmaValues.imag - SigmaValues.imag))
        print(f"std_real={diff_real}")
        print(f"std_imag={diff_imag}")
        diff = np.std(np.abs(NewSigmaValues - SigmaValues))
        print(f"std={diff}")
        if diff_real < self.Tolerance and diff_imag < self.Tolerance:
            return True
        else:
            return False


    #This function Gives  Green function
    def GiveFinalG(self):
        SigmaValues, SigmaIterations = self.Iterationprocess()
        #Sigma= self.InterpolateOfSigmaValues( SigmaValues,self.omegavalues)
        #self.TestingClass.CheckLeftSideIntegrand()
        #self.TestingClass.CheckDominantTerminD()
        #self.PF.PlotSigmaFunction(self.omegavalues, SigmaIterations)
        #self.PF.PlotGreenFunction(self.omegavalues, self.Mkf.G,  SigmaIterations,SigmaIterations.size)