import numpy as np
from MakeFunctions import MakeFunctions
from numpy.lib.scimath import sqrt as csqrt
import matplotlib.pyplot as plt
class MakeIterationFuntions:


  #Here is a construction of class variables
    def __init__(self, beta, x_values, t_value, omega_values, U, NumIterations, Tolerance, Rezolution):
        self.beta = beta            # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value #Hopping parameter
        self.omegavalues=omega_values            #omega for Greeen function                # Hopping parameter
        self.U = U   #ElectronIteraction
        self.Mkf= MakeFunctions(self.beta, self.x_values, self.t_value, self.omegavalues,self.U) # This calling of class is important for define another iterations
        self.NumIteration=NumIterations
        self.Rezolution=Rezolution
        self.Tolerance=Tolerance


    #This is a First equation for a
    def a (self, Sigma):
        #print(self.Mkf.Y(Sigma))
        a= 1.0+self.Mkf.Y(Sigma)
        return a

    #Another Iterations of Gamma
    def Gamma (self, Sigma):
        DotofItegrals= np.sqrt(2/(self.Mkf.D(Sigma)*self.U))
        Gamma= (self.U/(self.beta*np.pi))*DotofItegrals*(1/np.sqrt(self.a(Sigma)))
        return Gamma

    # Here is a Another itertions of Sigma function
    def NewSigma(self, omega, Sigma):
        Gamma= self.Gamma(Sigma)
        NewSigma=-self.Mkf.G(omega, Sigma)
        self.Mkf.PlotData("i-tá iterace", self.omegavalues,self.Mkf.G(omega, Sigma), "ω", "G(ω)", "Greenova funkce " )
        return NewSigma

    # First iteration of self-energy (complex)
    def FirstIterationofSelfEnergy(self):
        Sigma= 1j * self.Gamma_0() / (2 * self.t_value)
        SigmaValues= np.array([Sigma +i*0 for i in range (self.Rezolution)])
        return Sigma

        # Zeroth iteration of Gamma
    def Gamma_0(self):
        parameter = (2 * np.pi * self.t_value) / self.U
        gamma_0 = -(2 * self.t_value) ** 2 * np.exp(parameter)
        print(gamma_0)
        return gamma_0

    #Main Iteration Process I have upgraded and I added an array to collect values in every different iteration
    def Iterationprocess(self):
        Sigma=self.FirstIterationofSelfEnergy()
        SigmaIterations= []
        SigmaValues= np.full(len(self.omegavalues) , Sigma)
        SigmaIterations.append(SigmaValues)
        for i in range(self.NumIteration):
                NewSigma=self.NewSigma(self.omegavalues, Sigma)
                # Tolerance check:
                if np.std(np.abs(NewSigma - SigmaValues)) < self.Tolerance:
                    return SigmaValues, np.array(SigmaIterations)
                else:
                    SigmaValues = NewSigma,
        #returfinal Sigma value array values
        return Sigma, np.array(SigmaIterations)


    #Here I added ploting metod for every iteration  of Green function for the future
    # I will add the methods to plots every integrands of all functions
    def PlotGreenFunction(self, omega, Function, Sigma, NumberIterationToPlot):
        fig, ax = plt.subplots()
        for i in range(NumberIterationToPlot):
            ax.plot(omega, Function(omega, Sigma[i]),label= f"{i}-tá iterace funkce G(ω)"  )
        ax.set_xlabel("ω")
        ax.set_ylabel("G(ω)")
        ax.legend()
        plt.show()


    #This function Gives  Green function
    def GiveFinalG(self):
        Sigma, SigmaIterations = self.Iterationprocess()
        self.PlotGreenFunction(Sigma, self.Mkf.G,  SigmaIterations,SigmaIterations.size)






















