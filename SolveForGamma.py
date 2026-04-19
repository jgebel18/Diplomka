import numpy as np
from scipy.optimize import root_scalar,fsolve, minimize
from numericke_metody.nonlin_equations import bisection, newton, simple_iteration
from Class_for_testing_calculations import TestingCalculations
from Main_Terms_In_Equations import Main_Terms_Equations
from MakeFunction_Gamma import MakeFunctions_Gamma
from Operation_with_files import Operation_with_files
from PlotingFunctions import PlottigFunctions


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
        self.h5py_filename= 'Data.h5'
        self.Path_Files='Files'
        self.Path_Images='Images'
        self.Sigma_0=Sigma_0
        self.Gamma_values= np.linspace(-1, -0.30, self.num_steps+1)
        #self.Gamma_values_right= np.linspace(self.Gamma_0, self.Gamma_0+self.num_steps*self.step, self.num_steps+1)
        #np.unique(np.concatenate((self.Gamma_values_left,                    #                         self.Gamma_values_right)))
        self.TestingClass = TestingCalculations(self.beta, self.x_values, self.t_value, self.omegavalues, self.U,
                                                self.NumIteration, self.Tolerance,
                                                self.Rezolution, self.Gamma_0,
                                                self.Sigma_0)
        self.D2_values_Sigma= None

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
        Ow_files= Operation_with_files(self.h5py_filename)
        cals= np.array(['Y(Γ)', 'D(Γ)', 'a(Γ)', 'a_Nonliear(Γ)'])
        res = np.empty((4, self.beta.size, self.Gamma_values.size))
        for i, beta in enumerate(self.beta):
            Terms= Main_Terms_Equations(self.beta, self.x_values, self.t_value, self.omegavalues, self.U,
                                        self.NumIteration, self.Tolerance, self.Rezolution)
            mkf = MakeFunctions_Gamma(beta, self.x_values, self.t_value, self.omegavalues, self.U)
            for j, gamma in enumerate(self.Gamma_values):
                res[0, i, j] = mkf.Y(gamma)
                res[1, i, j] = mkf.D(gamma)*gamma**2
                # This is for saving the calculation times I decide to calculate Y
                # and D once I also created special functions in this class
                res[2, i, j] = 1+self.U*res[0, i, j]#self.Gamma_function(gamma, beta, mkf)
                res[3, i, j] = Terms.NonlinearEquationa(gamma, beta, res[2,i,j], res[1,i,j]/gamma**2)
        for i, nazev in enumerate(cals):
            Ow_files.Write_to_file(nazev, res[i] )
        #self.PF.Plot_Beta_Gamma_Dependence(res, self.beta, self.Gamma_values)
        #self.PF.Plot_Values_of_a_and_D(Y_values, D_values,Y_approx,D_approx,
        #                               self.Gamma_values, self.beta)


        #SigmaValues, SigmaIterations = self.Iterationprocess()
        # Sigma= self.InterpolateOfSiValues( SigmaValues,self.omegavalues)
        # self.TestingClass.CheckLeftSideIntegrand()
        # self.TestingClass.CheckDominantTerminD()
        # self.PF.PlotSigmaFunction(self.omegavalues, SigmaIterations)
        # self.PF.PlotGreenFunction(selfgma.omegavalues, self.Mkf.G,  SigmaIterations,SigmaIterations.size)

    #This method read required datasets from h5py file and plots it into graph.
    def Read_Data_and_plot(self):
        names= np.array(['Y(Γ)', 'D(Γ)', 'a(Γ)', 'a_Nonliear(Γ)'])
        res= np.empty((4, self.beta.size, self.Gamma_values.size))
        Ow_files= Operation_with_files(self.h5py_filename)
        for i , name in enumerate(names):
            res[i]= Ow_files.Read_file(name)
        self.PF.Plot_Beta_Gamma_Dependence(res, self.beta, self.Gamma_values)
        return res


    def AnalyzeTheSolution_read(self):
        x0, x1, x1_new= (np.full(self.beta.size,self.Gamma_values[0], dtype=float),
                         np.full((self.beta.size),self.Gamma_values[-1], dtype=float),
                         np.empty(self.beta.size))
        for i, beta in enumerate(self.beta):
            mkf=MakeFunctions_Gamma(beta, self.x_values, self.t_value, self.omegavalues, self.U)
            terms = Main_Terms_Equations(self.beta, self.x_values, self.t_value, self.omegavalues, self.U,
                                        self.NumIteration, self.Tolerance, self.Rezolution)
            solution_1 =root_scalar(f=terms.NonlinearEquationA,args=(beta ,mkf), method='brenth',maxiter=1000, xtol=1e-4, bracket=(x0[i],x1[i]))
            x1_new[i]= solution_1.root
        Ow_files= Operation_with_files(self.h5py_filename)
        Key_Values_a, Key_Values_gamma= Ow_files.Read_file('a_limit_values'), Ow_files.Read_file('Gamma_limit_values')
        values_complete = np.vstack((Key_Values_a,  Key_Values_gamma, x0, x1, x1_new))
        Result_values= self.Read_Data_and_plot()
        self.PF.Plot_Beta_Gamma_Dependence(Result_values, self.beta, self.Gamma_values, values_complete)
        return x0, x1, x1_new

    #For solving nonlinear equations using interval methods we needed to find specific interval limits
  # values x0,x1 to satisfy this condition
  # f(x0)f(x1)<0 and before using this interval method i decided to use correctly
  # directed Monte Carlo algorithm
    def GenerateIntervalMonteCarlo(self, function, search_interval=(-100, 0), limit=10):
        Num_Iteration=int(1e4)
        for _ in range(Num_Iteration):
            x0,x1=np.random.uniform(search_interval[0],
                                    search_interval[1], 2)
            print(f'čísla',(x0,x1))
            if x0>x1:
                x0,x1=x1,x0
            if (function(x0)*function(x1)<0 and
                    np.abs(x0-x1)<=50):
                print(f'Juchůůů, čísla nalezena')
                return x0,x1
        print('Chyba, čísla nenalezena')
        return None



    #I added this this method to make calculation and solving
  # final nonlinear equation
  # and write this data into h5py files
    def AnalyzeTheSolution_write(self):
        x0, x1, x2= (np.full(self.beta.size,self.Gamma_values[0], dtype=float),
                         np.full((self.beta.size),self.Gamma_values[-1], dtype=float),
                         np.empty(self.beta.size))
        #solutions= np.empty((self.beta.size))
        Key_Values_a= np.empty((3, self.beta.size))
        Key_Values_gamma= np.empty((3, self.beta.size))
        for i, beta in enumerate(self.beta):

            mkf=MakeFunctions_Gamma(beta, self.x_values, self.t_value, self.omegavalues, self.U)
            terms = Main_Terms_Equations(self.beta, self.x_values, self.t_value, self.omegavalues, self.U,
                                        self.NumIteration, self.Tolerance, self.Rezolution)

            solution_1 =root_scalar(f=terms.a,args=( mkf), method='brenth',maxiter=1000, xtol=1e-4, bracket=(x0[i],x1[i]))
            x2[i]= solution_1.root
            print(solution_1.root)
            Key_Values_a.T[i]= np.array([terms.a(x0[i],mkf), terms.a(x1[i],mkf), terms.a(x2[i],mkf)])
            Key_Values_gamma.T[i]= np.array([terms.Gamma_function(x0[i], beta, mkf),terms.Gamma_function(x1[i], beta, mkf),
                                        terms.Gamma_function(x2[i], beta, mkf)])
        values_complete = np.vstack((Key_Values_a,  Key_Values_gamma, x0, x1, x2))
        Ow_files= Operation_with_files(self.h5py_filename)
        Ow_files.Write_to_file('a_limit_values', Key_Values_a)
        Ow_files.Write_to_file('Gamma_limit_values', Key_Values_gamma)
        Result_values= self.Read_Data_and_plot()
        self.PF.Plot_Beta_Gamma_Dependence(Result_values, self.beta, self.Gamma_values, values_complete)
        return x0, x1, x2


    #THis method is so far unique particaly functable final nonlinear equation solver
    def SolveFinalEquation(self):
        solutions= np.empty((self.beta.size), )

        for i, beta in enumerate(self.beta):

            mkf=MakeFunctions_Gamma(beta, self.x_values, self.t_value, self.omegavalues, self.U)

            terms = Main_Terms_Equations(self.beta, self.x_values, self.t_value, self.omegavalues, self.U,
                                         self.NumIteration, self.Tolerance, self.Rezolution)
            function= lambda gamma : terms.NonlinearEquationA(gamma, beta,mkf)
            x0,x1=self.GenerateIntervalMonteCarlo(function)
            solution_1 =root_scalar(f=terms.NonlinearEquationA,args=( beta, mkf), method='toms748',
                                    maxiter=1000, xtol=1e-4, bracket=(x0,x1))
            solutions[i]= solution_1.root
            #self.PF.PlotGammaNonlinear(self.Gamma_values, self.GammaNonlinearFunction(self.Gamma_values,beta,mkf),beta)
            print(f'Řešení',solutions[i])
        return solutions