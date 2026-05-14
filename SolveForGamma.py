import numpy as np
from scipy.optimize import root_scalar,fsolve, minimize
from numericke_metody.nonlin_equations import bisection, newton, simple_iteration
from Class_for_testing_calculations import TestingCalculations
from Main_Terms_In_Equations import Main_Terms_Equations
from MakeFunction_Gamma import MakeFunctions_Gamma
from Operation_with_files import Operation_with_files
from PlotingFunctions import PlottigFunctions
import os

class MakeIterationProcessforGamma:
  #Here is a construction of class variables
    def __init__(self, beta, x_values, t_value, omega_values, U, U_values ,
                 NumIterations, Tolerance, Rezolution,Gamma_0, Sigma_0):

        self.beta = beta           # Inverse temperature
        self.x_values = x_values    # Array of x values
        self.t_value = t_value #Hopping parameter
        self.omegavalues=omega_values #omega for Greeen function                # Hopping parameter
        self.U= U #Electron interaction
        self.U_values = U_values    #Electron Iteractions array
         # This calling of class is important for define another iterations
        self.NumIteration=NumIterations #Num Iterations
        self.Rezolution=Rezolution #Rezolution
        self.Tolerance=Tolerance #Tolerance
        self.PF=PlottigFunctions() # Calling of class with plotiing methods
        self.Gamma_0=0
        self.step= 0.01
        self.num_steps=50
        self.h5py_filename= 'Data.h5'
        self.h5py_solution_filename='Solutions.h5'
        self.Path_Files='Files'
        self.Critial_Data_File_Name= 'Data_4.txt'
        self.Path_Images='Images'
        self.Sigma_0=Sigma_0
        self.Gamma_values= np.linspace(-2,-0.001, self.num_steps+1)
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


    #This is for findin key interval for applying the nonlinear equation solver
    def SpecificateInterval(self, f, x0, x1):
        step = 100
        max_kroku = 500  # Pojistka proti nekonečné smyčce

        # 1. Rovnou máme štěstí a interval jsme trefili
        if f(x0) * f(x1) <= 0:
            return min(x0, x1), max(x0, x1)

        pocitadlo = 0

        # 2. Obě hodnoty jsou záporné
        if f(x0) < 0 and f(x1) < 0:
            if abs(f(x1)) > abs(f(x0)):
                while f(x1) < 0 and pocitadlo < max_kroku:
                    x1 += step
                    pocitadlo += 1
            else:
                while f(x0) < 0 and pocitadlo < max_kroku:
                    x0 -= step
                    pocitadlo += 1

        # 3. Obě hodnoty jsou kladné
        else:
            if abs(f(x1)) > abs(f(x0)):
                while f(x0) > 0 and pocitadlo < max_kroku:
                    x0 -= step
                    pocitadlo += 1
            else:
                while f(x1) > 0 and pocitadlo < max_kroku:
                    x1 += step
                    pocitadlo += 1

        # Kontrola, jestli jsme nevyčerpali všechny kroky naslepo
        if pocitadlo >= max_kroku:
            print("Varování: Interval nenalezen, narazili jsme na limit kroků.")
            return None  # Nebo tady můžeš vyhodit výjimku (raise Exception)

        # Vrátíme hezky seřazený interval (menší, větší)
        return min(x0, x1), max(x0, x1)



# I modified this function to gives 3D numpy.array for every important termi in solved equations
  #
    def GiveFinalG(self):
            # Inicializace jednoho 3D pole: (4 metriky, počet_beta, počet_gamma)
        Ow_files= Operation_with_files(self.h5py_filename)
        cals= np.array(['Y(Γ)', 'D(Γ)', 'a(Γ)', 'a_Nonliear(Γ)'])
        res = np.empty((4, self.beta.size, self.Gamma_values.size))
        #self.PF.Plot_Fermi_function(self.x_values, self.beta)
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
                res[3, i, j] = Terms.NonlinearEquationa(gamma, beta, res[0,i,j], res[1,i,j]/gamma**2)
                #print(res[3,i,j])
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
        self.PF.Plot_Beta_Gamma_Dependence(res, self.beta, self.Gamma_values, self.U)
        return res

    #This is for  solution analyzing and read the h5py-file and plots in
  # graphs with limits of the intervals which are imput parameter to Nonlinear equation solver
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
        self.PF.Plot_Beta_Gamma_Dependence(Result_values, self.beta, self.Gamma_values, values_complete , self.U)
        return x0, x1, x1_new

    #For solving nonlinear equations using interval methods we needed to find specific interval limits
  # values x0,x1 to satisfy this condition
  # f(x0)f(x1)<0 and before using this interval method i decided to use correctly
  # directed Monte Carlo algorithm
    def GenerateIntervalMonteCarlo(self, function, search_interval=(-100, 0)):
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


    #This method is so far unique particaly functable final nonlinear equation solver
    def SolveFinalEquation(self):
        solutions = np.empty((self.beta.size, self.U_values.size), )
        current_U = {}
        Op_files = Operation_with_files(self.h5py_filename)
        # data=Op_files.Read_file('solutions')
        # for i,  beta  in enumerate(data):
        # print(f'Solutions-for-beta={self.beta[i]}:', data[i])
        tol = 1e-24
        # self.PF.Plot_Fermi_function(self.x_values, self.beta)
        for i, beta in enumerate(self.beta):
            for j, U in enumerate(self.U_values):
                mkf = MakeFunctions_Gamma(beta, self.x_values, self.t_value, self.omegavalues, U)
                terms = Main_Terms_Equations(beta, self.x_values, self.t_value, self.omegavalues, U,
                                             self.NumIteration, self.Tolerance, self.Rezolution)
                x0, x1 = -10, -1e-90  # self.GenerateIntervalMonteCarlo(function)
                f = lambda gamma: terms.NonlinearEquationA(gamma, beta, mkf)
                #print((f(x0),f(x1)))
                if f(x0) * f(x1) < 0:
                    solution_1 = root_scalar(f=terms.NonlinearEquationA, args=(beta, mkf), method='toms748',
                                             maxiter=1500, xtol=tol, rtol=2.22045e-16, bracket=(x0, x1))
                    # print(solution_1.root)
                    if np.abs(solution_1.root) > np.abs(tol):
                        solutions[i][j] = solution_1.root
                        print(f'a-solution-for-U={U}-beta={beta}', terms.a(solutions[i][j] , mkf)* beta ** 2)
                        print(f'Gamma-solution-U={U}-beta={beta}', solution_1.root)
                        print(f'U={U}-beta={beta}', solution_1.iterations)
                    else:
                        solutions[i][j] = np.nan
                else:
                    solutions[i][j] = np.nan
                current_U[f'for-U:{U}-beta-{beta}'] = solutions[i][j]
            print(f'Solutions-for-beta={beta}:', solutions[i])
            # self.PF.PlotGammaNonlinear(self.Gamma_values, self.GammaNonlinearFunction(self.Gamma_values,beta,mkf),beta)
        # self.PF.Plot_3D_solutions(self.beta, self.U_values, solutions)
        Op_files.Write_to_file(f'solutions-{(self.U_values.min, self.U_values.max)}-{(self.beta.min, self.beta.max)}', solutions)

    #I made this function to confim Šimon's theoretical concepts and plotting these integrals
  # and printing the a-function and this aproximation
    def PlotIntergrals(self, ):
        filtered = []
        Gamma_2 = np.linspace(-5, -0.1, 50)
        Mkf = MakeFunctions_Gamma(10, self.x_values, self.t_value, self.omegavalues, self.U)
        #Sigma= Mkf.Sigma(self.Gamma_0)
        for i, gamma in enumerate(Gamma_2):
            integral_Y= Mkf.Y(gamma)
            integral_D= Mkf.D(gamma)
            integral_Y_approx = Mkf.Y_Gamma_approx_integrand(gamma)
            integral_D_approx = Mkf.D_Gamma_approx_integrand(gamma)
            #Sigma = Mkf.Sigma(gamma)
            print(f'Integrals Y and D for gamma {gamma}', np.array([integral_Y, integral_D]))
            print(f'Integrals Y_aprox and D_approx for gamma {gamma}', np.array([integral_Y_approx, integral_D_approx]))
            #Mkf.CollectDataAndPlot(gamma)

    #I added this function for read gained data and get this data to plot for šimon


    def ReadFilewithData(self):
        filtered =  []
        a_critical_beta = np.zeros((len(self.beta), len(self.U_values)))
        roots= np.empty(a_critical_beta.shape)#a_critical_beta.copy()
        name = os.path.join(self.Path_Files, self.Critial_Data_File_Name)
        with open(name, "r", encoding="utf-8") as f:
            lines = f.readlines()
            for i, beta in enumerate(self.beta):
                for j, U in enumerate(self.U_values):
                    key = f'a-solution-for-U={U}-beta={beta}'
                    key_gamma = f'Gamma-solution-U={U}-beta={beta}'
                    for line in lines:
                        if line.startswith(key):
                            value = float(line.split()[-1])
                            a_critical_beta[i, j] = value
                        elif line.startswith(key_gamma):
                            #print('Gamma_found')
                            value=float(line.split()[-1])
                            roots[i,j]= value

            for i in range(self.beta.size):
                constant= roots[6][-1]
                condition=  (roots[i]!= constant)
                filtered_row = a_critical_beta[i][condition]
                filtered.append(filtered_row)


        self.PF.Plot_Final_Dependence(
            self.beta,
            self.U_values,
             filtered,
        )












