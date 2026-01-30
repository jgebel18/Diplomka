import matplotlib.pyplot as plt
import matplotlib
import numpy as np
matplotlib.use("TkAgg")

#Dear Sunil,
#I created this class because we had too many scattered plotting functions, so I decided to centralize them.

class PlottigFunctions:
    #Ploting Integrands for Values of Sigma,
    def PlotingItegrands(self, D_array, Y_array ,x_values, ):
        fig , axes = plt.subplots(2,1, figsize=(10,10))
        for i in range(D_array.shape[0]):
            axes[0].plot(x_values, D_array[i], label=f"{i}-Iteration")
            axes[1].plot(x_values, Y_array[i], label=f"{i}-Iteration")
            axes[0].set_title("D-Integrands-Iterations")
            axes[1].set_title("Y-Integrands-Iterations")
            axes[0].set_xlabel('x')
            axes[1].set_xlabel('x')
            axes[0].set_ylabel('D-Itegrands')
            axes[1].set_ylabel('Y-Itegrands')
            axes[0].legend()
            axes[1].legend()
            axes[0].grid(True)
            axes[1].grid(True)
        plt.savefig(f'Integrands_for_Interpolated_Values.png')
        plt.show()


    #Here is a function that can polt all in tegrands in different values of beta
    def PlotIntegrandsfordifferentbetavalues(self,beta,x_values, integrand_2, integrand_3):
        fig, ax = plt.subplots()
        for i in range(len(beta)):
            beta_val = beta[i]
            ax.plot(x_values, integrand_2[i], label=f"Integrand_lhs (β={beta_val})")
            ax.plot(x_values, integrand_3[i], label=f"Integrand_rhs (β={beta_val})")
        ax.set_title("Integrands for Different eslion Values")
        ax.set_xlabel('ε')
        ax.set_ylabel('Integrand Value')
        ax.legend()
        plt.tight_layout()
        plt.savefig('Integrands_for_beta_values.png')
        plt.show()

    def PlotResultsofIntegrals(self, X,Y,Y_2):
        plt.figure(figsize=(10,5))
        plt.plot(X, Y , label="D Integral Results for     different epsilon values")
        plt.plot(X, Y_2, label="Ideal beta =inf")
        plt.xlabel('ε')
        plt.ylabel('Integral results')
        plt.legend()
        plt.savefig('Results_for_epsilon_values.png')
        plt.show()



    #I added here function to plot all integrnds in one picture and every part which is important to define
    # D function
    def PlotIntegrands(self, nazev, DataX, DataY, nazev1, DataY2, nazev2, DataY3, nazev3):
        fig, axes = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
        # První subplot
        axes[0].plot(DataX, DataY, label=nazev1)
        axes[0].set_title(nazev1)
        axes[0].set_ylabel('Integrand')
        axes[0].legend()

        # Druhý subplot
        axes[1].plot(DataX, DataY2, label=nazev2, color='orange')
        axes[1].set_title(nazev2)
        axes[1].set_ylabel('Integrand')
        axes[1].legend()

        # Třetí subplot
        axes[2].plot(DataX, DataY3, label=nazev3, color='green')
        axes[2].set_title(nazev3)
        axes[2].set_xlabel('x')
        axes[2].set_ylabel('Integrand')
        axes[2].legend()
        fig.suptitle(nazev)
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # aby nebyl překryt hlavní titulek
        plt.savefig(nazev + '_subplots.png')
        plt.show()

    # Here is another plotting function which is ploting one plot in one picture
    def PlotData(self, dataX, dataY):
        plt.figure(figsize=(10, 10))
        plt.plot(dataX, dataY )
        plt.xlabel('ω')
        plt.ylabel('G(ω)')
        plt.savefig('Integrands.png')
        plt.show()

    def PlotTwoSubplots(self, Nazev1, Nazev2, dataX, dataY1, dataY2, Nazev_X, Nazev_Y, Title):
        plt.figure(figsize=(12, 6))  # širší formát na dva podgrafy vedle sebe

        # První subplot (levý)
        plt.subplot(1, 2, 1)
        plt.plot(dataX, dataY1, label=Nazev1, color='blue')
        plt.xlabel(Nazev_X)
        plt.ylabel(Nazev_Y)
        plt.title(Nazev1)
        plt.grid(True)

        # Druhý subplot (pravý)
        plt.subplot(1, 2, 2)
        plt.plot(dataX, dataY2, label=Nazev2, color='red')
        plt.xlabel(Nazev_X)
        plt.ylabel(Nazev_Y)
        plt.title(Nazev2)
        plt.grid(True)
        plt.suptitle(Title)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # místo pro hlavní titulek
        plt.savefig(Title + '.png')  # oprava přípony
        plt.show()

    # Here I added ploting metod for every iteration  of Green function for the future
    # I will add the methods to plots every integrands of all functions
    def PlotGreenFunction(self, omega, Function, Sigma, NumberIterationToPlot):
        fig, ax = plt.subplots()
        for i in range(NumberIterationToPlot):
            data= Function(omega, Sigma[i])
            ax.plot(omega, data.real, label=f"{i}-iterace funkce G(ω) reálná část")
            ax.plot(omega,data.imag, label=f"{i}-iterace funkce G(ω) imaginární část" )
        ax.set_xlabel("ω")
        ax.set_ylabel("G(ω)")
        ax.legend()
        plt.grid(True)
        plt.savefig('Greens_function_omega.png')
        plt.show()

 #This function plots every shape of Sigma Interpolated and Original Sigma values
    def PlotSigmaFunction(self, omega, SigmaInterpolated, SigmaOriginal, ):
        fig, ax = plt.subplots()
        newomega= np.linspace(omega.min(), omega.max(), num=100000) # This is our different rande to compare quality of interpolation
        ax.scatter(newomega, SigmaInterpolated(newomega).real,label=f"Interpolated  Σ(ω) reálná část")
        ax.scatter(newomega, SigmaInterpolated(newomega).imag,label=f"Interpolated Σ(ω) imaginární část")
        ax.scatter(omega, SigmaOriginal.real, label=f"Original  Σ(ω) reálná část")
        ax.scatter(omega, SigmaOriginal.imag, label=f"Original  Σ(ω) imaginární část")
        ax.set_xlabel("ω")
        ax.set_ylabel("Σ(ω)")
        ax.legend()
        plt.savefig("Interpolated_Sigma.png")
        plt.show()

# This is function for plting every root of quartic equation in the dependence of omega

    def PlotRootsofQuarticEquation(self, roots, omega):
        nazev = 'Self-Energy for every root of quartic equation'
        fig, axes = plt.subplots(2, 2, figsize=(10, 12), sharex=True)
        axes = axes.flatten()  # převede [ [a,b], [c,d] ] -> [a,b,c,d]

        for i in range(len(roots)):
            axes[i].scatter(omega, roots[i].real, label=f"{i}-th Root (Re)")
            axes[i].scatter(omega, roots[i].imag, label=f"{i}-th Root (Im)")
            axes[i].set_xlabel('ω')
            axes[i].set_ylabel('Σ(ω)')
            axes[i].set_title(f'{i}-th Root of Equation')
            axes[i].legend()
            axes[i].grid(True)

        fig.suptitle(nazev)
        plt.tight_layout(rect=[0, 0, 1, 0.96])  # aby nebyl překryt hlavní titulek
        plt.savefig(nazev + '_subplots.png', dpi=300)
        plt.show()

#This is modified function to plting the roots of both equation Nonlinear and Quartic
    def PlotRootsOfEquation(self, roots, omega):
        nazev = 'Self-Energy for every root of both equation'
        fig, axes = plt.subplots(3, 2, figsize=(10, 12), sharex=True)
        axes = axes.flatten()  # převede [ [a,b], [c,d] ] -> [a,b,c,d]
        for i in range(len(roots)):
            if i == 0:
                axes[i].scatter(omega, roots[i].real, label=f"Nonlinear equation Root (Re)")
                axes[i].scatter(omega, roots[i].imag, label=f"Nonlinear equation Root (Im)")
                axes[i].legend()
                axes[i].grid(True)
                axes[i].set_xlabel('ω')
                axes[i].set_ylabel('Σ(ω)')
            else:
                axes[i].scatter(omega, roots[i].real, label=f"{i}-th Root (Re)")
                axes[i].scatter(omega, roots[i].imag, label=f"{i}-th Root (Im)")
                axes[i].legend()
                axes[i].grid(True)
                axes[i].set_xlabel('ω')
                axes[i].set_ylabel('Σ(ω)')
        fig.suptitle(nazev)

        plt.tight_layout(rect=[0, 0, 1, 0.96])  # aby nebyl překryt hlavní titulek
        plt.savefig(nazev + '_subplots.png', dpi=300)
        plt.show()

    def PlotBothpartsofSigma(self, Sigma_values, omega):
        plt.figure(figsize=(10, 10))
        plt.scatter(omega, Sigma_values.real, label='Real Part')
        plt.scatter(omega,Sigma_values.imag, label='Imaginary Part')
        plt.legend()
        plt.grid(True)
        plt.xlabel('ω')
        plt.ylabel('Σ(ω)')
        plt.title('Complete Roots of Quartic equation')
        plt.show()

    def PlotRootOfEquation(self, root, omega):
        nazev = 'Self-Energy for root of nonlinear equation'
        plt.figure( figsize=(10, 12))
        plt.scatter(omega, root.real, label=f" Root (Re)")
        plt.scatter(omega, root.imag, label=f" Root (Im)")
        plt.xlabel('ω')
        plt.ylabel('Σ(ω)')
        plt.title(f'Root of Equation')
        plt.legend()
        plt.grid(True)
        plt.title(nazev)
        #plt.tight_layout(rect=[0, 0, 1, 0.96])  # aby nebyl překryt hlavní titulek

        plt.savefig(nazev + '_subplots.png', dpi=300)
        plt.show()

    def Plot_Values_of_a_and_D(self, values_a, Values_D,D_approx,
                               Gamma_Values):
        len = 5
        nazev='Values_of_Y_D.png'
        D_control= Values_D*Gamma_Values**2
        plt.figure(figsize=(10, 10))
        #plt.scatter(Gamma_Values, values_a, label= 'Values of Y(0,0) ')
        #plt.scatter(Gamma_Values, Values_D, label= 'Values of D')
        plt.scatter(Gamma_Values, D_approx ,label='D approx')
        #plt.plot(Gamma_Values[len], values_a[len], 'rv', label='MidPoint for a ')
        #plt.plot(Gamma_Values[len], Values_D[len], 'gv', label= 'Midpoint for D')
        plt.scatter(Gamma_Values, D_control, label= 'D control')
        plt.grid(True)
        plt.legend()
        plt.xlabel('Gamma')
        plt.ylabel('Values of Y and D ')
        plt.savefig(nazev)
        plt.show()
