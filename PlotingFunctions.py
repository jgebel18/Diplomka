import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")

#Dear Sunil,
#I created this class because we had too many scattered plotting functions, so I decided to centralize them.

class PlottigFunctions:
    #Here is a function that can polt all in tegrands in different values of beta
    def PlotIntegrandsfordifferentbetavalues(self,beta,x_values, integrand_2, integrand_3):
        fig, ax = plt.subplots()
        for i in range(len(beta)):
            beta_val = beta[i]
            ax.plot(x_values, integrand_2[i], label=f"Integrand_lhs (β={beta_val})")
            ax.plot(x_values, integrand_3[i], label=f"Integrand_rhs (β={beta_val})")
        ax.set_title("Integrands for Different β Values")
        ax.set_xlabel('x')
        ax.set_ylabel('Integrand Value')
        ax.legend()
        plt.tight_layout()
        plt.savefig('Integrands_for_beta_values.png')
        plt.show()

    def PlotResultsofIntegrals(self, X,Y):
        plt.figure(figsize=(10,5))
        plt.plot(X, Y , label="LHS Results for epsilon")
        #plt.plot(X, Y_2, label="Ideal beta =inf")
        plt.xlabel('epsilon')
        plt.ylabel('Results')
        plt.legend()
        plt.savefig('Results_for_beta_values.png')
        plt.savefig('Results_for_beta_values.pdf')
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
    def PlotData(self, Nazev, dataX, dataY, Nazev_X, Nazev_Y, Title):
        plt.figure(figsize=(10, 10))
        plt.plot(dataX, dataY, label=Nazev)
        plt.xlabel(Nazev_X)
        plt.ylabel(Nazev_Y)
        plt.title(Title)
        plt.savefig(Title + 'pgn')
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
        plt.show()
