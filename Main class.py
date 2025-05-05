import numpy as np
import matplotlib.pyplot as plt
from MakeFunctions import MakeFunctions
import matplotlib
matplotlib.use('TkAgg')
# Simulation parameters
t = 1.0
U = 4.0
beta = 500
x_values = np.linspace(-100, 100, 100)

# Class initialization
MkF = MakeFunctions(beta, x_values, t, U)

# Calculate the integrand (real part of complex function)
y = np.array([(1/np.pi) * MkF.FermiFunction(x) * MkF.ComplexFrac_in_Y(x) for x in x_values])

# Calculate the function a = 1 - U * Y(0,0)
a = 1 - U * MkF.Y()[0]  # [0] is the result of integration, [1] is the error estimate
print("a = ", a)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_values, y)
plt.title("Complex integrand (real part)")
plt.xlabel("x")
plt.ylabel("Integrand")
plt.grid(True)
plt.tight_layout()
plt.savefig("integrand_plot_different_Sigm.png")
plt.show()
