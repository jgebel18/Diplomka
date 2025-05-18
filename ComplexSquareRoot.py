import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.scimath import sqrt as csqrt

#Here is a testing code for confim the quastion about complex square root

def DefinineComplexSquareRoot(theta):
    z=np.cos(theta) + 1j*np.sin(theta)
    print(csqrt(1-1j))
    print(2**(0.25)*np.sqrt((1+1/np.sqrt(2))/2))
    print(2**(0.25)*np.sqrt((1-1/np.sqrt(2))/2))
    return csqrt(z)




theta_values= np.linspace(0,2*np.pi, 100)
ComplexSquareRoot= DefinineComplexSquareRoot(theta_values)

#
plt.figure(figsize=(8, 8))
plt.plot(ComplexSquareRoot.real, ComplexSquareRoot.imag, label="sqrt(z)")
plt.xlabel("Re")
plt.ylabel("Im")
plt.title("Complex Square Root on Unit Circle")
plt.grid(True)
plt.axis('equal')
plt.legend()
plt.savefig("ComplexSquareRoot.png")
plt.show()

