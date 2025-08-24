import numpy as np
import matplotlib.pyplot as plt



x_values= np.linspace(-2*np.pi, +2*np.pi, 100)


y=np.exp(1j*x_values)


plt.figure(figsize=(8,8))
plt.plot(x_values, np.abs(y))
#plt.plot(x_values, y.real)
plt.xlabel('x'
           )
plt.ylabel('y')
plt.show()