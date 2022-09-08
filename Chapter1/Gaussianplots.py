#Requirement install filterpy
#Plotting Gaussian distributions with \mu=23 and different variances \sigma^2
#GuelCortez 2022
from filterpy.stats import gaussian
import matplotlib.pyplot as plt
import numpy as np

xs = np.arange(15, 30, 0.05)
plt.plot(xs, gaussian(xs, 23, 0.2**2), label='$\sigma^2=0.2^2$')
plt.plot(xs, gaussian(xs, 23, .5**2), label='$\sigma^2=0.5^2$', ls=':')
plt.plot(xs, gaussian(xs, 23, 1**2), label='$\sigma^2=1^2$', ls='--')
plt.legend()
plt.xlabel("x")
plt.ylabel("p(x)")
plt.show()
