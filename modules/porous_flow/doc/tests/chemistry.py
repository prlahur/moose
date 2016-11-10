#!/usr/bin/env python

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc

#
# Reaction - convection - diffusion equation
#
# Read MOOSE simulation data
data = np.genfromtxt('../../tests/chemistry_advection/ad0_out_gc_0021.csv', delimiter = ',', names = True, dtype = float)

# The analytical solution is given by Lichtner, Continuum formulation of multcomponent
# multiphase reactive transport, Rev. Min. and Geochem., 1997
def analytical(x, t):
    porosity = 0.1
    ceq = 0.1
    v = 0.01
    D = 0.04
    ks = 0.1
    c0 = 1

    W = np.sqrt(1 + 4 * ks * porosity * D / v**2)
    q = v / (2 * porosity * D) * (W - 1)

    return ceq - 0.5 * (ceq - c0) * (np.exp(-q * x) * erfc((x - W * v * t / porosity) / (2 * np.sqrt(D * t))) +
      np.exp((1 + W) * (v * x / (2 * porosity * D))) * erfc((x + W * v * t / porosity) / (2 * np.sqrt(D * t))))

x = np.linspace(0, 2, 100)
t = 10

plt.figure(1)
plt.plot(x, analytical(x, t), label = 'Analytical')
plt.plot(data['x'], data['gc'], 'o', label = 'MOOSE')
plt.xlabel('x (m)')
plt.ylabel('Generalised concentration (-)')
plt.legend()
plt.title('Generalised concentration (t = 10 s)')
plt.ylim([-0.05,1])
plt.savefig("chemistry_advection_fig.pdf")

sys.exit(0)
