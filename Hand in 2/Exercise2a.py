import numpy as np
from Root import Root_finding as rf

k=1.38e-16 # erg/K
aB = 2e-13 # cm^3 / s
psi = 0.929
Tc = 1e4
Z = 0.015
A = 5e-10
xi = 1e-15


def equilibrium1(T, Z = Z, Tc = Tc, psi=psi):
    return psi*Tc*k - (0.684 - 0.0416 * np.log(T/(1e4 * Z*Z)))*T*k

xmin, xmax = 1, 1e7 # Bracket
root_finder = rf(xmin, xmax)

# Calculate the root using False position
root_eq1 = root_finder.false_position(100, 1, 0.1, equilibrium1)


print(f' The equilibrium temperature is {root_eq1[0]} K with a relative error of {root_eq1[1]} K ')
