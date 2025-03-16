import numpy as np

from Root import Root_finding as rf

k=1.38e-16 # erg/K
aB = 2e-13 # cm^3 / s
psi = 0.929
Tc = 1e4
Z = 0.015
A = 5e-10
xi = 1e-15

def equilibrium2(T, nH, Z = Z, Tc = Tc, psi = psi, A = A, xi = xi):
    return (psi*Tc - (0.684 - 0.0416 * np.log(T/(1e4 * Z*Z)))*T - .54 * ( T/1e4 )**.37 * T)*k*nH*aB + A*xi + 8.9e-26 * (T/1e4)

def equilibrium2_diff(T, nH, Z=Z, psi=psi, Tc=Tc, A = A, xi= xi):

    part1 = -0.684 + 0.0416 * (1 + np.log(T/(1e4 * Z*Z) ) )
    part2 = -0.54 * 1.37 * ( T/1e4 )**.37
    part3 = 8.9e-30
    return (part1 + part2) * k * nH * aB + part3

xmin, xmax = 1, 1e15
root_finder2 = rf(xmin, xmax)

n_list = [1e-4, 1, 10**4]
for i in range(len(n_list)):
    eq2 = lambda x: equilibrium2(x, n_list[i])
    eq2_diff = lambda x: equilibrium2_diff(x, n_list[i])
    
    root_eq2 = root_finder2.newton_raphson(1000000, 1, 1e-10, eq2, eq2_diff)
    print(f' The equilibrium temperature for a density of {n_list[i]} ' + r'$\mathrm{cm^{-3}}$'  + f' is {root_eq2[0]} K with a relative error of {root_eq2[1]} K ')

