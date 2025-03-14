import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar as rs

k=1.38e-16 # erg/K
aB = 2e-13 # cm^3 / s
psi = 0.929
Tc = 1e4
Z = 0.015
A = 5e-10
xi = 1e-15

def equilibrium1(T, Z = Z, Tc = Tc, psi=psi):
    return psi*Tc*k - (0.684 - 0.0416 * np.log(T/(1e4 * Z*Z)))*T*k

def equilibrium2(T, nH, Z = Z, Tc = Tc, psi = psi, A = A, xi = xi):
    return (psi*Tc - (0.684 - 0.0416 * np.log(T/(1e4 * Z*Z))) - .54 * ( T/1e4 )**.37 * T)*k*nH*aB + A*xi + 8.9e-26 * (T/1e4)

def equilibrium2_diff(T, nH, Z=Z, psi=psi, Tc=Tc, A = A, xi= xi): 
    part1 = -0.684 + 0.0416 * (1 + np.log(T/(1e4 * Z*Z) ) ) * k * nH * aB
    part2 = -0.54 * 1.37 * ( T/1e4 )**.37 * k * nH * aB
    part3 = 8.9e-30

    return part1 + part2 + part3

def bisection(xmin, xmax, function): 
    a, b = xmin, xmax
    
    iterations = int(np.log2( (b-a) / 0.0001))
    for i in range(iterations):
        c = (xmin + xmax)*0.5
        if function(a)*function(c) < 0: 
            xmax = c
        else:
            xmin = c
    return c

def secant(a, b, n, function):
    for i in range(n):
        b, a = b - (b -a) / (function(b) - function(a)) * function(b), b
        # print(np.abs(a-b)/b, np.abs(a-b))
        
    return a, np.abs(b-a)/b
            
def false_position(a, b, n, target, function): 
    
    best = 0
    for i in range(n): 
        new = b - (b -a) / (function(b) - function(a)) * function(b)
        if function(new)*function(a) < 0: 
            b  = new
        elif function(new) * function(b) < 0: 
            a = new   
            
        error = np.abs(best - new)
        best = new

        if error < target: 
            break
        
    return best, error

def Newton_Raphsen(x, n, target, function, diff_function): 
    best = x
    for i in range(n): 
        # print(i)
        new = best - function(best) / diff_function(best)
        # print(new)
        error = np.abs(best - new)
        # print(error)
        
        if error < target: 
            break
        best = new

        
    return best, error
        
def func_root(x): 
    return x**2 - 5


xmin, xmax = 1, 1e7

# root_Newton2 = Newton_Raphsen((xmax-xmin)*0.5, 100, 1e-12, 1e-4, equilibrium2, equilibrium2_diff)

array = np.arange(xmin, xmax)
plt.plot(array, equilibrium2(array, 1e-4))
plt.yscale('log')
plt.show()

root_eq1, error_eq1 = false_position(xmin, xmax, 100, 0.1, equilibrium1)
root_scipy = rs(equilibrium1, bracket = [xmin, xmax])

xmin, xmax = 1, 1e15

n_list = [1e-4, 1, 10**4]
for i in range(len(n_list)):
    eq2 = lambda x: equilibrium2(x, n_list[i]) * 1e15
    eq2_diff = lambda x: equilibrium2_diff(x, n_list[i])
    
    # false_pos_guess, error_guess = false_position(xmin, xmax, 1000, 1e-7, eq2)
    root_eq2 = Newton_Raphsen(5e14, 100, 1e-12, eq2, eq2_diff)

# root_eq2 = Newton_Raphsen((xmax-xmin)*0.5, 100, 1e-12, eq2, eq2_diff)
# root_scipy2 = rs(eq2, x0 = (xmax-xmin)*0.5, fprime = eq2_diff, method = 'newton')

# plt.axvline(root_eq1)
# plt.plot(array, equilibrium1(array, 1e-4))
# plt.yscale('log')
# plt.show()
