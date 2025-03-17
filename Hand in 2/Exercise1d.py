import numpy as np
import Integration

A= 1. # to be computed
Nsat=100
a=2.4
b=0.25
c=1.6

def func_integrate(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    return (A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)) * x**2 * 4 * np.pi

xmin, xmax = 10**-4, 5

integration_init = Integration.Integrator(func_integrate, xmin, xmax)
result, result_error = integration_init.romberg(15)
A = 1 / result[0]

def n(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    return A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)

def analytical_diff(x, A=A, Nsat=Nsat, a=a, b=b, c=c): 
    return ((A * Nsat * b**3 * (x/b)**a * np.exp(-(x/b)**c) * (a - c * (x/b)**c - 3))* x**-4)
    
def central_differentation(x, h, function): 
    if not callable(function): 
        raise ValueError
    
    df = ( function(x+h) - function(x-h) )/ (2*h) 
    return df

def Ridder(x, h, d, e, max_iters, function):

    Ridder_matrix = np.zeros([max_iters, len(x)])
    inv_d = 1/d  
    
    for i in range(max_iters): 
        Ridder_matrix[i,:] = central_differentation(x, h * inv_d ** i, function)
        
        
    last = Ridder_matrix[0].copy()
    previous_error = np.inf
    for i in range(max_iters - 1): 
        j = i+1
        Ridder_matrix[i] = (d**(2*j) * Ridder_matrix[i+1] - Ridder_matrix[i])/ (d**(2*j) - 1)
        error = np.abs(last - Ridder_matrix[i])
        # print(error)
        last = Ridder_matrix[0].copy()
        if np.any(error < e): 
            print(f'Error is reached at iteration {i+1}')
            break
        
        if np.any(error > previous_error): 
            print(f'Error is growing after iteration {i+1}')
            return Ridder_matrix[0]
        
        previous_error = error
    return Ridder_matrix[0]

def test(x):
    return x**3

result_diff = Ridder(np.array([1]), 0.1, 2, 1e-10, 100, n)[0] 
result_analytical_diff = analytical_diff(1)

print(result_diff)
print(result_analytical_diff)


