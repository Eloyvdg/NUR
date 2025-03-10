import numpy as np
import matplotlib.pyplot as plt 

def function2(x): 
    return x**2 * np.sin(x)

def analytical_diff(x): 
    return 2 * x * np.sin(x) + x**2 * np.cos(x)
    
def differentation(x, h, function): 
    if not callable(function): 
        raise ValueError
    
    df = ( function(x+h) - function(x-h) )/ (2*h) 
    return df

def Ridder(x, h, d, e, max_iters, function):
    
    Ridder_matrix = np.zeros([max_iters, len(x)])
    inv_d = 1/d  
    
    for i in range(max_iters): 
        Ridder_matrix[i] = differentation(x, h * inv_d ** i, function)
    
    last = Ridder_matrix[0].copy()
    
    for i in range(max_iters - 1): 
        j = i+1
        Ridder_matrix[i] = (d**(2*j) * Ridder_matrix[i+1] - Ridder_matrix[i])/ (d**(2*j) - 1)
        error = last - Ridder_matrix[0]
        last = Ridder_matrix[0].copy()
        if np.any(error < e): 
            print(f'break at {i}')
            break
        
    return Ridder_matrix[0]
        

array = np.linspace(0, 2*np.pi, 1000)
h = [0.1, 0.01, 0.001]

fig = plt.figure(figsize = (10,6))
for i in h:
    plt.plot(array, differentation(array, i, function2), label = f'Central difference, h = {i}')


plt.plot(array, analytical_diff(array), '--' ,label = 'Analytical derivative')
plt.show()

plt.plot(array, Ridder(array, 0.1, 2, 1e-8, 10, function2))
         