import numpy as np

def factorial(n): 
    if n == 0 or n == 1:
        return 1
    x = 1
    for i in range(2, n+1):
        x *= i    
    return x

def poisson(lamb, k):
    return (lamb**k * np.exp(-lamb))/factorial(k)


array_lambda = np.arange(0, 102)
array_k = np.arange(0, 102)
