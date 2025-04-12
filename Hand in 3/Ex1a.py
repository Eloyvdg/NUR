import numpy as np
from golden_section import golden_section

A= 256/(5*np.pi**(3/2))
Nsat=100
a=2.4
b=0.25
c=1.6
xmax = 5

golden_ratio = (1 + np.sqrt(5))/2
w = 1/(1+golden_ratio)

def N(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    return (A*((x/b)**(a-3))*np.exp(-(x/b)**c)) * x**2 * 4 * np.pi

def N_negative(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    return -(A*((x/b)**(a-3))*np.exp(-(x/b)**c)) * x**2

minimum = golden_section(1e-4, 1, 1e-5, 100, N_negative) # Calculate the minimum of the negative function
n_max = N(minimum) #Calculate the corresponding maximum value of N(x)

print(f'The maxmimum of N(x) is found at x = {minimum}')
print(f'The corresponding N(x) is {n_max}')
