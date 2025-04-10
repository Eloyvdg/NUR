import numpy as np
import matplotlib.pyplot as plt

A= 256/(5*np.pi**(3/2))
Nsat=100
a=2.4
b=0.25
c=1.6
xmax = 5

golden_ratio = (1 + np.sqrt(5))/2
w = 1/(1+golden_ratio)

def func_integrate(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    """
    Function to integrate
    """
    return (A*((x/b)**(a-3))*np.exp(-(x/b)**c)) * x**2 * 4 * np.pi


def func_integrate2(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    """
    Function to integrate
    """
    return -(A*((x/b)**(a-3))*np.exp(-(x/b)**c)) * x**2 * 4 * np.pi

def func(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    """
    Function to integrate
    """
    return A* Nsat * ( (x/b)**(a-3) ) *np.exp(-(x/b)**c)


def function_test(x): 
    return x * np.sin(x)

def bracketing(xmin, xmax, max_iter, func): 
    a = xmin
    b = xmax
    
    if func(b) > func(a): 
        a, b = b, a
        
    c = b + (b- a) * w

    for i in range(max_iter): 
        if func(c) > func(b):
            return [a, b, c]
        
        num = (b-a)**2 * (func(b) - func(c)) - (b-c)**2 * (func(b) - func(a))
        den = (b-a) * (func(b) - func(c)) - (b - c) * (func(b) - func(a))
        
        d = b - 0.5 * num / den
        
        if d > b and d < c: 
            if func(d) < func(c): 
                return [b, d, c] 
            elif func(d) > func(b): 
                return [a, b, d]
            else: 
                d = c + (c - b) * w
        elif d > c: 
            if np.abs(d - b) > 100 * np.abs(c - b):
                d = c + (c - b) * w
            else: 
                a, b, c = b, c, d
        else: 
            print('Initial bracket not able to find minumum')
            return None
         
        return None
        
    
    a, b, c = b, c, d         
    
    
def golden_section(xmin, xmax, target, max_iter, func): 
    # a, b, c = bracket[0], bracket[1], bracket[2]
    bracket = bracketing(xmin, xmax, max_iter, func)
    a, b, c = bracket[0], bracket[1], bracket[2]

    for i in range(max_iter): 
        if np.abs(c - b) < np.abs(b - a): 
            d = b + (a - b) * w
        else: 
            d = b + (c - b) * w

        if np.abs(c - a) < target: 
            if func(d) < func(b): 
                return d
            elif func(d) > func(b): 
                return b
        
        if func(d) < func(b): 
            if b < d < c:
                a, b = b, d
            elif a < d < b: 
                c, b = b, d
        else: 
            if b < d < c:
                c = d
            elif a < d < b: 
                a = d

xmin = 1e-4
xmax = 5

x = 10**np.linspace(np.log10(xmin), np.log10(xmax), 200) #replace!

plt.plot(x, func_integrate(x))
bracket = bracketing(1e-4, 1, 1000,func_integrate2)

minimum = golden_section(1e-4, 1, 1e-5, 100,func_integrate2)
# maximum = minimum - func_integrate(x)

# plt.plot(bracket, func_integrate(bracket), '.')
plt.plot(minimum, func_integrate(minimum), '.', label = 'Ja, we hebben hem!')
plt.yscale('log')
plt.ylim(10e-3, 10)
plt.xscale('log')
plt.legend()
plt.show()


    
