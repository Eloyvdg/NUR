import numpy as np
import matplotlib.pyplot as plt

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

# def false_position(): 
    
    
    
def func_root(x): 
    return 0.5*x - 3

root = bisection(0, 10, func_root)

plt.plot(np.arange(0,10), func_root(np.arange(0,10)))
plt.axhline(0)
plt.axvline(root)

