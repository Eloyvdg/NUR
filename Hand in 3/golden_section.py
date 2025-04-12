import numpy as np

golden_ratio = (1 + np.sqrt(5))/2
w = 1/(1+golden_ratio)


def bracketing(xmin, xmax, max_iter, func): 
    '''Function to find an initial bracket
        input: 
        xmin: initial minimum of bracket
        xmax: initial maximum of bracket
        max_iter: maximal amount of iterations before stopping
        func: function to find the bracket on
        '''
    a = xmin
    b = xmax
    
    if func(b) > func(a): # switch a and b if f(b) > f(a)
        a, b = b, a
        
    c = b + (b - a) * w # New point c

    for i in range(max_iter): 
        if func(c) > func(b): # Bracket found
            return [a, b, c]
        
        num = (b-a)**2 * (func(b) - func(c)) - (b-c)**2 * (func(b) - func(a))
        den = (b-a) * (func(b) - func(c)) - (b - c) * (func(b) - func(a))
        
        d = b - 0.5 * num / den # Propose new point d based on parabola
        
        if b < d < c: 
            if func(d) < func(c): 
                return [b, d, c] 
            elif func(d) > func(b): 
                return [a, b, d]
            else: 
                d = c + (c - b) * w
        elif d > c: 
            if np.abs(d - b) > 100 * np.abs(c - b):
                d = c + (c - b) * w # too far away -> propose new d

        else: 
            print('Initial bracket not able to find minumum')
            return None
         
        a, b, c = b, c, d
        return None
        
    
    a, b, c = b, c, d         
    
    
def golden_section(xmin, xmax, target, max_iter, func): 
    bracket = bracketing(xmin, xmax, max_iter, func)
    a, b, c = bracket[0], bracket[1], bracket[2]

    for i in range(max_iter): 
        if np.abs(c - b) < np.abs(b - a):  # Find the smallest interval
            d = b + (a - b) * w
        else: 
            d = b + (c - b) * w

        if np.abs(c - a) < target: # CHeck if target accuracy is reached
            if func(d) < func(b): 
                return d
            elif func(d) > func(b): 
                return b
        
        if func(d) < func(b): # Tighten the bracket
            if b < d < c:
                a, b = b, d
            elif a < d < b: 
                c, b = b, d
        else: 
            if b < d < c:
                c = d
            elif a < d < b: 
                a = d
