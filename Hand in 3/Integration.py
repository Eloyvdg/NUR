import numpy as np

        
def romberg(a, b, m, function): 
    """
    Romberg integration algorithm.
    m: Order of integration
    """
    # print(a)
    r = np.zeros((m, len(a)))
    
    # print(r)
    h = b - a
    # print(h)
    # Calculate r with decreasing stepsize h
    r[0,:] = 0.5 * h * (function(a) + function(b))
    # print(r[0,:], ' AHH')
    
    # print((r))
    N_p = 1
    for i in range(1,m): 
        r[i,:] = 0 
        delta = h
        h = 0.5*h
        x = a + h
        # print(x)
        # print(a, h)
        for j in range(N_p): 
            r[i,:] += function(x)
            x += delta
            
        r[i,:] = 0.5 * (r[i-1,:] + delta * r[i,:])
        
        N_p *= 2    
    
    # Combine the results to increase accuracy
    N_p = 1
    for i in range(1,m): 
        N_p *= 4
        for j in range(0, m-i): 
            r[j,:] = (N_p * r[j+1,:] - r[j,:]) / (N_p - 1)
                
    error = np.abs(r[0,:] - r[1,:])
    return r, error

def func(x):
    return x**2
