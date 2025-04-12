import numpy as np

        
def romberg(a, b, m, function): 
    """
    Romberg integration algorithm.
    input: 
    a: array of lower limits
    b: array of upper limits
    m: Order of integration
    function: function to integrate
    """
    r = np.zeros((m, len(a)))
    
    h = b - a
    # Calculate r with decreasing stepsize h for all limits in array
    r[0,:] = 0.5 * h * (function(a) + function(b))
    
    N_p = 1
    for i in range(1,m): 
        r[i,:] = 0 
        delta = h
        h = 0.5*h
        x = a + h
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
