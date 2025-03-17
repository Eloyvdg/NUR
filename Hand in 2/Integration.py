import numpy as np

class Integrator():
    """
    Class to integrate a function numerically
    """
    def __init__(self, function, a, b): 
        """
        Define the function and boundaries
        function: function to integrate
        a: lower boundary
        b: upper boundary
        """
        self.function = function  
        self.a = a
        self.b = b
        
    def romberg(self, m): 
        """
        Romberg integration algorithm.
        m: Order of integration
        """
        r = np.zeros(m)
        h = self.b - self.a
        # Calculate r with decreasing stepsize h
        r[0] = 0.5 * h * (self.function(self.a) + self.function(self.b))
        N_p = 1
        for i in range(1,m): 
            r[i] = 0 
            delta = h
            h = 0.5*h
            x = self.a + h
            for j in range(N_p): 
                r[i] += self.function(x)
                x += delta
                
            r[i] = 0.5 * (r[i-1] + delta * r[i])
            
            N_p *= 2    
        
        # Combine the results to increase accuracy
        N_p = 1
        for i in range(1,m): 
            N_p *= 4
            for j in range(0, m-i): 
                r[j] = (N_p * r[j+1] - r[j]) / (N_p - 1)
                    
        error = np.abs(r[0] - r[1])
        return r, error
