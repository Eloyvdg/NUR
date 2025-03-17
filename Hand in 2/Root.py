import numpy as np
import time

class Root_finding(): 
    ''' 
    Class for different root finders
    ''' 
    def __init__(self, xmin, xmax):
        """Define the bracket"""
        self.xmin = xmin
        self.xmax = xmax
        
    def bisection(self, function): 
        """
        Bisection method to find the root
        function: Function to find the root of
        """
        a, b = self.xmin, self.xmax
        
        iterations = int(np.log2( (b-a) / 0.0001))
        for i in range(iterations):
            c = (self.xmin + self.xmax)*0.5
            if function(a)*function(c) < 0: 
                self.xmax = c
            else:
                self.xmin = c
        return c

    def secant(self, n, function):
        """
        Secant method to find the root
        n: maximum number of iterations
        function: function to find the root of 
        """
        a, b = self.xmin, self.xmax
        for i in range(n):
            b, a = b - (b -a) / (function(b) - function(a)) * function(b), b
            
        return a, np.abs(b-a)/b
                
    def false_position(self, n, target, target_rel, function): 
        """
        False position method to find the root
        n: maximum number of iterations
        target: target absolute error
        target_rel: target relative error
        function: function to find the root of 
        """
        a, b = self.xmin, self.xmax
        best = 1
        iterations = 0
        start = time.time()
        for i in range(n):
            # Equation to find new best point
            new = b - (b -a) / (function(b) - function(a)) * function(b)
            # Calculate if a or b should be replaced with the new best point
            if function(new)*function(a) < 0: 
                b  = new
            elif function(new) * function(b) < 0: 
                a = new   
                
            error = np.abs(best - new)
            rel_error  = error/best
            best = new
            
            # Stop if the target errors are reached
            if (error < target) and (rel_error < target_rel): 
                end = time.time()
                print(f'It took {iterations} iterations and {end-start} seconds to find the root.')
                return best, rel_error
            
            iterations += 1
        return best, rel_error
          
    def newton_raphson(self, n, target, target_rel, function, diff_function): 
        """
        Newton-Raphson method to find the root
        n: maximum number of iterations
        target: target absolute error
        target_rel: target relative error
        function: function to find the root of 
        diff_function: Derivative of function
        """
        best = (self.xmax + self.xmin) * 0.5
        start = time.time()
        iterations = 0
        for i in range(n): 
            # Calculate the new best point
            new = best - function(best) / diff_function(best)

            error = np.abs(best - new)

            rel_error  = error/best

            # Stop if the target errors are reached
            if (error < target) and (rel_error < target_rel): 
                end = time.time()
                print(f'It took {iterations} iterations and {end-start} seconds to find the root.')
                return best, rel_error
            
            best = new
            iterations += 1
     
        return best, rel_error
