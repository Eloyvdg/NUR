import numpy as np

class Root_finding(): 
    
    def __init__(self, xmin, xmax):
        self.xmin = xmin
        self.xmax = xmax
        
    def bisection(self, function): 
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
        a, b = self.xmin, self.xmax
        for i in range(n):
            b, a = b - (b -a) / (function(b) - function(a)) * function(b), b
            
        return a, np.abs(b-a)/b
                
    def false_position(self, n, target, target_rel, function): 
        a, b = self.xmin, self.xmax
        best = 1
        for i in range(n): 
            new = b - (b -a) / (function(b) - function(a)) * function(b)
            if function(new)*function(a) < 0: 
                b  = new
            elif function(new) * function(b) < 0: 
                a = new   
                
            error = np.abs(best - new)
            rel_error  = error/best
            best = new

            if (error < target) and (rel_error < target_rel): 
                break
            
        return best, rel_error
          
    def newton_raphson(self, n, target, target_rel, function, diff_function): 
        best = (self.xmax + self.xmin) * 0.5
        for i in range(n): 
            new = best - function(best) / diff_function(best)

            error = np.abs(best - new)

            rel_error  = error/best

            if (error < target) and (rel_error < target_rel): 
                return best, rel_error
            
            best = new
     
        return best, rel_error