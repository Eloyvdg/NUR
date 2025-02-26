import numpy as np

class solver():   
    ''' Class to perform interpolation via LU decomposition. 
    Attributes: 
    x: list of real x values
    y: list of real y values '''
    
    def __init__(self, x, y): 
        self.x = x
        self.y = y
        self.matrix = solver._crout(self.x)
    
    @staticmethod
    def _vandermonde(x): 
        ''' Calculate the Vandermonde matrix given a list of x values'''
        size = len(x)
        matrix = np.zeros((size, size))

        for i in range(matrix.shape[1]):
            for j in range(matrix.shape[1]):
                matrix[i,j] = x[i]**j
        return matrix
    
    @staticmethod
    def _crout(x): 
        ''' Calculate the LU decomposition using Crout theorem. 
        Input: 
        x: list of real x values '''
        
        decomposition = solver._vandermonde(x).astype(float).copy()
        
        for j in range(decomposition.shape[1]): # Crout theorem
            for i in range(j+1): 
                for k in range(i):
                    decomposition[i,j] -= decomposition[i,k] * decomposition[k,j]     
                           
            for i in range(j+1, decomposition.shape[1]):
                for k in range(j): 
                    decomposition[i,j] -= decomposition[i,k] * decomposition[k,j]             
                decomposition[i,j] /= decomposition[j,j]
                
        return decomposition
    
    def solve(self):
        ''' Solve to find c in Vc = y using the LU decomposition, and forward and backward substitution'''
        c = self.y.copy()
        for i in range(len(c)): # Forward substitution
            for j in range(i): 
                c[i] -= self.matrix[i,j] * c[j]
                
        for i in range(len(c)-1, -1, -1): # Backward substitution
            for j in range(i+1, len(c)):
                c[i] -= self.matrix[i,j] * c[j]
            
            c[i] /= self.matrix[i,i]
        return c
    
    def iteration(self, iterations): 
        ''' Solve find c in Vc = y using the LU decomposition, and forward and backward substitution.
        In this definition, it is an iterated process to try improve the result.
        Input: 
        iterations: the amount of iterations wanted'''
        start = solver(self.x, self.y)
        b = start.solve()        
        original = solver._vandermonde(self.x)
        for i in range(iterations): 
            error = original @ b - self.y
            improved = solver(self.x, error)
            b -= improved.solve()
            
        return b
