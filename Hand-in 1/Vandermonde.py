import numpy as np

class solver(): 
    
    def __init__(self, x, y): 
        self.x = x
        self.y = y
    
    @staticmethod
    def _vandermonde(x): 
        size = len(x)
        matrix = np.zeros((size, size))

        for i in range(matrix.shape[1]):
            for j in range(matrix.shape[1]):
                matrix[i,j] = x[i]**j
        # self.matrix = matrix
        return matrix
    
    @staticmethod
    def _crout(x): 
        decomposition = solver._vandermonde(x).astype(float).copy()
        #decomposition= self.matrix.astype(float).copy()
        
        for j in range(decomposition.shape[1]): 
            for i in range(j+1): 
                    for k in range(i):
                        decomposition[i,j] -= decomposition[i,k] * decomposition[k,j]     
                           
            for i in range(j+1, decomposition.shape[1]):
                for k in range(j): 
                    decomposition[i,j] -= decomposition[i,k] * decomposition[k,j]             
                decomposition[i,j] /= decomposition[j,j]
                
        return decomposition
    
    def solve(self):
        matrix = solver._crout(x)
        b = self.y.copy()
        for i in range(len(b)): 
            for j in range(i): 
                b[i] -= matrix[i,j] * b[j]
        
        for i in range(len(b), 0):
            for j in range(i+1, len(b)):
                b[i] -= matrix[i,j] * b[j]
        return b
    

if __name__ == '__main__':
    data=np.genfromtxt(("Vandermonde.txt"),comments='#',dtype=np.float64)
    x=data[:,0]
    y=data[:,1]
    start = solver(x, y)
    
    LU = start._crout(x)
    b = start.solve()
