import numpy as np 

def crout(matrix): 
    
    decomposition= matrix.astype(float).copy()
    
    for j in range(decomposition.shape[1]): 
        for i in range(j+1): 
                for k in range(i):
                    print(j,i,k, 'beta')
                    decomposition[i,j] -= decomposition[i,k] * decomposition[k,j]     
                       
        for i in range(j+1, decomposition.shape[1]):
            for k in range(j): 
                print(j,i,k, 'alpha')
                decomposition[i,j] -= decomposition[i,k] * decomposition[k,j]             
            decomposition[i,j] /= decomposition[j,j]
        
    return decomposition


matrix_test = np.array([[1,2,3],
                        [4,5,6],
                        [7,8,9]])

test = crout(matrix_test)    
print(test)
