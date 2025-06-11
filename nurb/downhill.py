import numpy as np
from quicksort import quick_sort 

def downhill_simplex(simplex_start, func, max_iter, target): 
    '''Function to find the minimum using a downhill simplex.
        imput: 
        simplex_start: Simplex to start with, (N+1, N) dimensions, where N is the amount of parameters
        func: Function to find the minimum of
        max_iter: Maximal amount of iterations before stopping
        target: Target accuracy'''
    
    simplex = simplex_start.copy()
    y_best = []
    for i in range(max_iter): 
        y = np.zeros(simplex.shape[0])
        for j in range(simplex.shape[0]):# Calculate the values for the sets of parameters
            y[j] = func(simplex[j,:])
            
        y = np.array(y)  
        # print(y)
        sort_dict = {value: key for key, value in enumerate(y)}

        y = quick_sort(y)
        y_best += [y[0]]

        sorted_idx = [sort_dict[key] for key in y] # Sort the array such that f(x0) <= f(x1) ... <= f(xN)
        
        simplex = simplex[sorted_idx,:]
        centroid = np.sum(simplex[:-1,:], axis = 0)/(simplex.shape[1]) # Center of N-1 best points
        
        if np.abs(y[-1] - y[0]) / np.abs(0.5 * (y[-1] + y[0])) < target:
            return simplex[0,:], y[0], y_best # Target reached, return simplex and lowest value
            
        
        x_try = np.array([2 * centroid - simplex[-1,:]])[0] # Reflect the worst point
        
        func_xtry = func(x_try)
        
        if y[0] <= func_xtry < y[-1]:
            # print('I')
            simplex[-1,:] = x_try # Better, not the best, accept the point
            
        elif func_xtry < y[0]: # x_try is the best point, reflect again and check if it is even better
            # print('II')
            x_exp = 2 * x_try - centroid
            if func(x_exp) < func_xtry: 
                simplex[-1,:] = x_exp
            else: 
                simplex[-1,:] = x_try
                
        else:
            x_try = np.array([0.5 * (centroid + simplex[-1,:])])[0]
            func_xtry = func(x_try)
            if func_xtry < y[-1]: 
                # print('III')
                simplex[-1,:] = x_try # Contract worst point if this improves the result
            else: 
                for i in range(1, simplex.shape[0]):
                    # print('IV')
                    simplex[i, :] = 0.5 * (simplex[0,:] + simplex[i,:]) # Contract al worst points towards best point
        
            
    return simplex[0,:], y[0], y_best
    
    
