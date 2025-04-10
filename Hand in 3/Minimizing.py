import numpy as np
import matplotlib.pyplot as plt
from quicksort import quick_sort 

golden_ratio = (1 + np.sqrt(5))/2
w = 1/(1+golden_ratio)

def g_plotting(x, y):
    return -np.exp(x*x - y*y)

def g(point): 
    x, y = point
    return -np.exp(-x*x - y*y)
    
def downhill_simplex(simplex_start, func, max_iter, target): 
    
    simplex = simplex_start.copy()
    for i in range(max_iter): 
        y = []
        for j in range(simplex.shape[0]): 
            y += [func(simplex[j,:])]
            
        y = np.array(y)  
        # print(y)
        # y = func(simplex_start, simplex)
        y_sorted = quick_sort(y)
        sort_dict = {value: key for key, value in enumerate(y)}

        sorted_idx = [sort_dict[key] for key in y_sorted]
        
        simplex = simplex[sorted_idx,:]
        # print(simplex)
        print(simplex)
        # print(simplex)
        centroid = np.sum(simplex[:-1,:], axis = 0)/(simplex.shape[1])    
        
        if np.abs(y_sorted[-1] - y_sorted[0]) / np.abs(0.5 * (y_sorted[-1] + y_sorted[0])) < target:
            print('Target reached')
            return simplex[0,:]
            
        
        x_try = np.array([2 * centroid - simplex[-1,:]])[0]

        
        func_xtry = func(x_try)
        
        if y_sorted[0] <= func_xtry < y_sorted[-1]: # I
            print('I')
            simplex[-1,:] = x_try
            
        elif func_xtry < y_sorted[0]:   #II
            print('II')
            x_exp = 2 * x_try - centroid
            if func(x_exp) < func_xtry: 
                simplex[-1,:] = x_exp
            else: 
                simplex[-1,:] = x_try
                
        else: # func_xtry  >= y_sorted[-1]:    # III
            x_try = np.array([0.5 * (centroid + simplex[-1,:])])[0]
            func_xtry = func(x_try)
            if func_xtry < y_sorted[-1]: 
                print('III')
                simplex[-1,:] = x_try
            else: 
                print('IV')
                for i in range(1, simplex.shape[1]):
                    simplex[i, :] = 0.5 * (simplex[0,:] + simplex[i,:]) #IV
        # else: 
        #     print('wtf')
            
    return simplex[0,:]


vec = np.array((np.linspace(-2, 2, 10), np.linspace(-2, 2, 10))).T
grid_x, grid_y = np.meshgrid(vec[:,0], vec[:,1])

def bracketing(xmin, xmax, max_iter, func): 
    a = xmin
    b = xmax
    
    if func(b) > func(a): 
        a, b = b, a
        
    c = b + (b- a) * w

    for i in range(max_iter): 
        if func(c) > func(b):
            return [a, b, c]
        
        num = (b-a)**2 * (func(b) - func(c)) - (b-c)**2 * (func(b) - func(a))
        den = (b-a) * (func(b) - func(c)) - (b - c) * (func(b) - func(a))
        
        d = b - 0.5 * num / den
        
        if d > b and d < c: 
            if func(d) < func(c): 
                return [b, d, c] 
            elif func(d) > func(b): 
                return [a, b, d]
            else: 
                d = c + (c - b) * w
        elif d > c: 
            if np.abs(d - b) > 100 * np.abs(c - b):
                d = c + (c - b) * w
            else: 
                a, b, c = b, c, d
        else: 
            print('Initial bracket not able to find minumum')
            return None
         
        return None
        
    
    a, b, c = b, c, d         
    
    
def golden_section(xmin, xmax, target, max_iter, func): 
    # a, b, c = bracket[0], bracket[1], bracket[2]
    bracket = bracketing(xmin, xmax, max_iter, func)
    a, b, c = bracket[0], bracket[1], bracket[2]

    for i in range(max_iter): 
        if np.abs(c - b) < np.abs(b - a): 
            d = b + (a - b) * w
        else: 
            d = b + (c - b) * w

        if np.abs(c - a) < target: 
            if func(d) < func(b): 
                return d
            elif func(d) > func(b): 
                return b
        
        if func(d) < func(b): 
            if b < d < c:
                a, b = b, d
            elif a < d < b: 
                c, b = b, d
        else: 
            if b < d < c:
                c = d
            elif a < d < b: 
                a = d
                
# z = np.zeros((len(vec[:,0]), len(vec[:,1])))
# for i in vec[:,0]:
#     for j in vec[:,1]: 
# #         z[i,j] = g(np.array([[i],[j]]).T)

# z = g_plotting(grid_x, grid_y)
# grid = np.meshgrid((vec[:,0], vec[:,1]))
# plt.contourf(grid_x, grid_y, z, 20)
# plt.show()



if __name__ == '__main__': 
    test_vector = np.array([[9, 2], 
                            [5, 1.5],
                            [3, 10]])
    
    result = downhill_simplex(test_vector, g, 100, 1e-10)
    # print(result)
