import numpy as np
import matplotlib.pyplot as plt
from Integration import romberg
from Minimizing import golden_section
from Minimizing import downhill_simplex

xmin = 1e-4
xmax = 5

def n(x,A,Nsat,a,b,c):
    return A*((x/b)**(a-3))*np.exp(-(x/b)**c)

def N(x,A,Nsat,a,b,c):
    return Nsat * A*((x/b)**(a-3))*np.exp(-(x/b)**c) * np.pi * 4 * x**2

def normalize(A, Nsat, a, b, c, xmin=xmin, xmax = xmax): 
    '''Definitionto normalize n(x) * 4pi * x^2'''
    func_integrate = lambda x: 4* np.pi * x**2* n(x, A, Nsat, a, b, c)
    result, result_error = romberg(np.array([xmin]), np.array([xmax]), 10, func_integrate)
    return 1/result[0], func_integrate
    

def likelihood_poisson(A, Nsat, point, edges, bins, radii): 
    '''Function to calculate the Poisson likelihoof (without constant factors)'''
    A_new, func_integrate = normalize(A, Nsat, *point)
    result = A_new * func_integrate(radii)
    result = result[result != 0]
    likelihood = -np.sum(np.log(result), axis = 0)
    return likelihood
        
def readfile(filename):
    f = open(filename, 'r')
    data = f.readlines()[3:] #Skip first 3 lines 
    nhalo = int(data[0]) #number of halos
    radius = []
    
    for line in data[1:]:
        if line[:-1]!='#':
            radius.append(float(line.split()[0]))
    
    radius = np.array(radius, dtype=float)    
    f.close()
    return radius, nhalo #Return the virial radius for all the satellites in the file, and the number of halos


#Call this function as: 

filenames = ['satgals_m11.txt', 'satgals_m12.txt', 'satgals_m13.txt', 'satgals_m14.txt', 'satgals_m15.txt']

    
xmin, xmax = 1e-4, 5.
n_bins = 100 
edges = np.exp(np.linspace(np.log(xmin), np.log(xmax), n_bins+1))
ylim_min = [1e-8, 1e-12, 1e-5, 1e-3, 1e-2]
ylim_max = [10, 1e3    , 1e3   , 1e3 , 1e4]
fig1b, ax = plt.subplots(3,2,figsize=(6.4,8.0))


# Plot 1c (same code as above)
fig1c, ax = plt.subplots(3,2,figsize=(6.4,8.0))
for i in range(5):

    
    x_radii, Nhalo = readfile(filenames[i])
    Nsat = len(x_radii)/Nhalo # Calculate Nsat
    
    k = n_bins - 4
    
    binned_data=np.histogram(x_radii,bins=edges)[0]/Nhalo/np.diff(edges)
    simplex = np.array([[2.4, 0.2, 1.6], 
                            [1.0, 0.7, 1.9], 
                            [2.8, 1.8, 0.6], 
                            [0.5, 1.3, 2.9]])
    
    likelihood_calc = lambda point: likelihood_poisson(1, Nsat, point, edges, binned_data, x_radii)
    result_downhill, likelihood_result = downhill_simplex(simplex, likelihood_calc, 1000, 1e-10)
    best_a, best_b, best_c = result_downhill
    
    print(f''' The best fit paramaters for {filenames[i]} are:
          a = {result_downhill[0]}
          b = {result_downhill[1]}
          c = {result_downhill[2]}''')
          
    print(f''' 
          The corresponding log likelihood value is: {likelihood_result * Nhalo} \n
          ##################################################################\n''')
          
    func_integrate = lambda x: 4* np.pi *Nsat * x**2 * n(x, 1, Nsat, best_a, best_b, best_c) 
    result_A, result_error_A = romberg(np.array([xmin]), np.array([xmax]), 10, func_integrate) # Normalize again with bets found parameters
    
    x_array = np.logspace(-4, np.log10(5), 100)
    
    row=i//2
    col=i%2
    
    # Plot the results
    ax[row,col].step(edges[:-1], binned_data, where='post', label='binned data')
    ax[row,col].plot(x_array, N(x_array, 1/result_A[0], Nsat, best_a, best_b, best_c)*Nsat, label='best-fit profile')

    ax[row,col].set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
    ax[row,col].set_ylim(ylim_min[i], ylim_max[i])


ax[2,1].set_visible(False)
plt.tight_layout()
handles,labels=ax[2,0].get_legend_handles_labels()
plt.figlegend(handles, labels, loc=(0.65,0.15))
plt.savefig('my_solution_1c.png', dpi=600)
