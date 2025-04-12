import numpy as np
import matplotlib.pyplot as plt
from Integration import romberg
from downhill import downhill_simplex
import scipy.special as sc


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

def likelihood(A, Nsat, point, edges, bins): 
    '''Definition to calculate the likelihood of a Chi^ (without constant factors)'''
    A_new, error = romberg(np.array([xmin]), np.array([xmax]), 10, lambda x: N(x, A, Nsat, *point))
    A_new = 1/A_new[0] # Find the normalization factor A of the function

    edges = np.exp(np.linspace(np.log(1e-4), np.log(5), 100+1))
    variances, error = romberg(edges[:-1], edges[1:], 10, lambda x: A_new * N(x, A, Nsat, *point)) # Calculate variance and mean for every bin
    variances = variances[0,:]
    
    likelihood = np.sum((bins - variances)**2/variances)
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

def g_test(observed, expected): 
    idx_nonzero = [observed != 0][0]
    return 2 * np.sum(observed[idx_nonzero] * np.log(observed[idx_nonzero]/expected[idx_nonzero]))

def significance(x, k): 
    P = sc.gammainc(k*0.5, x*0.5) / sc.gamma(k/2)
    print(P, 'P')
    return 1 - P


filenames = ['satgals_m11.txt', 'satgals_m12.txt', 'satgals_m13.txt', 'satgals_m14.txt', 'satgals_m15.txt']

xmin, xmax = 1e-4, 5.
n_bins = 100
edges = np.exp(np.linspace(np.log(xmin), np.log(xmax), n_bins+1))
ylim_min = [1e-8, 1e-12, 1e-5, 1e-3, 1e-2] # Define the limits for the plots
ylim_max = [10, 1e3    , 1e3   , 1e3 , 1e4]
fig1b, ax = plt.subplots(3,2,figsize=(6.4,8.0))

for i in range(5):
    x_radii, Nhalo = readfile(filenames[i])
    Nsat = len(x_radii)/Nhalo # Calculate Nsat
    
    simplex = np.array([[2.4,  0.2,  1.6 ],
                        [2.45, 0.2,  1.6 ],
                        [2.4,  0.25, 1.6 ],
                        [2.4,  0.2,  1.65]]) # Define the initial simplex
    
    binned_data=np.histogram(x_radii,bins=edges)[0]/Nhalo
    
    likelihood_calc = lambda point: likelihood(1, Nsat, point, edges, binned_data)
    result_downhill, likelihood_result = downhill_simplex(simplex, likelihood_calc, 1500, 1e-10) # Perform downhill simplex on likelihood
    k = len(x_radii) - 4

    best_a, best_b, best_c = result_downhill
    print(f''' The best fit paramaters for {filenames[i]} are:
          a = {result_downhill[0]}
          b = {result_downhill[1]}
          c = {result_downhill[2]}''')
          
    print(f''' 
          The corresponding Chi^2 value is: {likelihood_result * Nhalo} \n
          The corresponding Chi^2/k value is {likelihood_result * Nhalo/k}
          ##################################################################\n''')
    
    func_integrate = lambda x: 4* np.pi * x**2 * n(x, 1, Nsat, best_a, best_b, best_c)
    result_A, result_error_A = romberg(np.array([xmin]), np.array([xmax]), 10, func_integrate) # Normalize again with bets found parameters
        
    x_array = np.logspace(-4, np.log10(5), 100)
    x_centers = 0.5 * (edges[1:] + edges[:-1])
    
    Ntilda, Ntilda_error = romberg(edges[:-1], edges[1:], 10, lambda x: 1/result_A[0] * N(x, 1, Nsat, best_a, best_b, best_c)) # Integrate over bins
    Ntilda= Ntilda[0,:]
    
    result_G_test = g_test(binned_data*Nhalo, Ntilda*Nhalo)

    row=i//2
    col=i%2
    ax[row,col].step(edges[:-1], binned_data, where='post', label='binned data')
    ax[row,col].step(x_centers, Ntilda, where='post', label='best-fit profile')

    ax[row,col].set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
    ax[row,col].set_ylim(ylim_min[i], ylim_max[i])
ax[2,1].set_visible(False)
plt.tight_layout()
handles,labels=ax[2,0].get_legend_handles_labels()
plt.figlegend(handles, labels, loc=(0.65,0.15))
plt.savefig('my_solution_1b.png', dpi=600)
