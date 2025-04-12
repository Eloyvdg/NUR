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
    func_integrate = lambda x: 4* np.pi * x**2* n(x, A, Nsat, a, b, c)
    result, result_error = romberg(np.array([xmin]), np.array([xmax]), 10, func_integrate)
    # result, result_error = romberg(xmin, xmax, 10, N(A, Nsat, a, b, c))
    return 1/result[0], func_integrate

def likelihood(A, Nsat, point, edges, bins): 

    # A_new, func_integrate = normalize(A, Nsat, *point)
    A_new, error = romberg(np.array([xmin]), np.array([xmax]), 10, lambda x: N(x, A, Nsat, *point))
    A_new = 1/A_new[0]
    # print(A_new)
    # A_new = romberg(xmin, xmax, 10, N)
    # likelihood = 0
    # variances = np.zeros(len(edges)-1)
    # print(edges)
    edges = np.exp(np.linspace(np.log(1e-4), np.log(5), 100+1))
    variances, error = romberg(edges[:-1], edges[1:], 10, lambda x: A_new * N(x, A, Nsat, *point))
    variances = variances[0,:]
    # for i in range(len(edges)-1):
    #     # variance, error = romberg(edges[i], edges[i+1], 10, lambda x: A_new * func_integrate(x))
    #     variance, error = romberg(edges[i], edges[i+1], 10, lambda x: A_new * N(x, A, Nsat, *point))
    #     variances[i] = variance[0]

        # bins_var = bins.copy()
        # bins_var[bins_var == 0] = 1e-5
        # likelihood += (bins[i] - variance[0])**2/(variance[0])   
    
    likelihood = np.sum((bins - variances)**2/variances)
    return likelihood

def likelihood_poisson(A, Nsat, point, edges, bins, radii): 
    
    A_new, func_integrate = normalize(A, Nsat, *point)
    likelihood = -np.sum(np.log(A_new * func_integrate(radii)), axis = 0)
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

# def g_test(): 
#     2 * np.sum(observed * np.log(observed/expected))

#Call this function as: 
Nsat_calc = []
radii = []
Nhalo = []
filenames = ['satgals_m11.txt', 'satgals_m12.txt', 'satgals_m13.txt', 'satgals_m14.txt', 'satgals_m15.txt']
for i in filenames: 
    radius, nhalo = readfile(i)
    radii += [radius]
    Nhalo += [nhalo]
    Nsat_calc += [len(radius)/nhalo]


# Plot of binned data with the best fit (question 1b.4 and 1c)
# As always, feel free to replace by your own plotting routines if you want
xmin, xmax = 1e-4, 5. # replace by your choices
n_bins = 100 # replace by your binning
edges = np.exp(np.linspace(np.log(xmin), np.log(xmax), n_bins+1))
ylim_min = [1e-8, 1e-12, 1e-5, 1e-3, 1e-2]
ylim_max = [10, 1e3    , 1e3   , 1e3 , 1e4]
fig1b, ax = plt.subplots(3,2,figsize=(6.4,8.0))
for i in range(5):
    Nsat = Nsat_calc[i] # replace by actual appropriate number for mass bin i
    x_radii = radii[i] # replace by actual data for mass bin i
    
    # simplex = np.array([[2.4, 0.2, 1.6], 
    #                     [2.3, 0.1, 1.2], 
    #                     [2.0, 0.5, 1.4], 
    #                     [2.6, 1.1, 1.8]])
    
    # simplex = np.array([[2.4, 0.2, 1.6], 
    #                     [1.0, 0.7, 1.9], 
    #                     [2.8, 1.8, 0.6], 
    #                     [0.5, 1.3, 2.9]])
    
    # simplex = np.array([[2.4, 0.2, 1.6], 
    #                     [2.9, 0.7, 2.1], 
    #                     [3.6, 1.4, 2.8], 
    #                     [3.8, 1.6, 3.0]])
    
    simplex = np.array([[2.4,  0.2,  1.6 ],
                        [2.45, 0.2,  1.6 ],
                        [2.4,  0.25, 1.6 ],
                        [2.4,  0.2,  1.65]])
    
    # simplex = np.array([[2.4, 0.2, 1.6], 
    #                         [1.0, 0.7, 1.9], 
    #                         [2.8, 1.8, 0.6], 
    #                         [0.5, 1.3, 2.9]])
    
    
    # simplex = np.array([[1.31420016, 1.11502032, 3.13383544], 
    #                         [1.31420016, 1.11502032, 3.13383544], 
    #                         [1.31420016, 1.11502032, 3.13383544], 
    #                         [1.31420016, 1.11502032, 3.13383544]]) * np.random.uniform(0.95, 1.05, size=(4,3))
    
    # simplex = np.array([[1.5, 2.5, 0.2],
    #                         [1, 0.5, 1.8],
    #                         [1.8, 1.6, 1],
    #                         [0.8, 2.3, 5.5]]) * 10
    
    print('AHHHH')
    
    binned_data=np.histogram(x_radii,bins=edges)[0]/Nhalo[i]#/np.diff(edges)
    
    # print(edges)
    likelihood_calc = lambda point: likelihood(1, Nsat, point, edges, binned_data)
    result_downhill = downhill_simplex(simplex, likelihood_calc, 1000, 1e-10)
    print(result_downhill)
    best_a, best_b, best_c = result_downhill #[1.53874748,   0.47773451,   1.19624515]
    
    func_integrate = lambda x: 4* np.pi * x**2 * n(x, 1, Nsat, best_a, best_b, best_c)
    result_A, result_error_A = romberg(np.array([xmin]), np.array([xmax]), 10, func_integrate)
        
    x_array = np.logspace(-4, np.log10(5), 100)
    x_centers = 0.5 * (edges[1:] + edges[:-1])
    
    # Ntilda = lambda x: 4*np.pi * x**2 * n(x, 1/result_A[0], Nsat, best_a, best_b, best_c) * Nsat * np.diff(edges)
    Ntilda = N(x_centers, 1/result_A[0], Nsat, best_a, best_b, best_c) * np.diff(edges)

    row=i//2
    col=i%2
    ax[row,col].step(edges[:-1], binned_data, where='post', label='binned data')
    ax[row,col].step(x_centers, Ntilda, where='post', label='best-fit profile')
    # ax[row,col].step(x_array, N(x_array, 1/result_A[0], Nsat, best_a, best_b, best_c), where='post', label='best-fit profile')


    # ax[row,col].step(edges[:-1], Ntilda, where='post', label='best-fit profile')
    ax[row,col].set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
    ax[row,col].set_ylim(ylim_min[i], ylim_max[i])
ax[2,1].set_visible(False)
plt.tight_layout()
handles,labels=ax[2,0].get_legend_handles_labels()
plt.figlegend(handles, labels, loc=(0.65,0.15))
plt.savefig('my_solution_1b.png', dpi=600)

# Plot 1c (same code as above)
# fig1c, ax = plt.subplots(3,2,figsize=(6.4,8.0))
# for i in range(5):
#     print('AHHH')
#     Nsat = Nsat_calc[i] # replace by actual appropriate number for mass bin i
#     x_radii = radii[i] # replace by actual data for mass bin i
    
#     binned_data=np.histogram(x_radii,bins=edges)[0]/Nhalo[i]/np.diff(edges)
    
#     # simplex = np.array([[2.4, 0.2, 1.6], 
#     #                     [2.3, 0.1, 1.2], 
#     #                     [2.0, 0.5, 1.4], 
#     #                     [2.6, 1.1, 1.8]])
#     simplex = np.array([[2.4, 0.2, 1.6], 
#                             [1.0, 0.7, 1.9], 
#                             [2.8, 1.8, 0.6], 
#                             [0.5, 1.3, 2.9]])
    
#     likelihood_calc = lambda point: likelihood_poisson(1, Nsat, point, edges, binned_data, x_radii)
#     result_downhill = downhill_simplex(simplex, likelihood_calc, 1000, 1e-10)
#     best_a, best_b, best_c = result_downhill
    
#     print(result_downhill)

#     func_integrate = lambda x: 4* np.pi *Nsat * x**2 * n(x, 1, Nsat, best_a, best_b, best_c)
#     result_A, result_error_A = romberg(np.array([xmin]), np.array([xmax]), 10, func_integrate)
    
#     x_array = np.logspace(-4, np.log10(5), 100)
#     row=i//2
#     col=i%2
#     ax[row,col].step(edges[:-1], binned_data, where='post', label='binned data')
#     ax[row,col].plot(x_array, N(x_array, 1/result_A[0], Nsat, best_a, best_b, best_c)*Nsat, label='best-fit profile')

#     # ax[row,col].step(edges[:-1], Ntilda, where='post', label='best-fit profile')
#     ax[row,col].set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
#     ax[row,col].set_ylim(ylim_min[i], ylim_max[i])


# ax[2,1].set_visible(False)
# plt.tight_layout()
# handles,labels=ax[2,0].get_legend_handles_labels()
# plt.figlegend(handles, labels, loc=(0.65,0.15))
# plt.savefig('my_solution_1c.png', dpi=600)
