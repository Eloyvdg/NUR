import numpy as np
import matplotlib.pyplot as plt
from Integration import romberg
from Minimizing import golden_section
from Minimizing import downhill_simplex

xmin = 1e-4
xmax = 5

def n(x,A,Nsat,a,b,c):
    return A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)

def N(x,A,Nsat,a,b,c):
    return A*((x/b)**(a-3))*np.exp(-(x/b)**c) * np.pi * 4 * x**2

def likelihood(A, Nsat, point, edges, bins): 
    a, b, c = point
    func_integrate = lambda x: 4*np.pi *x**2 * n(x, A, Nsat, a, b, c)
    result, result_error = romberg(xmin, xmax, 10, func_integrate)
    A_new = 1/(result[0])
    likelihood = 0
    for i in range(len(edges)-1):
        variance, error = romberg(edges[i], edges[i+1], 5, lambda x: A_new * func_integrate(x))
        bins_var = bins.copy()
        bins_var[bins_var == 0] = 1
        likelihood += (bins[i] - variance[0])**2/(bins_var[i])   
    
    return likelihood

def likelihood_poison(A, Nsat, point, edges, bins, radii): 
    a, b, c = point
    func_integrate = lambda x: 4*np.pi *x**2 * n(x, A, Nsat, a, b, c)
    result, result_error = romberg(xmin, xmax, 10, func_integrate)
    A_new = 1/(result[0])
    
    likelihood = 0
    for i in range(len(radii)): 
        likelihood -= np.log(N(radii[i], A_new, Nsat, a, b, c))
        

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

# def likelihood(A, Nsat, a_par, b_par, c_par, edges, bins):
#     likelihood_total = []
#     for i in range(len(a_par)):
#         for j in range(len(b_par)):
#             for k in range(len(c_par)):
#                 a, b, c = a_par[i], b_par[j], c_par[k]
        
#                 func_integrate = lambda x: 4* np.pi * x**2 * n(x, A, Nsat, a, b, c)
#                 result, result_error = romberg(xmin, xmax, 10, func_integrate)
#                 norm = 1/result[0]
#                 func_integrate = lambda x: 4* np.pi * x**2 * n(x, norm, Nsat, a, b, c)
#                 likelihood = 0
#                 for u in range(len(edges)-1): 
#                     var_integrate, error = romberg(edges[u], edges[u+1], 10, func_integrate)
#                     model = func_integrate(x = (edges[u] + edges[u+1])*0.5)
#                     likelihood += 0.5 * ((bins[u] - model)**2/var_integrate[0])
                
#                 likelihood_total += [likelihood]
#     return likelihood_total 

#Call this function as: 
Nsat_calc = []
radii = []
Nhalo = []
filenames = ['satgals_m11.txt', 'satgals_m12.txt', 'satgals_m13.txt', 'satgals_m14.txt', 'satgals_m15.txt']
for i in filenames: 
    radius, nhalo = readfile(i)
    radii += [radius]
    Nhalo += [nhalo]
    Nsat_calc += [nhalo/len(radius)]
    


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
    
    simplex = np.array([[2.4, 0.2, 1.6], 
                        [2.3, 0.1, 1.2], 
                        [2.0, 0.5, 1.4], 
                        [2.6, 1.1, 1.8]])
    
    # simplex = np.array([[1.31420016, 1.11502032, 3.13383544], 
    #                         [1.31420016, 1.11502032, 3.13383544], 
    #                         [1.31420016, 1.11502032, 3.13383544], 
    #                         [1.31420016, 1.11502032, 3.13383544]]) * np.random.uniform(0.95, 1.05, size=(4,3))
    
    # simplex = np.array([[1.5, 2.5, 0.2],
    #                         [1, 0.5, 1.8],
    #                         [1.8, 1.6, 1],
    #                         [0.8, 2.3, 5.5]]) * 10
    
    print('AHHHH')
    
    # Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for mass bin i integrated per radial bin
    binned_data=np.histogram(x_radii,bins=edges)[0]/Nhalo[i]/np.diff(edges)
    
    likelihood_calc = lambda point: likelihood(1, Nsat, point, edges, binned_data)
    result_downhill = downhill_simplex(simplex, likelihood_calc, 1000, 1e-10)
    print(result_downhill)
    best_a, best_b, best_c = result_downhill #[1.53874748,   0.47773451,   1.19624515]
    
    func_integrate = lambda x: 4* np.pi * x**2 * n(x, 1, Nsat, best_a, best_b, best_c)
    result_A, result_error_A = romberg(xmin, xmax, 10, func_integrate)
        
    x_array = np.logspace(-4, np.log10(5), 100)

    # print(result)
    
    row=i//2
    col=i%2
    ax[row,col].step(edges[:-1], binned_data, where='post', label='binned data')
    ax[row,col].plot(x_array, N(x_array, 1/result_A[0], Nsat, best_a, best_b, best_c), label='best-fit profile')

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
#     Nsat = Nsat_calc[i] # replace by actual appropriate number for mass bin i
#     x_radii = radii[i] # replace by actual data for mass bin i
    
#     # Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for mass bin i integrated per radial bin
#     binned_data=np.histogram(x_radii,bins=edges)[0]/Nsat
    
#     simplex = np.array([[2.4, 0.2, 1.6], 
#                         [2.3, 0.1, 1.2], 
#                         [2.0, 0.5, 1.4], 
#                         [2.6, 1.1, 1.8]])
    
#     likelihood_calc = lambda point: likelihood_poison(1, Nsat, point, edges, binned_data, x_radii)
#     result_downhill = downhill_simplex(simplex, likelihood_calc, 1000, 1e-10)
#     best_a, best_b, best_c = result_downhill

#     func_integrate = lambda x: 4* np.pi * x**2 * n(x, 1, Nsat, best_a, best_b, best_c)
#     result_A, result_error_A = romberg(xmin, xmax, 10, func_integrate)
    
#     x_array = np.logspace(-4, np.log10(5), 100)


    

#     row=i//2
#     col=i%2
#     ax[row,col].step(edges[:-1], binned_data, where='post', label='binned data')
#     ax[row,col].plot(x_array, N(x_array, 1/result_A[0], Nsat, best_a, best_b, best_c), label='best-fit profile')

#     # ax[row,col].step(edges[:-1], Ntilda, where='post', label='best-fit profile')
#     ax[row,col].set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{{11+i}}} M_{{\\odot}}/h$")
#     ax[row,col].set_ylim(ylim_min[i], ylim_max[i])


# ax[2,1].set_visible(False)
# plt.tight_layout()
# handles,labels=ax[2,0].get_legend_handles_labels()
# plt.figlegend(handles, labels, loc=(0.65,0.15))
# plt.savefig('my_solution_1c.png', dpi=600)


# BONUS: Monte Carlo resampled fits (1e)
# num_samples = 10 #replace by how many resamplings you can draw/fit in reasonable time
# fig1e, ax = plt.subplots()
# Nsat = 100 # replace by actual appropriate number for mass bin i
# x_radii = np.random.rand(10000) * (xmax-xmin) # replace by actual data for chosen mass bin
# binned_data=np.histogram(x_radii,bins=edges)[0]/Nsat
# ax.step(edges[:-1], binned_data, where='post', label='binned data')
# Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for chosen mass bin integrated per radial bin
# ax.step(edges[:-1], Ntilda, where='post', label='best-fit profiles', color="C1")
# for i in range(num_samples):
#     Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for chosen mass bin integrated per radial bin
#     ax.step(edges[:-1], Ntilda, where='post', color="C1")
# # Also plot the mean or median fitted profile here
# ax.set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{...}} M_{{\\odot}}/h$")
# plt.legend()
# plt.savefig('my_solution_1e_chisq.png', dpi=600)

# num_samples = 10 #replace by how many resamplings you can draw/fit in reasonable time
# fig1e, ax = plt.subplots()
# Nsat = 100 # replace by actual appropriate number for mass bin i
# x_radii = np.random.rand(10000) * (xmax-xmin) # replace by actual data for chosen mass bin
# binned_data=np.histogram(x_radii,bins=edges)[0]/Nsat
# ax.step(edges[:-1], binned_data, where='post', label='binned data')
# Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for chosen mass bin integrated per radial bin
# ax.step(edges[:-1], Ntilda, where='post', label='best-fit profiles', color="C2")
# for i in range(num_samples):
#     Ntilda = np.sort(np.random.rand(n_bins)) * (xmax-xmin) # replace by fitted model for chosen mass bin integrated per radial bin
#     ax.step(edges[:-1], Ntilda, where='post', color="C2")
# # Also plot the mean or median fitted profile here
# ax.set(yscale='log', xscale='log', xlabel='x', ylabel='N', title=f"$M_h \\approx 10^{{...}} M_{{\\odot}}/h$")
# plt.legend()
# plt.savefig('my_solution_1e_poisson.png', dpi=600)
