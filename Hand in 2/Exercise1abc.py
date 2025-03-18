import numpy as np
import matplotlib.pyplot as plt
import time
import RNG
import Integration

A= 1. # to be computed
Nsat=100
a=2.4
b=0.25
c=1.6

def quick_sort(array): 
    """
    Quicksort algorithm to sort an array of numbers
    array: array with numbers to sort
    """
    middle_idx = len(array)//2
    first, last, middle = array[0], array[-1], array[middle_idx]
    
    # Sort the first, middle and last numbers of the array
    if first >= middle:
        array[0], array[middle_idx], array[middle_idx], array[0]
    if middle >= last: 
        array[-1], array[middle_idx] = array[middle_idx], array[-1]
    if first >= last: 
        array[0], array[-1] = array[-1], array[0]
        
    
    pivot = array[middle_idx]
    
    if len(array) < 3: 
        # Array is already sorted 
        return
       
    i_flag, j_flag = False, False
    j = len(array) - 1
    i = 0
    
    # Loop over the array from the left and right and switch if needed
    while j > i: 
        
        if i_flag == False:
            if array[i] >= pivot:
                i_switch = i
                i_flag = True
            else:
                i += 1
        
        if j_flag == False:
            if array[j] <= pivot:
                j_switch = j
                j_flag = True
            else:
                j -= 1
            
        if i_flag * j_flag:
            if array[i_switch] != array[j_switch]:
                array[i_switch], array[j_switch] = array[j_switch], array[i_switch]
            else: 
                i += 1
            i_flag, j_flag = False, False
    
    loc_pivot = np.where(array == pivot)[0][0]

    
    # Sort the subarrays left and right from pivot
    if len(array[:loc_pivot]) > 1: 
        quick_sort(array[:loc_pivot])
        
    if len(array[loc_pivot:]) > 1: 
        quick_sort(array[loc_pivot+1:])
        
    return array

 

def func_integrate(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    """
    Function to integrate
    """
    return (A*((x/b)**(a-3))*np.exp(-(x/b)**c)) * x**2 * 4 * np.pi

xmin, xmax = 10**-4, 5 # boundaries of integral

integration_init = Integration.Integrator(func_integrate, xmin, xmax)
result, result_error = integration_init.romberg(15)
A = 1 / result[0] # Normalization factor
print(f'The normalization factor A is: {A}')

N_generate = 10000

relative_radius = 10**np.linspace(np.log10(xmin), np.log10(xmax), 200) #replace!
analytical_function = func_integrate(relative_radius, A=A)

accepted_x = []
accepted_len = 0
it = 0

new_state1 = time.time()
new_state2 = time.time() 

initial = RNG.RNG()

# Rejection sampling till N_generate samples are accepted 
while accepted_len < N_generate: 
    it += 1
    if it > 1e8: 
        break
    rand_x = np.array(initial.random_numbers(xmin, xmax))
    rand_y = np.array(initial.random_numbers(0, np.max(analytical_function)))
    if rand_y < func_integrate(rand_x, A = A): # Accept if smaller than analytical P(x)
        accepted_x += [rand_x]
        accepted_len += 1
    else:
        continue

# Plot the results of the rejections sampling
#21 edges of 20 bins in log-space
edges = 10**np.linspace(np.log10(xmin), np.log10(xmax), 21)
hist = np.histogram(accepted_x, bins=edges)[0]
hist_scaled = hist/len(accepted_x)/np.diff(edges)

fig1b, ax = plt.subplots()
ax.stairs(hist_scaled, edges=edges, fill=True, label='Satellite galaxies') #just an example line, correct this!
plt.plot(relative_radius, analytical_function, 'r-', label='Analytical solution') #correct this according to the exercise!

ax.set(xlim=(xmin, xmax), ylim=(10**(-3), 10), yscale='log', xscale='log',
       xlabel='Relative radius', ylabel='Number of galaxies')
ax.legend()
plt.savefig('my_solution_1b.png', dpi=600)


def fisher_yates(array, N):
    """
    Shuffle the array using Fisher-Yates algorithm
    array: array to be shuffled
    N: amount of samples to select
    """
    for i in range(len(array))[::-1]: 
        rand_int = int(RNG.RNG().random_numbers(0,i))
        array[i], array[rand_int] = array[rand_int], array[i]
    return array[:N]

selection = fisher_yates(np.array(accepted_x), 100)
sorted_selection = quick_sort(selection) # Sort using quicksort

# Plot the results.
#Cumulative plot of the chosen galaxies (1c)
chosen = xmin + np.sort(np.random.rand(Nsat))*(xmax-xmin) #replace!
fig1c, ax = plt.subplots()
ax.plot(sorted_selection, np.arange(100))
ax.set(xscale='log', xlabel='Relative radius', 
       ylabel='Cumulative number of galaxies',
       xlim=(xmin, xmax), ylim=(0, 100))
plt.savefig('my_solution_1c.png', dpi=600)
