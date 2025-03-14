import numpy as np
import matplotlib.pyplot as plt
import time
import RNG
a1 = np.uint64(21)
a2 = np.uint64(35)
a3 = np.uint64(4)
a = np.uint64(4294957665)

seed = time.time() 

def XOR(x, a1=a1, a2=a2, a3=a3): 
    x = np.uint64(x)
    x = x ^ (x >> a1)
    x = x ^ (x << a2)
    x = x ^ (x >> a3)
    return x

def Carry(x, a=a):
    x = np.uint64(x)
    x = a * (x & (np.uint64(2**(32)) - np.uint64(1))) + (x >> np.uint64(32))
    return np.uint32(x)

def random_numbers( n, a_min, a_max, seed1 = time.time(), seed2 = time.time()): 
    state1 = np.uint64(seed1)
    state2 = np.uint64(seed2)
    
    rand_numbers = []
    for i in range(n):
        state1 = XOR(state1)
        state2 = Carry(state2)
        rand_num = a_min + (a_max - a_min) * np.uint32(state1 ^ state2) * 2**-32
    seed = rand_num
    rand_numbers += rand_num
    
    
    return rand_num

def trapezoid(h, array, function):   
    result = 0
    for i in range(len(array) -1): 
        result += h * 0.5 * (function(array[i]) + function(array[i+1]))
    return result

def simpson(array, h, function): 
    result = 0
    inv_3 = 1/3
    for i in range(1, len(array) -1):
        result += h * inv_3 * (function(array[i-1]) + 4*function(array[i]) + function(array[i+1]))


def Romberg(a, b, m, function):
    r = np.zeros(m)
    h = b - a
    r[0] = 0.5 * h * (function(a) + function(b))
    N_p = 1
    for i in range(1,m): 
        r[i] = 0 
        delta = h
        h = 0.5*h
        x = a + h
        for j in range(N_p): 
            r[i] += function(x)
            x += delta
            
        r[i] = 0.5 * (r[i-1] + delta * r[i])
        
        N_p *= 2    
    print(r)
    N_p = 1
    for i in range(1,m): 
        N_p *= 4
        for j in range(0, m-i): 
            r[j] = (N_p * r[j+1] - r[j]) / (N_p - 1)
                
    error = np.abs(r[0] - r[1])
    return r, error

def func(x):
    return x**2

A= 1. # to be computed
Nsat=100
a=2.4
b=0.25
c=1.6

def n(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    return A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)

def func_integrate(x, A=A, Nsat=Nsat, a=a, b=b, c=c):
    return (A*Nsat*((x/b)**(a-3))*np.exp(-(x/b)**c)) * x**2 * 4 * np.pi

result, result_error = Romberg(1e-6,5, 15, func_integrate)
A = 1 / result[0]

#Plot of histogram in log-log space with line (question 1b)
xmin, xmax = 10**-4, 5
N_generate = 10000

relative_radius = 10**np.linspace(np.log10(xmin), np.log10(xmax), 200) #replace!
analytical_function = func_integrate(relative_radius, A=A) #replace


# rand_x = 10**np.array(random_numbers(int(N_generate), -4, np.log10(5)))
# rand_y = np.array(random_numbers(int(N_generate), 0, np.max(analytical_function)))

accepted_x = []
accepted_len = 0
it = 0

new_state1 = time.time()
new_state2 = time.time() 

initial = RNG.RNG()

while accepted_len < N_generate: 
    it += 1
    if it > 1e8: 
        break
    # rand_x = 10**np.array(initial.random_numbers(-4, np.log10(xmax)))
    rand_x = np.array(initial.random_numbers(xmin, xmax))
    rand_y = np.array(initial.random_numbers(0, np.max(analytical_function)))
    # rand_x = 10**np.array(random_numbers(1, -4, np.log10(xmax), seed = new_state1))
    # rand_y = np.array(random_numbers(1, 0, np.max(analytical_function), seed = new_state1))
    if rand_y < func_integrate(rand_x, A = A): 
        accepted_x += [rand_x]
        accepted_len += 1
    else:
        continue

#21 edges of 20 bins in log-space
edges = 10**np.linspace(np.log10(xmin), np.log10(xmax), 21)
hist = np.histogram(accepted_x, bins=edges)[0] #replace!
hist_scaled = hist/len(accepted_x)/np.diff(edges) #replace; this is NOT what you should be plotting, this is just a random example to get a plot with reasonable y values (think about how you *should* scale hist)

fig1b, ax = plt.subplots()
ax.stairs(hist_scaled, edges=edges, fill=True, label='Satellite galaxies') #just an example line, correct this!
plt.plot(relative_radius, analytical_function, 'r-', label='Analytical solution') #correct this according to the exercise!


# plt.hist(accepted_x, density= True, bins = edges)

ax.set(xlim=(xmin, xmax), ylim=(10**(-3), 10), yscale='log', xscale='log',
       xlabel='Relative radius', ylabel='Number of galaxies')
ax.legend()
plt.savefig('my_solution_1b.png', dpi=600)

#Cumulative plot of the chosen galaxies (1c)
chosen = xmin + np.sort(np.random.rand(Nsat))*(xmax-xmin) #replace!
fig1c, ax = plt.subplots()
ax.plot(chosen, np.arange(100))
ax.set(xscale='log', xlabel='Relative radius', 
       ylabel='Cumulative number of galaxies',
       xlim=(xmin, xmax), ylim=(0, 100))
plt.savefig('my_solution_1c.png', dpi=600)
