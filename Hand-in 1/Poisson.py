import numpy as np

def poisson(lamb, k): # Calculate Poisson in logarithmic space
    result = k * np.log(lamb).astype(np.float32) - lamb
    log_fact = 0
    for i in range(1, k + 1):
        log_fact += np.log(i).astype(np.float32)
    return np.exp((result - log_fact).astype(np.float32))

lamb = np.array([1, 5, 3, 2.6, 100, 101]).astype(np.float32)
k = np.array([0, 10, 21, 40, 5, 200]).astype(np.int32)
P_k = []

for j in range(len(k)): 
    P_k += [poisson(lamb[j], k[j])]

# Print the results
for i in range(len(P_k)):
    print(f'For $\\lambda$ = {lamb[i]:.1e} and k = {int(k[i])}, P = {P_k[i]:.6e}')
