import numpy as np

def poisson_log(lamb, k):
    result = k * np.log(lamb) - lamb
    log_fact = 0
    for i in range(1, k + 1):
        log_fact += np.log(k)
    return np.exp(result - log_fact)

lamb = np.array([1, 5, 3, 2.6, 100, 101]).astype(np.float32)
k = np.array([0, 10, 21, 40, 5, 200]).astype(np.int32)
P = []
for j in range(len(k)): 
    P += [poisson_log(lamb[j], k[j])]
