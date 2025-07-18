import numpy as np
import matplotlib.pyplot as plt

def FFT(array): 
    # Function to calculare the FFT of an array
    array = np.asarray(array, dtype = np.complex64)
    N = len(array)
    
    half_N = N // 2
    N_inv = 1/N
    
    if N > 2: 
        # Split the array in even and odd
        even = FFT(array[::2])
        odd = FFT(array[1::2])
        array[:half_N] = even
        array[half_N:] = odd

        
    comb = np.zeros(N, dtype = np.complex64)
    
    for k in range(0, half_N):
        # Calculate the FFT elements
        t = array[k]
        factor = np.exp(2j*np.pi*k*N_inv) * array[k + N //2]
        comb[k] = t + factor
        comb[k + half_N] = t - factor
            
    return comb

def inv_FFT(array): 
    # Function to calculare the inverse FFT of an array
    array = np.asarray(array, dtype = np.complex64)
    N = len(array)
    half_N = N // 2
    N_inv = 1/N
    
    if N > 2: 
        # Split the array in even and odd
        even = inv_FFT(array[::2])
        odd = inv_FFT(array[1::2])
        array[:half_N] = even
        array[half_N:] = odd

        # array = np.concatenate((even, odd))
        
    comb = np.zeros(N, dtype = np.complex64)
    
    for k in range(0, half_N):
        # Calculate the iFFT elements
        t = array[k]
        factor = np.exp(-2j*np.pi*k*N_inv) * array[k + N //2]
        comb[k] = (t + factor) * 0.5
        comb[k + half_N] = (t - factor) * 0.5
            
    return comb
