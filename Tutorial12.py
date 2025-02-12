import numpy as np
import matplotlib.pyplot as plt
from Interpolator import Interpolation
import timeit

## Question 1
def factorial(n): 
    if n == 0 or n == 1:
        return 1
    x = 1
    for i in range(2, n+1):
        x *= i    
    return x
    
def sinc(x):
    return np.sin(x)/x

def sinc_2(x, n):
    tot = 0
    for i in range(n):
        tot += ( (-1)**i * x**(2*i) )/ (factorial(2*i + 1))
    return tot  


x_range = np.arange(0.1, 0.5, 0.01)

plt.plot(x_range, sinc(x_range))
plt.plot(x_range, sinc_2(x_range, 5), '--')
plt.show()

## Question 2 
G = 6.67e-11
c = 2.98e3
c_inv = 1/c
c_inv2 = c_inv * c_inv

def Rs(M):
    return (2 * G * M)/ (c**2)

def Rs2(M):
    return (2 * G * M) * c_inv2

mean, std = 10**6, 10**5
masses = np.random.normal(mean, std, 10000)

radius = timeit.timeit(lambda: Rs(masses), number = 10000)
radius2 = timeit.timeit(lambda: Rs2(masses), number = 10000)
print(radius, radius2)

## Question 3

from matplotlib.image import imread

image = imread('M42128.jpg')
row = image[0]
interp_x = np.linspace(0,128, 201)

# def bisection(xdata, x, n):
#     low = xdata[0]
#     high = xdata[-1]
#     k = low-1
#     for i in xdata:
#         if i > k:
#             k = i
#         else: 
#             print('Data is not monotonic')
    
#     while (high - low) > n:
#         midpoint = int((high + low)*0.5)
#         if x > midpoint: 
#             low = midpoint
#         else: 
#             high = midpoint
#     return int(low), (high)

x_test = np.arange(0, 7)
x_real = np.linspace(0, 6, 200)

x_real_2 = np.linspace(0,6, 200)
y_test = np.sin(x_test)
y_real = np.sin(x_real)
y_real_2 = np.sin(x_real_2)

interp = Interpolation(x_test, y_test)

interp_values = interp.linear(x_real_2)

plt.plot(x_test, y_test, '.')
plt.plot(x_real_2, interp_values)
plt.plot(x_real_2, y_real_2)
plt.show()

    
    
