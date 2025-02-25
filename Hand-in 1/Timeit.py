mysetup = '''import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from Decomposition import solver
import Interpolator

data=np.genfromtxt(os.path.join(sys.path[0],"Vandermonde.txt"),comments='#',dtype=np.float64)

x=data[:,0]
y=data[:,1]
xx=np.linspace(x[0],x[-1],1001)
yya=np.zeros(1001,dtype=np.float64) #replace!
ya=np.zeros(len(x),dtype=np.float64) #replace! 
yyb=yya.copy() #replace!
yb=ya.copy() #replace!
yyc10=yya.copy() #replace!
yc10=ya.copy() #replace!'''

mycode_a = '''start = solver(x,y)
LU = start._crout(x)
b = start.solve()

for i in range(len(xx)): 
    yya[i] = np.dot(b, xx[i]**np.arange(0,20)) '''
    
mycode_b = '''
interp = Interpolator.Interpolation(x,y)
for i in range(len(xx)):
    yyb[i] = interp.neville(xx[i], 19) '''
    
mycode_c = '''start = solver(x,y)
b_improved_10 = start.iteration(10)
for i in range(len(xx)): 
    yyc10[i] = np.dot(b_improved_10, xx[i]**np.arange(0,20))'''


timeit_a = timeit.timeit(setup = mysetup, stmt = mycode_a, number = 100)
timeit_b = timeit.timeit(setup = mysetup, stmt = mycode_b, number = 100)
timeit_c = timeit.timeit(setup = mysetup, stmt = mycode_c, number = 100)

print(f'The duration of the the LU decomposition is {timeit_a} seconds')
print(f'The duration of the the Neville interpolation is {timeit_b} seconds')
print(f'The duration of the the LU decomposition with 10 iterations is {timeit_c} seconds')
