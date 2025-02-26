import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from Solve_matrix import solver
import Interpolator
import copy

#This script is to get you started with reading the data and plotting it
#You are free to change whatever you like/do it completely differently,
#as long as the results are clearly presented

data=np.genfromtxt(os.path.join(sys.path[0],"Vandermonde.txt"),comments='#',dtype=np.float64)

x=data[:,0]
y=data[:,1]
xx=np.linspace(x[0],x[-1],1001) #x values to interpolate at

#Insert your own code to calculate interpolated y values here!
#The plotting code below assumes you've given the interpolated
#values for 2a suffix "a", those for 2b "b", and those for 2c 
#"c1" and "c10" – feel free to change these
#Note that you interpolate to xx for the top panel but to
#x for the bottom panel (since we can only compare the accuracy
#to the given points, not in between them)
yya=np.zeros(1001,dtype=np.float64) #replace!
ya=np.zeros(len(x),dtype=np.float64) #replace!

# Use a single LU decomposition to find the 19th order polynomial
start = solver(x,y)
b = start.solve()
print('The resulting vector after solving with LU decomposition is' + r'\textbf{c}' + f'is {b}')

for i in range(len(xx)): 
    yya[i] = np.dot(b, xx[i]**np.arange(0,20))
    
for i in range(len(x)): 
    ya[i] = np.dot(b, x[i]**np.arange(0,20))

yyb=yya.copy() #replace!
yb=ya.copy() #replace!

# Perform Neville's interpolation to find the 19th order polynomial
interp = Interpolator.Interpolation(x,y)
for i in range(len(xx)):
    yyb[i] = interp.neville(xx[i], 19)
    
for i in range(len(x)): 
    yb[i] = interp.neville(x[i],19)

yyc1=yya.copy() #replace!
yc1=ya.copy() #replace!
yyc10=yya.copy() #replace!
yc10=ya.copy() #replace!

# Calculate the 19th order polynomial using LU decomposition, with 1 and 10 iterations respectively 
b_improved_1 = start.iteration(1)
b_improved_10 = start.iteration(10)

for i in range(len(xx)): 
    yyc1[i] = np.dot(b_improved_1, xx[i]**np.arange(0,20))
    
for i in range(len(x)): 
    yc1[i] = np.dot(b_improved_1, x[i]**np.arange(0,20))

for i in range(len(xx)): 
    yyc10[i] = np.dot(b_improved_10, xx[i]**np.arange(0,20))
    
for i in range(len(x)): 
    yc10[i] = np.dot(b_improved_10, x[i]**np.arange(0,20))

#Don't forget to output the coefficients you find with your LU routine
print('The vector for 1 iteration with LU decompoisition is ' + r'\textbf{c}' + f'is {b_improved_1}')
print('The vector for 10 iteration with LU decompoisition is ' + r'\textbf{c}' + f'is {b_improved_10}')

#Plot of points with absolute difference shown on a log scale (question 2a)
fig=plt.figure()
gs=fig.add_gridspec(2,hspace=0,height_ratios=[2.0,1.0])
axs=gs.subplots(sharex=True,sharey=False)
axs[0].plot(x,y,marker='o',linewidth=0)
plt.xlim(-1,101)
axs[0].set_ylim(-400,400)
axs[0].set_ylabel('$y$')
# axs[1].set_ylim(1e-16,1e1)
axs[1].set_ylabel('$|y-y_i|$')
axs[1].set_xlabel('$x$')
axs[1].set_yscale('log')
line,=axs[0].plot(xx,yya,color='orange')
line.set_label('Via LU decomposition')
axs[0].legend(frameon=False,loc="lower left")
axs[1].plot(x,abs(y-ya),color='orange')
plt.savefig('plots/my_vandermonde_sol_2a.png',dpi=600)

#For questions 2b and 2c, add this block
line,=axs[0].plot(xx,yyb,linestyle='dashed',color='green')
line.set_label('Via Neville\'s algorithm')
axs[0].legend(frameon=False,loc="lower left")
axs[1].plot(x,abs(y-yb),linestyle='dashed',color='green')
plt.savefig('plots/my_vandermonde_sol_2b.png',dpi=600)

#For question 2c, add this block too
line,=axs[0].plot(xx,yyc1,linestyle='dotted',color='red')
line.set_label('LU with 1 iteration')
axs[1].plot(x,abs(y-yc1),linestyle='dotted',color='red')
line,=axs[0].plot(xx,yyc10,linestyle='dashdot',color='purple')
line.set_label('LU with 10 iterations')
axs[1].plot(x,abs(y-yc10),linestyle='dashdot',color='purple')
axs[0].legend(frameon=False,loc="lower left")
plt.savefig('plots/my_vandermonde_sol_2c.png',dpi=600)

#Don't forget to caption your figures to describe them/
#mention what conclusions you draw from them!
