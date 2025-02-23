import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread

    
class Interpolation:
    
    def __init__(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata
        
    def _bisection(self, x): 
        low = 0
        high = len(self.xdata) - 1
        k = low-1
        for i in self.xdata:
            if i > k:
                k = i
            else: 
                print('Data is not monotonic')
        
        while (high - low) > 1:
            midpoint = (low + high) // 2
            if x > self.xdata[midpoint]: 
                low = midpoint
            else: 
                high = midpoint
        return low, high
    
    def linear(self, xdata):
        interp_values = []
        for i in xdata:
            idx1, idx2 = self._bisection(i)


            ylow, yhigh = self.ydata[idx1], self.ydata[idx2]
            interp_values += [ylow + (yhigh-ylow)/(self.xdata[idx2]-self.xdata[idx1]) * (i - self.xdata[idx1])]
        return interp_values
    
    def neville(self, x, order):
        
        M = order + 1
        idx1, idx2 = self._bisection(x)

        split = M // 2        

        if idx1 < split: 
            idx2 = M
            idx1 = 0
        elif idx2 > len(self.xdata) - split: 
            idx1 = len(self.xdata) - M 
            idx2 = len(self.xdata)
        else: 
            idx1 -= split
            idx2 += split
    
        P = self.ydata[idx1:idx2].copy()
        for k in range(1,M):
            for i in range(M-k):
                j = i + k
                P[i] = ( (self.xdata[idx1:idx2][j] - x) * P[i] + (x - self.xdata[idx1:idx2][i]) * P[i+1]) / (self.xdata[idx1:idx2][j] - self.xdata[idx1:idx2][i])
        
        return P[0]
    
if __name__ == '__main__':
    x_sin = np.arange(0,12)
    y_sin = np.sin(x_sin)
    
    x_sin_interp = np.linspace(0,11, 200)
    interp_sin = Interpolation(x_sin, y_sin)
    y_sin_interp = interp_sin.linear(x_sin_interp)
    
    P = []
    for i in x_sin_interp: 
        P += [interp_sin.neville(i, 3)]
    
    plt.figure(figsize = (8,6))
    plt.plot(x_sin, y_sin, '.', label = 'Data', zorder = 3)
    plt.plot(x_sin_interp, y_sin_interp, label = 'Linear interpolation', zorder = 2)
    plt.plot(x_sin_interp, np.sin(x_sin_interp), label = 'Analytical', zorder = 1)
    plt.plot(x_sin_interp, P, label = 'Neville')
    plt.legend()
    plt.show()
    
    data=np.genfromtxt(("Vandermonde.txt"),comments='#',dtype=np.float64)
    x=data[:,0]
    y=data[:,1]
    xx=np.linspace(x[0],x[-1],1001) #x values to interpolate at

    interp = Interpolation(x,y)
    P = []
    for i in xx:
        P += [interp.neville(i, 5)]

    plt.figure(figsize = (8,6))
    plt.plot(x, y, '.')
    plt.plot(xx, P)  
