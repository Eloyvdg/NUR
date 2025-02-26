import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import imread

    
class Interpolation:
    ''' Class to perform linear or Neville's interpolation
    Attributes: 
    xdata: List of known x data
    ydata: List of corresponding y data '''
    
    def __init__(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata
        
    def _bisection(self, x): 
        ''' Performs a bisection to a given x.
        Input: 
        x: x value to perform bisection at
        Output: 
        low, high: The closest points left and right of x ''' 
        
        low = 0
        high = len(self.xdata) - 1
        
        k = low-1
        for i in self.xdata: # Check if data is monotonic
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
        ''' Perform linear interpolation.
        Input: 
        xdata: all the x values to interpolate at
        Output: 
        interp_values: The corresponding interpolated y values '''
        
        interp_values = []
        for i in xdata:
            idx1, idx2 = self._bisection(i)
            ylow, yhigh = self.ydata[idx1], self.ydata[idx2]
            interp_values += [ylow + (yhigh-ylow)/(self.xdata[idx2]-self.xdata[idx1]) * (i - self.xdata[idx1])]
        return interp_values
    
    def neville(self, x, order):
        '''
        Perform interpolation according to Neville's theorem.
        Input: 
        x: all the x values to interpolate at
        Output: 
        The corresponding interpolated y values '''

        M = order + 1
        idx1, idx2 = self._bisection(x)

        split = M // 2        

        # Check if there are enough real points left or right
        if idx1 < split: 
            idx2 = M
            idx1 = 0
        elif idx2 > len(self.xdata) - split: 
            idx1 = len(self.xdata) - M 
            idx2 = len(self.xdata)
        else: 
            idx1 -= split
            idx2 += split

        # Calculate the interpolated values
        P = self.ydata[idx1:idx2].copy()
        for k in range(1,M):
            for i in range(M-k):
                j = i + k
                P[i] = ( (self.xdata[idx1:idx2][j] - x) * P[i] + (x - self.xdata[idx1:idx2][i]) * P[i+1]) / (self.xdata[idx1:idx2][j] - self.xdata[idx1:idx2][i])
        
        return P[0]    
