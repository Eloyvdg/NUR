import numpy as np
import matplotlib.pyplot as plt

def bisection(xdata, x, n):
    low = xdata[0]
    high = xdata[-1]
    k = low-1
    for i in xdata:
        if i > k:
            k = i
        else: 
            print('Data is not monotonic')
    
    while (high - low) > n:
        midpoint = int((high + low)*0.5)
        if x > midpoint: 
            low = midpoint
        else: 
            high = midpoint
    return low, high
    
class Interpolation:
    def __init__(self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata
    
    def linear(self, xdata):
        interp_values = []
        for i in xdata:
            idx1, idx2 = bisection(np.arange(0, len(self.xdata)), i, 1)
            ylow, yhigh = self.ydata[idx1], self.ydata[idx2]
            interp_values += [ylow + (yhigh-ylow)/(idx2-idx1) * (i - idx1)]
        return interp_values
    
    # x_test = np.arange(0, 7)
    # x_real = np.linspace(0, 6, 200)

    # x_real_2 = np.linspace(0,6, 200)
    # y_test = np.sin(x_test)
    # y_real = np.sin(x_real)
    # y_real_2 = np.sin(x_real_2)
    # interp_values = []

    # for i in x_real: 
    #     pix_low, pix_high = bisection(np.arange(0, len(x_test)), i, 1)
    #     interp_values += [y_test[pix_low] + ((y_test[pix_high] - y_test[pix_low])/(pix_high - pix_low)) * (i - pix_low)]
    #     # interp_xvalues += [np.arange(pix_low, pix_high)]



    # plt.plot(x_test, y_test, '.')
    # plt.plot(x_real, interp_values)
    # plt.plot(x_real_2, y_real_2)
    # plt.show()
