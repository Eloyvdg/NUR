import numpy as np
import matplotlib.pyplot as plt
import time

a1 = np.uint64(21)
a2 = np.uint64(35)
a3 = np.uint64(4)
a = np.uint64(4294957665)

class RNG():
    def __init__(self):
        self.state1 = time.time()
        self.state2 = time.time()

    @staticmethod
    def _XOR(x, a1=a1, a2=a2, a3=a3):
        x = np.uint64(x)
        x = x ^ (x << a1)
        x = x ^ (x >> a2)
        x = x ^ (x << a3)
        # self.state1 = x
        return x

    @staticmethod
    def _Carry(x, a=a):
        x = np.uint64(x)
        x = a * (x & (np.uint64(2**(32)) - np.uint64(1))) + (x >> np.uint64(32))
        # self.state2 = x
        return np.uint32(x)

    def random_numbers(self, a_min, a_max):
        state1 = np.uint64(self.state1)
        state2 = np.uint64(self.state2)

        rand_numbers = []
        # for i in range(n):
        self.state1 = RNG._XOR(state1)
        self.state2 = RNG._Carry(state2)

        rand_num = a_min + (a_max - a_min) * np.uint32(state1 ^ state2) * 2**-32
        # seed = rand_num
        # rand_numbers += rand_num
        return rand_num
    
    @staticmethod
    def _fisher_yates(array, N):
        for i in range(len(array))[::-1]: 
            rand_int = RNG.random_numbers(0,i).astype(int)
        # 
            array[i], array[rand_int] = array[rand_int], array[i]
            
        return array[:N]
        

if __name__ == '__main__':
    initial = RNG()
    rand_number = []
    for i in range(10000):
        rand_number += [initial.random_numbers(10, 40)]

    plt.hist(rand_number, cumulative = True, density = True)
    plt.show()
