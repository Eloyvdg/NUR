import numpy as np
import matplotlib.pyplot as plt
import time

a1 = np.uint64(21)
a2 = np.uint64(35)
a3 = np.uint64(4)
a = np.uint64(4294957665)

class RNG():
    """Class for a random number generator"""
    def __init__(self):
        """Define the seed for XOR and Multiply With Carry"""
        self.state1 = time.time()
        self.state2 = time.time()

    @staticmethod
    def _XOR(x, a1=a1, a2=a2, a3=a3):
        """Use XOR to calculate a random number"""
        x = np.uint64(x)
        x = x ^ (x << a1)
        x = x ^ (x >> a2)
        x = x ^ (x << a3)
        return x

    @staticmethod
    def _Carry(x, a=a):
        """Use Multiply with Carry to calculate a random number"""
        x = np.uint64(x)
        x = a * (x & (np.uint64(2**(32)) - np.uint64(1))) + (x >> np.uint64(32))
        return np.uint32(x)

    def random_numbers(self, a_min, a_max):
        """Combine XOR and MWC to calculate a new random number
        a_min: lower limit for the random numbers
        a_max: upper limit for the random nmbers"""
        state1 = np.uint64(self.state1)
        state2 = np.uint64(self.state2)

        self.state1 = RNG._XOR(state1)
        self.state2 = RNG._Carry(state2)

        # Equation to give a lower and upper limit for the random numbers
        rand_num = a_min + (a_max - a_min) * np.uint32(state1 ^ state2) * 2**-32

        return rand_num
        
if __name__ == '__main__':
    # Check if the numbers are random
    initial = RNG()
    rand_number = []
    for i in range(10000):
        rand_number += [initial.random_numbers(10, 40)]

    plt.hist(rand_number, cumulative = True, density = True)
    plt.show()
