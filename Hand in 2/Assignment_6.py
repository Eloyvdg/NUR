import numpy as np

def selection_sort(array): 
    
    for i in range(len(array)-1):
        i_min = i
        for j in range(i+1, len(array)):
            if array[j] < array[i_min]: 
                i_min = j 
                
        if i_min != i:
                array[i_min], array[i] = array[i], array[i_min]
    
    return array

def quick_sort(array): 
    
    middle_idx = len(array)//2
    first, last, middle = array[0], array[-1], array[middle_idx]
    
    if first >= middle:
        array[0], array[middle_idx], array[middle_idx], array[0]
    if middle >= last: 
        array[-1], array[middle_idx] = array[middle_idx], array[-1]
    if first >= last: 
        array[0], array[-1] = array[-1], array[0]
        
    pivot = array[middle_idx]
    
    if len(array) < 3: 
        return
    
    
    i_flag, j_flag = False, False
    j = len(array) - 1
    i = 0
    
    while j > i: 
        
        if i_flag == False:
            if array[i] > pivot:
                i_switch = i
                i_flag = True
            else:
                i += 1
        
        if j_flag == False:
            if array[j] < pivot:
                j_switch = j
                j_flag = True
            else:
                j -= 1
            
        if i_flag * j_flag:
            array[i_switch], array[j_switch] = array[j_switch], array[i_switch]
            i_flag, j_flag = False, False
    
    loc_pivot = np.where(array == pivot)[0][0]


    if len(array[:loc_pivot]) > 1: 
        quick_sort(array[:loc_pivot])
        
    if len(array[loc_pivot:]) > 1: 
        quick_sort(array[loc_pivot+1:])
        
    return array

    
array_test = np.array([2, 257, 6, 6, 3890, 890, 0, 45, 4])

array_test2 = np.arange(0,1e6)
np.random.shuffle(array_test2)

print(array_test)
sort = quick_sort(array_test)
print(sort)

check = all(sort[i] <= sort[i+1] for i in range(len(sort) - 1))

print('Sorted?', check)
