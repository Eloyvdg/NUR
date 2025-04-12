import numpy as np

def quick_sort(array): 
    """
    Quicksort algorithm to sort an array of numbers
    array: array with numbers to sort
    """
    array = np.array(array)
    middle_idx = len(array)//2
    first, last, middle = array[0], array[-1], array[middle_idx]
    
    # Sort the first, middle and last numbers of the array
    if first >= middle:
        array[0], array[middle_idx], array[middle_idx], array[0]
    if middle >= last: 
        array[-1], array[middle_idx] = array[middle_idx], array[-1]
    if first >= last: 
        array[0], array[-1] = array[-1], array[0]
        
    
    pivot = array[middle_idx]
    
    if len(array) < 3: 
        # Array is already sorted 
        return
       
    i_flag, j_flag = False, False
    j = len(array) - 1
    i = 0
    
    # Loop over the array from the left and right and switch if needed
    while j > i: 
        
        if i_flag == False:
            if array[i] >= pivot:
                i_switch = i
                i_flag = True
            else:
                i += 1
        
        if j_flag == False:
            if array[j] <= pivot:
                j_switch = j
                j_flag = True
            else:
                j -= 1
            
        if i_flag * j_flag:
            if array[i_switch] != array[j_switch]:
                array[i_switch], array[j_switch] = array[j_switch], array[i_switch]
            else: 
                i += 1
            i_flag, j_flag = False, False
    
    loc_pivot = np.where(array == pivot)[0][0]

    
    # Sort the subarrays left and right from pivot
    if len(array[:loc_pivot]) > 1: 
        quick_sort(array[:loc_pivot])
        
    if len(array[loc_pivot:]) > 1: 
        quick_sort(array[loc_pivot+1:])
        
    return array
