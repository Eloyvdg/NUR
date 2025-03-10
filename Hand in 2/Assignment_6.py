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
    first, last, middle = array[0], array[-1], array[len(array)//2]
    print(first, middle, last, 'locs')
    pivot = np.median([first, middle, last]).astype(array.dtype)
    
    loc_pivot = np.where(array == pivot)[0][0]
    array[len(array)//2], array[loc_pivot] = pivot, array[len(array)//2]
    i_flag, j_flag = False, False
    
    print(pivot, 'pivot')
    print(array, 'pivot replacing')
    
    for i,j in zip(range(len(array)), range(len(array))[::-1]):
        
        if j <= i: 
            break
        
        if array[i] >= pivot:
            i_switch = i
            i_flag = True
            
        if array[j] <= pivot:
            j_switch = j
            j_flag = True
            
        if i_flag * j_flag:
            array[i_switch], array[j_switch] = array[j_switch], array[i_switch]
            print(array, i_switch, j_switch, 'switching')
            i_flag, j_flag = False, False
            #break
            
        # if j <= i: 
        #     # if i_flag == True: 
        #     #     array[i_switch], array[loc_pivot] = array[loc_pivot], array[i_switch]
        #     # if j_flag == True: 
        #     #     array[j_switch], array[loc_pivot] = array[loc_pivot], array[j_switch]
        #     break
    
    loc_pivot = np.where(array == pivot)[0][0]

    if len(array[:loc_pivot]) > 1: 
        quick_sort(array[:loc_pivot])
        
    if len(array[:loc_pivot+1]) > 1: 
        quick_sort(array[loc_pivot+1:])
        
        
    # print(array)
    return array
        # if array[i] > pivot:
        #     print(array[i])
        #     for j in range(len(array)-1, loc_pivot-1)[::-1]:
        #         if array[j] < pivot:
        #             array[j], array[i] = array[i], array[j]
                
    # print(array)
    # print(first, middle, last)
    # print(pivot)
    
array_test = np.array([2, 257, 6, 3890, 890, 0, 45, 4])
# np.random.shuffle(array_test)

# print(array_test)
# print(selection_sort(array_test))
array_test2 = np.arange(0,10)
np.random.shuffle(array_test2)

print(array_test2)
sort = quick_sort(array_test)
print(sort)

check = all(sort[i] <= sort[i+1] for i in range(len(sort) - 1))

print('Sorted?', check)