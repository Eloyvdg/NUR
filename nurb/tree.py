import numpy as np
from quicksort import quick_sort
import matplotlib.patches as patches
import matplotlib.pyplot as plt

class Node(): 
    def __init__(self, coord, depth, axis, left_child, right_child, box_min, box_max):
        self.axis = axis
        self.depth = depth
        self.coord = coord
        self.left_child = left_child
        self.right_child = right_child
        self.box_min = box_min # lower left and upper right
        self.box_max = box_max
    
    def __repr__(self): 
        return f'depth = {str(self.depth)}\n axis = {str(self.axis)}\n left child = {str(self.left_child)}\n median = {str(self.coord)}\n'
        
        
def build_tree(array, depth, max_depth, box_min, box_max): 
    if depth >= max_depth or len(array) == 0: 
        return None
    
    
    # box_min = np.min(array, axis = 0)
    # box_max = np.max(array, axis = 0)
    # print(box_min, box_max)
    
    
    k = array.shape[1]
    axis = (depth) % k
    # print(axis)

    N = len(array)
    array_1d = array[:,axis]
    array_sorted = quick_sort(array_1d)
    
    # median_coord = array_sorted[N//2]
    
    midpoint = 0.5 * (box_max[axis] + box_min[axis])
    # print(midpoint)
    
    # if axis == 0: 
    #     box_min = box_min
    #     box_max = box_max - 0.5 * (box_max - box_min)
    # elif axis == 1: 
    #     box_min = box_min + 0.5 * (box_max - box_min)
    #     box_max = box_max
        
    # print(box_min, box_max)
        
    # print(median_coord)
    # loc_coord = np.where(array_1d == median_coord)[0][0]
    mask = array[:,axis] < midpoint
    median_coord = 0
    
    split_left, split_right = array[mask, :], array[~mask,:]
    # split_left, split_right = array[:N//2, :], array[N//2:,:]
    # print(split_l1eft)

    
    # axis = (axis + 1) % k

    # depth += 1
    
    # print(axis)
    
    box_min_left = box_min.copy()
    box_max_left = box_max.copy()
    box_max_left[axis] = midpoint
    print(box_min_left, box_max_left, midpoint)
    
    box_min_right = box_min.copy()
    box_max_right = box_max.copy()
    box_min_right[axis] = midpoint
    
    node1 = build_tree(split_left, depth + 1, max_depth, box_min_left, box_max_left)
    node2 = build_tree(split_right, depth + 1, max_depth, box_min_right, box_max_right)
    
    # if axis == 0:
    #     box_min_new = box_min.copy()
    #     box_max_new = box_max.copy()
    #     box_min_new[1] = box_min_new[1] + midpoint
    #     box_max_new[1] = box_max_new[1] - midpoint
    #     print(box_min_new)
    #     node1 = build_tree(split_left, depth + 1, max_depth, box_min_new, box_max)
    #     node2 = build_tree(split_right, depth + 1, max_depth, box_min, box_max_new)
    # elif axis == 1: 
    #     box_min_new = box_min.copy()
    #     box_max_new = box_max.copy()
    #     box_min_new[0] = box_min[0] + 0.5 * (box_max[0] - box_min[0])
    #     box_max_new[0] = box_max[0] - 0.5 * (box_max[0] - box_min[0])
    #     node1 = build_tree(split_left, depth + 1, max_depth, box_min, box_max_new)
    #     node2 = build_tree(split_right, depth + 1, max_depth, box_min_new, box_max)
    
    

    node = Node(median_coord, depth, axis, node1, node2, box_min, box_max)
    # print(node)
    
    return node
    
 

fig, ax = plt.subplots(figsize = (6,6))

def plot_tree(node):
    try: 
        height = node.box_max[1] - node.box_min[1]
        width = node.box_max[0] - node.box_min[0]
    
        # print(height, width)
        # print(node.box_min[0], node.box_min[1])
        # print()
        rect = patches.Rectangle((node.box_min[0], node.box_min[1]), width, height, linewidth=1, edgecolor='black', facecolor='none')
        ax.add_patch(rect)
    
        plot_tree(node.left_child)
        plot_tree(node.right_child)
    except AttributeError: 
            pass
        
    

        
    
if __name__ == '__main__': 
    np.random.seed (42)
    n_particles = 100
    dim = 2
    positions = np.random.rand(n_particles, dim)
    
    tree = build_tree(positions, 0, 4, np.array([0,0]), np.array([1,1]))
    # tree = build_tree(positions, 0, 4, np.min(positions, axis = 0), np.max(positions, axis = 0))

    plot_tree(tree)
    plt.plot(positions[:,0], positions[:,1], '.')
    # print(tree)
    
    # start = Node(positions)
    # median = start._find_median(positions, axis = 0)
