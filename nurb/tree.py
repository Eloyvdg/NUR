import numpy as np
from quicksort import quick_sort
import matplotlib.patches as patches
import matplotlib.pyplot as plt

class Node(): 
    def __init__(self, depth, axis, child1, child2, child3, child4, box_min, box_max):
        self.axis = axis
        self.depth = depth
        self.child1 = child1
        self.child2 = child2
        self.child3 = child3
        self.child4 = child4
        self.box_min = box_min # lower left and upper right
        self.box_max = box_max
    
    def __repr__(self): 
        return f'depth = {str(self.depth)}\n axis = {str(self.axis)}\n left child = {str(self.left_child)}\n'
        
        
def build_tree(array, depth, max_depth, box_min, box_max): 
    if depth >= max_depth or len(array) < 3: 
        return None
    
    k = array.shape[1]
    axis = (depth) % k

    # N = len(array)
    # array_1d = array[:,axis]
    # array_sorted = quick_sort(array_1d)
        
    midpoint = 0.5 * (box_max[axis] + box_min[axis])

    mask = array[:,axis] < midpoint
    
    split_left, split_right = array[mask, :], array[~mask,:]
    
    box_min_left = box_min.copy()
    box_max_left = box_max.copy()
    box_max_left[axis] = midpoint
    
    box_min_right = box_min.copy()
    box_max_right = box_max.copy()
    box_min_right[axis] = midpoint
    
    node1 = build_tree(split_left, depth + 1, max_depth, box_min_left, box_max_left)
    node2 = build_tree(split_right, depth + 1, max_depth, box_min_right, box_max_right)
    node3 = build_tree(split_left, depth + 1, max_depth, box_min_left, box_min_right)
    node4 = build_tree(split_right, depth + 1, max_depth, box_max_left, box_max_right)

    node = Node(depth, axis, node1, node2, node3, node4, box_min, box_max)
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
    
        plot_tree(node.child1)
        plot_tree(node.child2)
        plot_tree(node.child3)
        plot_tree(node.child4)
    except AttributeError: 
            pass
        
    

        
    
if __name__ == '__main__': 
    np.random.seed (42)
    n_particles = 100
    dim = 2
    positions = np.random.rand(n_particles, dim)
    
    tree = build_tree(positions, 0, 4, np.array([0.0,0.0]), np.array([1.0,1.0]))
    # tree = build_tree(positions, 0, 4, np.min(positions, axis = 0), np.max(positions, axis = 0))

    plot_tree(tree)
    plt.plot(positions[:,0], positions[:,1], '.')
    # print(tree)
    
    # start = Node(positions)
    # median = start._find_median(positions, axis = 0)
    
    
