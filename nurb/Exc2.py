import numpy as np
import h5py
import matplotlib.pyplot as plt

# Question 2: Calculating potentials

# with h5py.File("/disks/cosmodm/DMO_a0.1_256.hdf5","r") as handle:
#     pos=handle["Position"][...] #particle positions, shape (Np,3), comoving
    #vel=handle["Velocity"][...] #particle velocities, shape (Np,3), comoving <-- not used, but if you're interested
    
with h5py.File("DMO_a0.1_256.hdf5","r") as handle:
    pos=handle["Position"][...] #particle positions, shape (Np,3), comoving

Np=np.int64(256)**3 #number of particles
mp=np.float32(3.64453e10) #particle mass in Msun; all 32-bit to save memory
G=np.float32(4.3009e-9) #gravitational constant in Mpc*(km/s)^2/Msun
h=np.float32(0.3755) #Hubble parameter (this is a Einstein-de Sitter universe with Omega_m=1)
L=np.float32(250.0) #side length of periodic cubic simulation volume
scale_factor=np.float32(0.1) #scale factor a
redshift=1.0/scale_factor-1
rho_mean=Np*mp/L**3 #mean density in Msun/Mpc^3 (comoving, matches 3*H_0^2/(8*pi*G))

# Question 2a: using Barnes-Hut [note: not actually calculating a potential, unless you do the bonus question]

# TO DO: build an octree (use a class for a node, so it can also refer to child nodes; avoid using lists for anything or the memory will balloon)

# Plotting the mass distribution for a slice
        
class Node(): 
    def __init__(self, depth, child1, child2, child3, child4, child5, child6, child7, child8, 
                 box_center, length, start_index, mass):
        self.depth = depth
        self.child1 = child1
        self.child2 = child2
        self.child3 = child3
        self.child4 = child4
        self.child5 = child5
        self.child6 = child6
        self.child7 = child7
        self.child8 = child8
        
        self.mass = mass
        self.length = length
        self.start_index = start_index
        self.box_center = box_center
        
index_array = np.arange(0,Np)
    
def build_tree(array, depth, max_depth, box_center, particles, start_index, length): 
    if depth > max_depth or len(array) == 0: 
        return None
    
    start_index = int(start_index)
    length = int(length)
    
    start_indices = np.zeros(8)
    start_indices[0] = start_index
    len_array = np.zeros(8)
    octant_indices = [None]*8
    
    indices_check = index_array[start_index:start_index+length]
    particles_check = particles[indices_check,:]

    
    mask = particles_check >= box_center
    mask = np.where(mask, 1, 0)
    
    index_octant = (mask[:,0] << 2) | (mask[:,1]) << 1| (mask[:,2] << 0)
    
    # new_indices = np.zeros(length)
    # ind = 0
    for j in range(8): 
        mask2 = index_octant == j
        new_ind = indices_check[mask2]
        len_array[j] = len(new_ind)
        octant_indices[j] = new_ind
        # start_indices[j] = int(start_index + np.sum(len_array[:ind]))
        # index_array[int(start_indices[j]): int(start_indices[j] + len_array[j])] = new_ind
        # ind+=1
        
    for k in range(1,8): 
        start_indices[k] = start_indices[k-1] + len_array[k-1]
        
    for l in range(8): 
        index_array[int(start_indices[l]):int(start_indices[l]+len_array[l])] = octant_indices[l]

    
    size = L/(2**depth)
    offsets = [-size/4, size/4]
    new_midpoints = np.zeros((3,8))
    ind = 0

    # new_indices = np.zeros(length)
    for i in range(2):
        for j in range(2): 
            for k in range(2):
                box_center_new = box_center + np.array([offsets[i], offsets[j], offsets[k]])
                new_midpoints[:,ind] = box_center_new
                ind += 1
        
    
    node1 = build_tree(array, depth + 1, max_depth, new_midpoints[:,0], particles, start_indices[0], len_array[0])
    node2 = build_tree(array, depth + 1, max_depth, new_midpoints[:,1], particles, start_indices[1], len_array[1])
    node3 = build_tree(array, depth + 1, max_depth, new_midpoints[:,2], particles, start_indices[2], len_array[2])
    node4 = build_tree(array, depth + 1, max_depth, new_midpoints[:,3], particles, start_indices[3], len_array[3])
    node5 = build_tree(array, depth + 1, max_depth, new_midpoints[:,4], particles, start_indices[4], len_array[4])
    node6 = build_tree(array, depth + 1, max_depth, new_midpoints[:,5], particles, start_indices[5], len_array[5])
    node7 = build_tree(array, depth + 1, max_depth, new_midpoints[:,6], particles, start_indices[6], len_array[6])
    node8 = build_tree(array, depth + 1, max_depth, new_midpoints[:,7], particles, start_indices[7], len_array[7])

    node = Node(depth, node1, node2, node3, node4, node5, node6, node7, node8, box_center, length, start_index, length*mp)
    
    return node


def traverse_tree(tree, depth, x_level): 
    
    pixels = 2**depth
    mapping = np.zeros((pixels, pixels))
    center = tree.box_center
    
    x_middle = center[0]/pixels
    # print(x_middle)
    def traverse(n): 
        if n == None: 
            return
        
        if n.depth == depth: 

            if n.box_center[0] == (x_middle + (2 * (x_level-1) * x_middle)):
                y = n.box_center[1]
                z = n.box_center[2]
                ind_z = int((z+x_middle)/(x_middle*2) -1)
                ind_y = int((y+x_middle)/(x_middle*2) -1)
                mapping[ind_y, ind_z] = n.mass
                
                return
        
        if n.depth <= depth:
            for child in [n.child1, n.child2, n.child3, n.child4, n.child5, n.child6, n.child7, n.child8]: 
                traverse(child)
            
    traverse(tree)
    return mapping
                
                
    
    

half_L = 0.5 * L
center_start = np.array([half_L, half_L, half_L]) 
tree =  build_tree(np.zeros(8), 0, 7, center_start, pos, 0, Np)

first_octant = tree.child4
start_ind = first_octant.start_index
length = first_octant.length
pos_first = pos[index_array[start_ind:start_ind+length],:]

for level in [3,5,7]: #feel free to change any of this code
    pixels=2**level
    size = L/pixels
    massmap=np.zeros((4,pixels,pixels),dtype=np.float32)
    for i in range(1,5): 
        massmap[i-1,:,:] = traverse_tree(tree, level, i)
    
    # TO DO: traverse the octree, fill map massmap[0,:,:] with the masses of nodes at depth 3 and x_index=x_0,
    #        massmap[1,:,:] with the masses of nodes at depth 3 and x_index=x_1, etc; then plot these slices;
    #        then do the same for levels 5 and 7

    fig, ax = plt.subplots(2,2, figsize=(10,8))
    pcm = ax[0,0].pcolormesh(np.linspace(0,L, pixels),np.linspace(0,L, pixels), massmap[0,:,:])
    ax[0,0].set(ylabel='z (Mpc)', title=f'x [0-{size}]')
    fig.colorbar(pcm, ax=ax[0,0], label='Total mass inside node')
    pcm =ax[0,1].pcolormesh(np.linspace(0,L, pixels), np.linspace(0,L, pixels), massmap[1,:,:])
    ax[0,1].set(xlabel='y (Mpc)', title=f'x [{size}-{2*size}]')
    fig.colorbar(pcm, ax=ax[0,1], label='Total mass inside node')
    pcm = ax[1,0].pcolormesh(np.linspace(0,L, pixels), np.linspace(0,L, pixels), massmap[2,:,:])
    ax[1,0].set(ylabel='z (Mpc)', title=f'x [{2*size}-{3*size}]')
    fig.colorbar(pcm, ax=ax[1,0], label='Total mass inside node')
    pcm = ax[1,1].pcolormesh(np.linspace(0,L, pixels), np.linspace(0,L, pixels), massmap[3,:,:])
    ax[1,1].set(xlabel='y (Mpc)', title=f'x [{3*size}-{4*size}]')
    fig.colorbar(pcm, ax=ax[1,1], label='Total mass inside node')
    ax[0,0].set_aspect('equal', 'box')
    ax[0,1].set_aspect('equal', 'box')
    ax[1,0].set_aspect('equal', 'box')
    ax[1,1].set_aspect('equal', 'box')
    # plt.savefig(f"fig2a_level{level}.png",dpi=300)
    plt.show()
    plt.close()

# # Question 2b: using the FFT

def FFT(array): 
    array = np.complex64(array)
    N = len(array)   

    if (N & (N-1) != 0) and N != 0:
        print((N & (N-1) == 0), N!= 0)
        raise Exception('Array does not have a length of power 2')
    
    def FFT_recursive(array):
        N = len(array)
        if N > 2: 
            even = array[::2]
            uneven = array[1::2]
            
            even = FFT_recursive(even)
            uneven = FFT_recursive(uneven)
            
            # array = np.concatenate((even, uneven))
            
        N_inv = 1/N
        half_N = int(0.5*N)

        comb = np.zeros(N, dtype = np.complex64)
        for k in range(0, half_N -1):
            # print(k)
            t = even[k]
            factor = np.exp(2j*np.pi*k*N_inv) * uneven[k]
            comb[k] = t + factor
            comb[k + half_N] = t - factor
        
        return comb
            
    array = FFT_recursive(array)
    
    # half_N = int(0.5 * len(array))
    # first_half = array[:half_N]
    # second_half = array[half_N:]

    # FFT_result = np.concatenate((second_half, first_half))
    
    return array

def inv_FFT(array): 
    array = np.complex64(array)
    N = len(array)   

    if (N & (N-1) != 0) and N != 0:
        print((N & (N-1) == 0), N!= 0)
        raise Exception('Array does not have a length of power 2')
    
    def FFT_recursive(array):
        N = len(array)
        if N > 2: 
            even = array[::2]
            uneven = array[1::2]
            
            even = FFT_recursive(even)
            uneven = FFT_recursive(uneven)
            
            # array = np.concatenate((even, uneven))
            

        N_inv = 1/N
        half_N = int(0.5*N)

        comb = np.zeros(N, dtype = np.complex64)
        for k in range(0, half_N -1):
            # print(k)
            t = even[k]
            factor = np.exp(-2j*np.pi*k*N_inv) * uneven[k] / 2
            comb[k] = t + factor
            comb[k + half_N] = t - factor
        
        return comb
            
    array = FFT_recursive(array)
    
    # half_N = int(0.5 * len(array))
    # first_half = array[:half_N]
    # second_half = array[half_N:]

    # FFT_result = np.concatenate((second_half, first_half))
    
    return array
    

def grid(tree, depth): 
    
    pixels = 2**depth
    grid = np.zeros((pixels, pixels, pixels))
    center = tree.box_center
    
    x_middle = center[0]/pixels
    def traverse(n): 
        if n == None: 
            return
        if n.depth == depth:
            x = n.box_center[0]
            y = n.box_center[1]
            z = n.box_center[2]
            ind_z = int((z+x_middle)/(x_middle*2) -1)
            ind_y = int((y+x_middle)/(x_middle*2) -1)
            ind_x = int((x+x_middle)/(x_middle*2) -1)
            grid[ind_x, ind_y, ind_z] = n.mass /(x_middle**3)
            
            return
        
        if n.depth <= depth:
            for child in [n.child1, n.child2, n.child3, n.child4, n.child5, n.child6, n.child7, n.child8]: 
                traverse(child)
            
    traverse(tree)
    return grid

Ngrid = 128
densgrid=grid(tree, 7)
potential= densgrid.copy()
potential = np.complex64(potential)

for i in range(densgrid.shape[0]): 
    for j in range(densgrid.shape[0]): 
        FFT_array = FFT(potential[i,j,:])
        potential[i,j,:] = FFT_array

for i in range(densgrid.shape[0]): 
    for j in range(densgrid.shape[0]): 
        FFT_array = FFT(potential[:,i,j])
        potential[:,i,j] = FFT_array
        
for i in range(densgrid.shape[0]): 
    for j in range(densgrid.shape[0]): 
        FFT_array = FFT(potential[i,:,j])
        potential[i,:,j] = FFT_array
        
val = 1.0/(Ngrid * 250/Ngrid)
N = (Ngrid-1)//2 + 1
k = np.zeros(Ngrid)
k[:N] = np.arange(0, N)*2*np.pi
k[N:] = np.arange(-(Ngrid//2), 0)*2*np.pi

kx, ky, kz = np.meshgrid(k, k, k)
k_squared = kx**2 + ky**2 + kz**2
k_squared[0,0,0] = 100
potential = potential/k_squared 

for i in range(densgrid.shape[0]): 
    for j in range(densgrid.shape[0]): 
        FFT_array = inv_FFT(potential[i,j,:])
        potential[i,j,:] = FFT_array

for i in range(densgrid.shape[0]): 
    for j in range(densgrid.shape[0]): 
        FFT_array = inv_FFT(potential[:,i,j])
        potential[:,i,j] = FFT_array
        
for i in range(densgrid.shape[0]): 
    for j in range(densgrid.shape[0]): 
        FFT_array = inv_FFT(potential[i,:,j])
        potential[i,:,j] = FFT_array
   
potential = -G/np.pi * potential.real

# # TO DO: assign particle masses to densgrid, convert to density, and calculate potentials from it

# # Plotting four slices of a grid

# potential = np.fft.fftn(densgrid)
# potential = -G / np.pi * np.fft.ifftn(potential/k_squared)
# potential = potential.real

fig, ax = plt.subplots(2,2, figsize=(10,8))
pcm = ax[0,0].pcolormesh(np.linspace(0,L, Ngrid), np.linspace(0,L, Ngrid),potential[0,:,:])
ax[0,0].set(ylabel='z (Mpc)', title=r'$\mathrm{x_0}$')
fig.colorbar(pcm, ax=ax[0,0], label='Potential')
pcm = ax[0,1].pcolormesh(np.linspace(0,L, Ngrid), np.linspace(0,L, Ngrid), potential[16,:,:])
ax[0,1].set(ylabel='z (Mpc)', title=r'$\mathrm{x_{16}}$')
fig.colorbar(pcm, ax=ax[0,1], label='Potential')
pcm = ax[1,0].pcolormesh(np.linspace(0,L, Ngrid), np.linspace(0,L, Ngrid), potential[32,:,:])
ax[1,0].set(xlabel='y (Mpc)', title=r'$\mathrm{x_{32}}$')
fig.colorbar(pcm, ax=ax[1,0], label='Potential')
pcm = ax[1,1].pcolormesh(np.linspace(0,L, Ngrid), np.linspace(0,L, Ngrid), potential[64,:,:])
ax[1,1].set(xlabel='y (Mpc)', title=r'$\mathrm{x_{64}}$')
fig.colorbar(pcm, ax=ax[1,1], label='Potential')
ax[0,0].set_aspect('equal', 'box')
ax[0,1].set_aspect('equal', 'box')
ax[1,0].set_aspect('equal', 'box')
ax[1,1].set_aspect('equal', 'box')
# plt.savefig("fig2b.png",dpi=300)
plt.show()
plt.close()
