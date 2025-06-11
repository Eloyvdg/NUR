import numpy as np
import h5py
import matplotlib.pyplot as plt
from Fourier import FFT, inv_FFT

# Question 2: Calculating potentials

#with h5py.File("/disks/cosmodm/DMO_a0.1_256.hdf5","r") as handle:
#    pos=handle["Position"][...] #particle positions, shape (Np,3), comoving
#    #vel=handle["Velocity"][...] #particle velocities, shape (Np,3), comoving <-- not used, but if you're interested

with h5py.File("/net/vdesk/data2/daalen//DMO_a0.1_256.hdf5","r") as handle:
    pos=handle["Position"][...] #particle positions, shape (Np,3), comoving
    #vel=handle["Velocity"][...] #particle velocities, shape (Np,3), comoving <-- not used, but if you're interested

Np=np.int64(256)**3 #number of particles
mp=np.float32(3.64453e10) #particle mass in Msun; all 32-bit to save memory
G=np.float32(4.3009e-9) #gravitational constant in Mpc*(km/s)^2/Msun
h=np.float32(0.3755) #Hubble parameter (this is a Einstein-de Sitter universe with Omega_m=1)
L=np.float32(250.0) #side length of periodic cubic simulation volume
scale_factor=np.float32(0.1) #scale factor a
redshift=1.0/scale_factor-1
rho_mean=Np*mp/L**3 #mean density in Msun/Mpc^3 (comoving, matches 3*H_0^2/(8*pi*G))

# Question 2a: using Barnes-Hut [note: not actually calculating a potential, unless you do the bonus question]
        
class Node(): 
    '''Class that stores all required information of a node in an octree'''
    def __init__(self, depth, child1, child2, child3, child4, child5, child6, child7, child8, 
                 box_center, length, start_index, mass, com):
        self.depth = depth # Depth of the node
        
        # The 8 child nodes of the node
        self.child1 = child1
        self.child2 = child2
        self.child3 = child3
        self.child4 = child4
        self.child5 = child5
        self.child6 = child6
        self.child7 = child7
        self.child8 = child8
        
        self.mass = mass # Total mass in the node
        self.length = length # Length of the index array
        self.start_index = start_index # Start index of the index array
        self.box_center = box_center # Geometrical center of the node 
        self.com = com # Center of mass of the node
        
index_array = np.arange(0,Np) # Global index array with the indices of all the particles
    
def build_tree(array, depth, max_depth, box_center, particles, start_index, length): 
    # Check if the maximum dpeth is reached or if the node is empty (no particles)
    if depth > max_depth or len(array) == 0: 
        return None
    
    
    start_index = int(start_index)
    length = int(length)
    
    start_indices = np.zeros(8) # array for new start indices for 8 child nodes
    start_indices[0] = start_index # For the first child node, it is the same as the previous node
    len_array = np.zeros(8) # Array for new array lengths (number of particles) for 8 child nodes
    octant_indices = [None]*8
    
    # Only select the particles in the specific node
    indices_check = index_array[start_index:start_index+length]
    particles_check = particles[indices_check,:]

    # Check in which child node the particles are located
    mask = particles_check >= box_center
    # Use the bits to give an index to each particle, correpsonding to a child node
    mask = np.where(mask, 1, 0)
    index_octant = (mask[:,0] << 2) | (mask[:,1]) << 1| (mask[:,2] << 0)

    for j in range(8): 
        # Append the indices of the particles to the correct child node
        mask2 = index_octant == j
        new_ind = indices_check[mask2]
        len_array[j] = len(new_ind)
        octant_indices[j] = new_ind
        
    for k in range(1,8): 
        # Define the new start indices for each child node
        start_indices[k] = start_indices[k-1] + len_array[k-1]
        
    for l in range(8): 
        # Sort the global index array again
        index_array[int(start_indices[l]):int(start_indices[l]+len_array[l])] = octant_indices[l]

    
    size = L/(2**depth)
    offsets = [-size/4, size/4]
    new_midpoints = np.zeros((3,8))
    ind = 0

    # Calculate the new geometrical midpoints for each child node
    for i in range(2):
        for j in range(2): 
            for k in range(2):
                box_center_new = box_center + np.array([offsets[i], offsets[j], offsets[k]])
                new_midpoints[:,ind] = box_center_new
                ind += 1
        
    # Calculate the center of mass of the node
    pos_node = pos[index_array[start_index:start_index + length],:]
    com = np.sum(pos_node * mp, axis = 0) / (length * mp)
    
    # Define the new child nodes
    node1 = build_tree(array, depth + 1, max_depth, new_midpoints[:,0], particles, start_indices[0], len_array[0])
    node2 = build_tree(array, depth + 1, max_depth, new_midpoints[:,1], particles, start_indices[1], len_array[1])
    node3 = build_tree(array, depth + 1, max_depth, new_midpoints[:,2], particles, start_indices[2], len_array[2])
    node4 = build_tree(array, depth + 1, max_depth, new_midpoints[:,3], particles, start_indices[3], len_array[3])
    node5 = build_tree(array, depth + 1, max_depth, new_midpoints[:,4], particles, start_indices[4], len_array[4])
    node6 = build_tree(array, depth + 1, max_depth, new_midpoints[:,5], particles, start_indices[5], len_array[5])
    node7 = build_tree(array, depth + 1, max_depth, new_midpoints[:,6], particles, start_indices[6], len_array[6])
    node8 = build_tree(array, depth + 1, max_depth, new_midpoints[:,7], particles, start_indices[7], len_array[7])

    node = Node(depth, node1, node2, node3, node4, node5, node6, node7, node8, box_center, length, start_index, length*mp, com)
    
    return node

def grid(tree, depth): 
    # Function to create a mass grid for the tree for a specific depth
    pixels = 2**depth
    grid = np.zeros((pixels, pixels, pixels))
    center = tree.box_center
    
    x_middle = center[0]/pixels
    def traverse(n): 
        if n == None: 
            return
        
        if n.depth == depth:
            # IF depth is reached, append the mass to the grid at the rught location
            x = n.box_center[0]
            y = n.box_center[1]
            z = n.box_center[2]
            ind_z = int((z+x_middle)/(x_middle*2) -1)
            ind_y = int((y+x_middle)/(x_middle*2) -1)
            ind_x = int((x+x_middle)/(x_middle*2) -1)
            grid[ind_x, ind_y, ind_z] = n.mass
            
            return
        
        if n.depth <= depth:
            # Go one level lower if depth is not reached yet
            for child in [n.child1, n.child2, n.child3, n.child4, n.child5, n.child6, n.child7, n.child8]: 
                traverse(child)
    traverse(tree)
    return grid
                

half_L = 0.5 * L
center_start = np.array([half_L, half_L, half_L]) 
tree =  build_tree(np.zeros(8), 0, 7, center_start, pos, 0, Np)

for level in [3,5,7]: #feel free to change any of this code
    pixels=2**level
    size = L/pixels
    # Calculate the mass map for the depth
    mass_grid = grid(tree, level)
    massmap=np.zeros((4,pixels,pixels),dtype=np.float32)
    # Append the first four slices in the x-direction to the massmap
    for i in range(4): 
        massmap[i,:,:] = mass_grid[i,:,:]

    # Plot the results
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
    plt.savefig(f"fig2a_level{level}.png",dpi=300)
    plt.close()

# Question 2b: using the FFT

Ngrid = 128
densgrid= grid(tree, 7) / ((L/128)**3) # Change the massmap to a density map by dividing by the leaf size
potential= densgrid.copy()
potential = np.complex64(potential)

# Perform a FFT in all three directions, after each other
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
     
# Calculate kx, ky, kz
k = np.zeros(Ngrid)
N_half = Ngrid // 2
k[:N_half] = np.arange(0, N_half)
k[N_half:] = np.arange(-N_half, 0)

d = Ngrid/L
k *= 2 * np.pi / (Ngrid**2)

kx, ky, kz = np.meshgrid(k, k, k)

# Calculate k^2
k_squared = kx**2 + ky**2 + kz**2 

k_squared[0,0,0] = 100 # To avoid division by 0
potential = potential/k_squared # Divide the FFT of the potential by k^2

# Perform an inverse FFT on the result, again along all directions
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
   
potential = -G/np.pi * potential.real # Calculate the potential

# Plot the result at different slices in the x-direction
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
plt.savefig("fig2b.png",dpi=300)
plt.close()
