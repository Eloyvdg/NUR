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
    if depth >= max_depth or len(array) == 0: 
        return None
    
    start_index = int(start_index)
    length = int(length)
    
    new_midpoints = np.zeros((3,8))
    start_indices = np.zeros(8)
    start_indices[0] = start_index
    len_array = np.zeros(8)
    
    indices_check = index_array[start_index:start_index+length]
    particles_check = particles[indices_check,:]

    
    size = L/(2**depth)
    offsets = [-size/4, size/4]
    
    ind = 0
    
    new_indices = np.zeros(length)
    for i in range(2):
        for j in range(2): 
            for k in range(2):
                box_center_new = box_center + np.array([offsets[i], offsets[j], offsets[k]])
                new_midpoints[:,ind] = box_center_new
                
                
                box_upper = box_center_new + np.array([offsets[1], offsets[1], offsets[1]])
                box_lower = box_center_new - np.array([offsets[1], offsets[1], offsets[1]])
                
                mask = (
                    (particles_check[:, 0] >= box_lower[0]) & (particles_check[:, 0] <= box_upper[0]) &
                    (particles_check[:, 1] >= box_lower[1]) & (particles_check[:, 1] <= box_upper[1]) &
                    (particles_check[:, 2] >= box_lower[2]) & (particles_check[:, 2] <= box_upper[2]))

                indices = indices_check[mask]  # map relative indices back to global

                len_indices = int(len(indices))
                child_start = int(start_index + np.sum(len_array[:ind]))
                new_indices[child_start:child_start + len_indices] = indices
                
                start_indices[ind] = child_start
                len_array[ind] = len_indices
                                
                ind += 1
                    
    index_array[start_index:start_index+length] = new_indices
    
    for i in range(len(start_indices)-1):
        start_indices[i+1] = int(start_indices[i] + len_array[i])
        
    # print(np.sum(len_array))
        
    
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

half_L = 0.5 * L
center_start = np.array([half_L, half_L, half_L]) 
tree =  build_tree(np.zeros(8), 0, 2, center_start, pos, 0, Np)

first_octant = tree.child1
start_ind = first_octant.start_index
length = first_octant.length
pos_first = pos[index_array[start_ind:start_ind+length],:]

plt.scatter(pos_first[:,0], pos_first[:,1])


# loc = np.where((pos[index_array,0] < 125.0) & (pos[index_array,0] > 0.0) & 
#                 (pos[index_array,1] < 125.0) & (pos[index_array,1] > 0.0) &
#                 (pos[index_array,2] < 125.0) & (pos[index_array,2] > 0.0))

# for level in [3,5,7]: #feel free to change any of this code
#     pixels=2**level
#     massmap=np.zeros(4,(pixels,pixels),dtype=np.float32)
#     # TO DO: traverse the octree, fill map massmap[0,:,:] with the masses of nodes at depth 3 and x_index=x_0,
#     #        massmap[1,:,:] with the masses of nodes at depth 3 and x_index=x_1, etc; then plot these slices;
#     #        then do the same for levels 5 and 7

#     fig, ax = plt.subplots(2,2, figsize=(10,8))
#     pcm = ax[0,0].pcolormesh(np.arange(pixels), np.arange(pixels), massmap[0,:,:])
#     #ax[0,0].set(ylabel='...', title='...')
#     fig.colorbar(pcm, ax=ax[0,0], label='Total mass inside node')
#     pcm =ax[0,1].pcolormesh(np.arange(pixels), np.arange(pixels), massmap[1,:,:])
#     #ax[0,1].set(title='...')
#     fig.colorbar(pcm, ax=ax[0,1], label='Total mass inside node')
#     pcm = ax[1,0].pcolormesh(np.arange(pixels), np.arange(pixels), massmap[2,:,:])
#     #ax[1,0].set(ylabel='...', xlabel='...', title='...')
#     fig.colorbar(pcm, ax=ax[1,0], label='Total mass inside node')
#     pcm = ax[1,1].pcolormesh(np.arange(pixels), np.arange(pixels), massmap[3,:,:])
#     #ax[1,1].set(xlabel='...', title='...')
#     fig.colorbar(pcm, ax=ax[1,1], label='Total mass inside node')
#     ax[0,0].set_aspect('equal', 'box')
#     ax[0,1].set_aspect('equal', 'box')
#     ax[1,0].set_aspect('equal', 'box')
#     ax[1,1].set_aspect('equal', 'box')
#     plt.savefig(f"fig2a_level{level}.png",dpi=300)
#     plt.close()

# # Question 2b: using the FFT

# Ngrid=np.int64(128)
# densgrid=np.zeros((Ngrid,Ngrid,Ngrid),dtype=np.float32)
# potential=np.zeros((Ngrid,Ngrid,Ngrid),dtype=np.float32)
# # TO DO: assign particle masses to densgrid, convert to density, and calculate potentials from it

# # Plotting four slices of a grid

# fig, ax = plt.subplots(2,2, figsize=(10,8))
# pcm = ax[0,0].pcolormesh(np.arange(Ngrid), np.arange(Ngrid), potential[0,:,:])
# #ax[0,0].set(ylabel='...', title='...')
# fig.colorbar(pcm, ax=ax[0,0], label='Potential')
# pcm =ax[0,1].pcolormesh(np.arange(Ngrid), np.arange(Ngrid), potential[16,:,:])
# #ax[0,1].set(title='...')
# fig.colorbar(pcm, ax=ax[0,1], label='Potential')
# pcm = ax[1,0].pcolormesh(np.arange(Ngrid), np.arange(Ngrid), potential[32,:,:])
# #ax[1,0].set(ylabel='...', xlabel='...', title='...')
# fig.colorbar(pcm, ax=ax[1,0], label='Potential')
# pcm = ax[1,1].pcolormesh(np.arange(Ngrid), np.arange(Ngrid), potential[64,:,:])
# #ax[1,1].set(xlabel='...', title='...')
# fig.colorbar(pcm, ax=ax[1,1], label='Potential')
# ax[0,0].set_aspect('equal', 'box')
# ax[0,1].set_aspect('equal', 'box')
# ax[1,0].set_aspect('equal', 'box')
# ax[1,1].set_aspect('equal', 'box')
# plt.savefig("fig2b.png",dpi=300)
# plt.close()
