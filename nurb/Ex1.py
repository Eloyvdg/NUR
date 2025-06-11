# we need to import the time module from astropy
from astropy.time import Time
# import some coordinate things from astropy
from astropy.coordinates import solar_system_ephemeris
from astropy.coordinates import get_body_barycentric_posvel
from astropy import units as u
from astropy import constants as c
import numpy as np
import matplotlib.pyplot as plt

# Question 1: Simulating the solar system

# Change G to correct units
G = c.G.to(u.AU**3 * u.d**-2 * u.kg**-1)
  
def specific_acceleration(r1, r2): 
    # Calculates the specific acceleration 
    return -G.value / ( (np.sqrt(np.sum((r2 - r1)**2)) ) **3) * (r2 - r1)
 
def acc_tensor(position, masses): 
    # Calculates the acceleration tensor
    x, y, z = position[:,0], position[:,1], position[:,2] # Define x, y, z
    tensor = np.zeros((x.shape[0],x.shape[0],3))
    for i in range(tensor.shape[0]):
        for j in range(i + 1, tensor.shape[0]):
            
            pos_i = np.array([x[i], y[i], z[i]]) 
            pos_j = np.array([x[j], y[j], z[j]]) 
            
            # Calculate the specific acceleration for the positions of 2 objects
            acc = specific_acceleration(pos_i, pos_j) 
            
            # Update the tensor
            tensor[i,j,:] = -acc * masses[j]
            tensor[j,i,:] = acc  * masses[i]
    return np.sum(tensor, axis = 1)

# pick a time (please use either this or the current time)
t = Time("2021-12-07 10:00")
delta_t = 12 * 3600 # Half a day

with solar_system_ephemeris.set('jpl'):
    sun = get_body_barycentric_posvel('sun', t)
    mercury = get_body_barycentric_posvel('mercury', t)
    venus = get_body_barycentric_posvel('venus', t)
    earth = get_body_barycentric_posvel('earth', t)
    mars = get_body_barycentric_posvel('mars', t)
    jupiter = get_body_barycentric_posvel('jupiter', t)
    saturn = get_body_barycentric_posvel('saturn', t)
    uranus = get_body_barycentric_posvel('uranus', t)
    neptune = get_body_barycentric_posvel('neptune', t)
    
masses = np.array([0.0553, 0.815, 1, 0.1075, 317.8, 95.2, 14.6, 17.2]) * 5.9722e24 # masses in kg
masses = np.append(np.array([1.989e30]), masses) # Solar mass

# Problem 1.a
positions = [sun[0], mercury[0], venus[0], earth[0], mars[0], jupiter[0], saturn[0], uranus[0], neptune[0]]
velocities = [sun[1], mercury[1], venus[1], earth[1], mars[1], jupiter[1], saturn[1], uranus[1], neptune[1]]

pos = np.zeros((9,3))
vel = np.zeros((9,3))

# Get the positions and velocities in correct units
for i, obj in enumerate(positions): 
    pos[i, 0] = obj.x.to_value(u.AU)
    pos[i, 1] = obj.y.to_value(u.AU)
    pos[i, 2] = obj.z.to_value(u.AU)
for i, obj in enumerate(velocities):
    vel[i, 0] = obj.x.to_value(u.AU/u.d)
    vel[i, 1] = obj.y.to_value(u.AU/u.d)
    vel[i, 2] = obj.z.to_value(u.AU/u.d)
 
x, y, z = pos[:,0], pos[:,1], pos[:,2]

# Plot the initial positions of the the planets and the Sun
names = np.array(['Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune'])
fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].scatter(x[i], y[i], label=obj)
    ax[1].scatter(x[i], z[i], label=obj)
ax[0].set_aspect('equal', 'box')
ax[1].set_aspect('equal', 'box')
ax[0].set(xlabel='X [AU]', ylabel='Y [AU]')
ax[1].set(xlabel='X [AU]', ylabel='Z [AU]')
plt.legend(loc=(1.05,0))
plt.savefig("fig1a.png", dpi = 300)
plt.close()

# Problem 1.b
time = np.linspace(0, 200*365, 200*365*2 + 1)
delta_t = np.diff(time)[0]

positions = [sun[0], mercury[0], venus[0], earth[0], mars[0], jupiter[0], saturn[0], uranus[0], neptune[0]]
velocities = [sun[1], mercury[1], venus[1], earth[1], mars[1], jupiter[1], saturn[1], uranus[1], neptune[1]]

pos_l = np.zeros((9,len(time),3))
vel_l = np.zeros((9,len(time),3))

# Get the initial positions and velocities
for i, obj in enumerate(positions): 
    pos_l[i, 0, 0] = obj.x.to_value(u.AU)
    pos_l[i, 0, 1] = obj.y.to_value(u.AU)
    pos_l[i, 0, 2] = obj.z.to_value(u.AU)
for i, obj in enumerate(velocities):
    vel_l[i, 0, 0] = obj.x.to_value(u.AU/u.d)
    vel_l[i, 0, 1] = obj.y.to_value(u.AU/u.d)
    vel_l[i, 0, 2] = obj.z.to_value(u.AU/u.d)
    

# Start of the leapfrog algorithm
# Calculate the first half step of the algorithm
vel_half = vel_l[:,0,:] + 0.5 * acc_tensor(pos_l[:,0,:], masses) * delta_t

for t in range(len(time) - 1): 
    # Calculate the new position
    pos_l[:,t+1,:] = pos_l[:,t,:] + vel_half * delta_t
    
    # Calculate the acceleration tensor corresponding to this new position
    acc = acc_tensor(pos_l[:,t+1,:], masses)
    
    # Update the half velocity step and calculate the velocity at the new position
    vel_l[:,t+1,:] = vel_half + 0.5 * acc * delta_t
    vel_half = vel_half + acc * delta_t


    
x_l, y_l, z_l = pos_l[:,:,0], pos_l[:,:,1], pos_l[:,:,2]
time = x_l.copy()*0 +np.linspace(0,200 * 365, 200*365*2 +1)

# Plot the results of the leapfrog algorithm
fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(x_l[i,:], y_l[i,:], label=obj)
    ax[1].plot(time[i,:]*u.d.to(u.yr), z_l[i,:], label=obj)
ax[0].set_aspect('equal', 'box')
ax[0].set(xlabel='X [AU]', ylabel='Y [AU]')
ax[1].set(xlabel='Time [yr]', ylabel='Z [AU]')
plt.legend(loc=(1.05,0))
plt.savefig("fig1b.png", dpi = 300)
plt.close()


# Problem 1.c
time = np.linspace(0, 200*365, 200*365*2 + 1)
delta_t = np.diff(time)[0]

positions = [sun[0], mercury[0], venus[0], earth[0], mars[0], jupiter[0], saturn[0], uranus[0], neptune[0]]
velocities = [sun[1], mercury[1], venus[1], earth[1], mars[1], jupiter[1], saturn[1], uranus[1], neptune[1]]

pos_r = np.zeros((9,len(time),3))
vel_r = np.zeros((9,len(time),3))

# Get the initial positions and velocities
for i, obj in enumerate(positions): 
    pos_r[i, 0, 0] = obj.x.to_value(u.AU)
    pos_r[i, 0, 1] = obj.y.to_value(u.AU)
    pos_r[i, 0, 2] = obj.z.to_value(u.AU)
for i, obj in enumerate(velocities):
    vel_r[i, 0, 0] = obj.x.to_value(u.AU/u.d)
    vel_r[i, 0, 1] = obj.y.to_value(u.AU/u.d)
    vel_r[i, 0, 2] = obj.z.to_value(u.AU/u.d)


# RK4 algorithm
for t in range(len(time) -1):
    vel_start = vel_r[:,t,:]
    acc_start = acc_tensor(pos_r[:,t,:], masses)
    
    # Calculate k1
    k1_vel = acc_tensor(pos_r[:,t,:], masses) * delta_t
    k1_pos = vel_start * delta_t
    
    # Calculate k2
    k2_vel = acc_tensor(pos_r[:,t,:] + 0.5 * k1_pos, masses) * delta_t
    k2_pos = (vel_start + 0.5 * k1_vel)* delta_t
    
    # Calculate k3
    k3_vel = acc_tensor(pos_r[:,t,:] + 0.5 * k2_pos, masses) * delta_t
    k3_pos = (vel_start + 0.5 * k2_vel) * delta_t
    
    # Calculate k4
    k4_vel = acc_tensor(pos_r[:,t,:] + k3_pos, masses) * delta_t
    k4_pos = (vel_start + k3_vel) * delta_t
    
    # Calculate the new positions and velocities with a weighted average
    vel_r[:,t+1,:] = vel_start + 1/6 * (k1_vel + k4_vel + 2 * (k2_vel + k3_vel))
    pos_r[:,t+1,:] = pos_r[:,t,:] + 1/6 * (k1_pos + k4_pos + 2 * (k2_pos + k3_pos))


x_r, y_r, z_r = pos_r[:,:,0], pos_r[:,:,1], pos_r[:,:,2]
time = x_r.copy()*0 +np.linspace(0,200 * 365, 200*365*2 +1)

# Plot the results
fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(x_r[i,:], y_r[i,:], label=obj)
    ax[1].plot(time[i,:]*u.d.to(u.yr), z_r[i,:], label=obj)
ax[0].set_aspect('equal', 'box')
ax[0].set(xlabel='X [AU]', ylabel='Y [AU]')
ax[1].set(xlabel='Time [yr]', ylabel='Z [AU]')
plt.legend(loc=(1.05,0))
plt.savefig("fig1c.png")
plt.close()

# Plot both the orbits of leapfrog and RK4
fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(x_l[i,:], y_l[i,:], label=obj)
    ax[1].plot(x_r[i,:], y_r[i,:], label=obj)
ax[0].set(xlabel='X [AU]', ylabel='Y [AU]', title='Leapfrog')
ax[1].set(xlabel='X [AU]', ylabel='Y [AU]', title='RK4')
plt.legend(loc=(1.05,0))
plt.savefig("fig1c1.png")
plt.close()

# Plot the z positions for both leapfrog and RK4, and the absolute difference in x-positions between the two algorithm
fig, ax = plt.subplots(1,3, figsize=(18,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(time[i,:]*u.d.to(u.yr), z_l[i,:], label=obj)
    ax[1].plot(time[i,:]*u.d.to(u.yr), z_r[i,:], label=obj)
    ax[2].plot(time[i,:]*u.d.to(u.yr), np.abs(x_l[i,:] - x_r[i,:]), label=obj)
ax[0].set(xlabel='Time [yr]', ylabel='Z [AU]', title='Leapfrog')
ax[1].set(xlabel='Time [yr]', ylabel='Z [AU]', title='RK4')
ax[2].set(xlabel='Time [yr]', ylabel='X [AU]', title='Absolute difference in x-position')
plt.legend(loc=(1.05,0))
plt.savefig("fig1c2.png")
plt.close()

# Plot both the orbits of leapfrog and RK4
fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(x_l[i,:], y_l[i,:], label=obj)
    ax[1].plot(x_r[i,:], y_r[i,:], label=obj)
ax[0].set(xlabel='X [AU]', ylabel='Y [AU]', title='Leapfrog')
ax[1].set(xlabel='X [AU]', ylabel='Y [AU]', title='RK4')
ax[0].set_xlim(-1.8, 1.8)
ax[0].set_ylim(-1.8, 1.8)
ax[1].set_xlim(-1.8, 1.8)
ax[1].set_ylim(-1.8, 1.8)
plt.legend(loc=(1.05,0))
plt.savefig("fig1c3.png")
plt.close()
