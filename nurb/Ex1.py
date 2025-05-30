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

G = c.G.to(u.AU**3 * u.d**-2 * u.kg**-1)

def acc_tensor(position, masses): 
    x, y, z = position[:,0], position[:,1], position[:,2]
    tensor = np.zeros((x.shape[0],x.shape[0],3))
    for i in range(tensor.shape[0]):
        for j in range(i + 1, tensor.shape[0]):
            
            pos_i = np.array([x[i], y[i], z[i]]) 
            pos_j = np.array([x[j], y[j], z[j]]) 
            
            acc = specific_accelaration(pos_i, pos_j) 
            
            tensor[i,j,:] = -acc * masses[j]
            tensor[j,i,:] = acc  * masses[i]

        
    # tensor_sum = np.sum(tensor, axis = 1)
    # for i, mass in enumerate(masses): 
    #     tensor_sum[i,:] = tensor_sum[i,:] * mass

    # print(tensor)
    return np.sum(tensor, axis = 1)

def specific_accelaration(r1, r2): 
    return -G.value / (np.linalg.norm(r2 - r1)**3) * (r2 - r1)

# pick a time (please use either this or the current time)
t = Time("2021-12-07 10:00")
delta_t = 12 * 3600

# initialize the planets; Mars is shown as an example


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
    
masses = np.array([0.0553, 0.815, 1, 0.1075, 317.8, 95.2, 14.6, 17.2]) * 5.9722e24
masses = np.append(np.array([1.989e30]), masses)

# calculate the x position in AU
# print(marsposition.x.to_value(u.AU))

# calculate the v_x velocity in AU/day
# print(marsvelocity.x.to_value(u.AU/u.d))

# Problem 1.a
positions = [sun[0], mercury[0], venus[0], earth[0], mars[0], jupiter[0], saturn[0], uranus[0], neptune[0]]
velocities = [sun[1], mercury[1], venus[1], earth[1], mars[1], jupiter[1], saturn[1], uranus[1], neptune[1]]

pos = np.zeros((9,3))
vel = np.zeros((9,3))

for i, obj in enumerate(positions): 
    pos[i, 0] = obj.x.to_value(u.AU)
    pos[i, 1] = obj.y.to_value(u.AU)
    pos[i, 2] = obj.z.to_value(u.AU)
for i, obj in enumerate(velocities):
    vel[i, 0] = obj.x.to_value(u.AU/u.d)
    vel[i, 1] = obj.y.to_value(u.AU/u.d)
    vel[i, 2] = obj.z.to_value(u.AU/u.d)
 
x, y, z = pos[:,0], pos[:,1], pos[:,2]

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
# plt.savefig("fig1a.png")
# plt.close()
plt.show()

# Problem 1.b
# For visibility, you may want to do two versions of this plot: 
# one with all planets, and another zoomed in on the four inner planets

time = np.linspace(0, 200*365, 200*365*2 + 1)
delta_t = np.diff(time)[0]

positions = [sun[0], mercury[0], venus[0], earth[0], mars[0], jupiter[0], saturn[0], uranus[0], neptune[0]]
velocities = [sun[1], mercury[1], venus[1], earth[1], mars[1], jupiter[1], saturn[1], uranus[1], neptune[1]]

pos_l = np.zeros((9,len(time),3))
vel_l = np.zeros((9,len(time),3))

for i, obj in enumerate(positions): 
    pos_l[i, 0, 0] = obj.x.to_value(u.AU)
    pos_l[i, 0, 1] = obj.y.to_value(u.AU)
    pos_l[i, 0, 2] = obj.z.to_value(u.AU)
for i, obj in enumerate(velocities):
    vel_l[i, 0, 0] = obj.x.to_value(u.AU/u.d)
    vel_l[i, 0, 1] = obj.y.to_value(u.AU/u.d)
    vel_l[i, 0, 2] = obj.z.to_value(u.AU/u.d)
    
vel_half = vel_l[:,0,:] + 0.5 * acc_tensor(pos_l[:,0,:], masses) * delta_t

for t in range(len(time) - 1): 
    pos_l[:,t+1,:] = pos_l[:,t,:] + vel_half * delta_t
    
    acc = acc_tensor(pos_l[:,t+1,:], masses)
    vel_half = vel_half + acc * delta_t
    
    vel_l[:,t+1,:] = vel_l[:,t,:] + 0.5 * vel_half * delta_t

    
x_l, y_l, z_l = pos_l[:,:,0], pos_l[:,:,1], pos_l[:,:,2]
time = x_l.copy()*0 +np.linspace(0,200 * 365, 200*365*2 +1)

fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(x_l[i,:], y_l[i,:], label=obj)
    ax[1].plot(time[i,:]*u.d.to(u.yr), z_l[i,:], label=obj)
ax[0].set_aspect('equal', 'box')
ax[0].set(xlabel='X [AU]', ylabel='Y [AU]')
ax[1].set(xlabel='Time [yr]', ylabel='Z [AU]')
plt.legend(loc=(1.05,0))
# plt.savefig("fig1b.png")
# plt.close()
plt.show()


# Problem 1.c
time = np.linspace(0, 200*365, 200*365*2 + 1)
delta_t = np.diff(time)[0]

positions = [sun[0], mercury[0], venus[0], earth[0], mars[0], jupiter[0], saturn[0], uranus[0], neptune[0]]
velocities = [sun[1], mercury[1], venus[1], earth[1], mars[1], jupiter[1], saturn[1], uranus[1], neptune[1]]

pos_r = np.zeros((9,len(time),3))
vel_r = np.zeros((9,len(time),3))

for i, obj in enumerate(positions): 
    pos_r[i, 0, 0] = obj.x.to_value(u.AU)
    pos_r[i, 0, 1] = obj.y.to_value(u.AU)
    pos_r[i, 0, 2] = obj.z.to_value(u.AU)
for i, obj in enumerate(velocities):
    vel_r[i, 0, 0] = obj.x.to_value(u.AU/u.d)
    vel_r[i, 0, 1] = obj.y.to_value(u.AU/u.d)
    vel_r[i, 0, 2] = obj.z.to_value(u.AU/u.d)


for t in range(len(time) -1):
    vel_start = vel_r[:,t,:]
    acc_start = acc_tensor(pos_r[:,t,:], masses)
    
    k1_vel = acc_tensor(pos_r[:,t,:], masses) * delta_t
    k1_pos = vel_start * delta_t
    
    k2_vel = acc_tensor(pos_r[:,t,:] + 0.5 * k1_pos, masses) * delta_t
    k2_pos = (vel_start + 0.5 * k1_vel)* delta_t
    
    k3_vel = acc_tensor(pos_r[:,t,:] + 0.5 * k2_pos, masses) * delta_t
    k3_pos = (vel_start + 0.5 * k2_vel) * delta_t
    
    k4_vel = acc_tensor(pos_r[:,t,:] + k3_pos, masses) * delta_t
    k4_pos = (vel_start + k3_vel) * delta_t
    
    vel_r[:,t+1,:] = vel_start + 1/6 * (k1_vel + k4_vel + 2 * (k2_vel + k3_vel))
    pos_r[:,t+1,:] = pos_r[:,t,:] + 1/6 * (k1_pos + k4_pos + 2 * (k2_pos + k3_pos))

x_r, y_r, z_r = pos_r[:,:,0], pos_r[:,:,1], pos_r[:,:,2]
time = x_r.copy()*0 +np.linspace(0,200 * 365, 200*365*2 +1)

fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(x_r[i,:], y_r[i,:], label=obj)
    ax[1].plot(time[i,:]*u.d.to(u.yr), z_r[i,:], label=obj)
ax[0].set_aspect('equal', 'box')
ax[0].set(xlabel='X [AU]', ylabel='Y [AU]')
ax[1].set(xlabel='Time [yr]', ylabel='Z [AU]')
plt.legend(loc=(1.05,0))
# plt.savefig("fig1b.png")
# plt.close()
plt.show()

fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(x_l[i,:], y_l[i,:], label=obj)
    ax[1].plot(x_r[i,:], y_r[i,:], label=obj)
ax[0].set(xlabel='X [AU]', ylabel='Y [AU]', title='Leapfrog')
ax[1].set(xlabel='X [AU]', ylabel='Y [AU]', title='RK4')
plt.legend(loc=(1.05,0))
plt.show()

fig, ax = plt.subplots(1,3, figsize=(18,5), constrained_layout=True)
for i, obj in enumerate(names):
    ax[0].plot(time[i,:]*u.d.to(u.yr), z_l[i,:], label=obj)
    ax[1].plot(time[i,:]*u.d.to(u.yr), z_r[i,:], label=obj)
    ax[2].plot(time[i,:]*u.d.to(u.yr), np.abs(x_l[i,:] - x_r[i,:]), label=obj)
ax[0].set(xlabel='Time [yr]', ylabel='Z [AU]', title='Leapfrog')
ax[1].set(xlabel='Time [yr]', ylabel='Z [AU]', title='RK4')
ax[2].set(xlabel='Time [yr]', ylabel='X [AU]', title='Absolute difference in x-position')
plt.legend(loc=(1.05,0))
plt.show()

# fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
# for i, obj in enumerate(names):
#     ax[0].plot(time[i,:]*u.d.to(u.yr), z_l[i,:], label=obj)
#     ax[1].plot(time[i,:]*u.d.to(u.yr), z_r[i,:], label=obj)
# ax[0].set(xlabel='Time [yr]', ylabel='Z [AU]', title='Leapfrog')
# ax[1].set(xlabel='Time [yr]', ylabel='Z [AU]', title='RK4')
# plt.legend(loc=(1.05,0))
# plt.show()

# fig, ax = plt.subplots(1,2, figsize=(12,5), constrained_layout=True)
# for i, obj in enumerate(names):
#     ax[0].scatter(x[i,:], y_r[i,:], label=obj)
#     ax[1].scatter(x[i,:], z_r[i,:], label=obj)
# ax[0].set_aspect('equal', 'box')
# ax[1].set_aspect('equal', 'box')
# ax[0].set(xlabel='X [AU]', ylabel='Y [AU]')
# ax[1].set(xlabel='X [AU]', ylabel='Z [AU]')
# plt.legend(loc=(1.05,0))
# # plt.savefig("fig1a.png")
# # plt.close()
# plt.show()



# plt.savefig("fig1c.png")
# plt.close()
