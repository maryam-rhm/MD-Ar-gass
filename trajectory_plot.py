import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import seaborn as sns

data = pd.read_csv('result.csv')
steps = sorted(data['Step'].unique())

x_positions = [data[data['Step'] == step][' x'].values for step in steps]
y_positions = [data[data['Step'] == step][' y'].values for step in steps]

fig, ax = plt.subplots()
scat = ax.scatter(x_positions[0], y_positions[0])

ax.set_xlim(0, 20)
ax.set_ylim(0, 20)
ax.set_xlabel('X Position (angstrom)')
ax.set_ylabel('Y Position (angstrom)')
ax.set_title('Particle Trajectories Over Time')

def update(frame):
    scat.set_offsets(list(zip(x_positions[frame], y_positions[frame])))
    return scat,
    
ani = animation.FuncAnimation(fig, update, frames=len(steps), interval=50, blit=True)
ani.save('t_animation.mp4', writer='ffmpeg')
plt.show()

#### Heatmap ####
heatmap, xedges, yedges = np.histogram2d(data[' x'], data[' y'], bins=(100, 100), range=[[0, 20], [0, 20]])

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(heatmap.T, cmap="viridis", ax=ax)
ax.set_xlabel('X Position (A)')
ax.set_ylabel('Y Position (A)')
ax.set_title('Heatmap of Particles Position')
ax.set_xticks([])
ax.set_yticks([])

plt.savefig('p_heatmap.png')
plt.show()

### Plot Energies over Time ###
pot = sorted(data[' Potential Energy(kJ/mol)'].unique())
kin = sorted(data[' Kinetic Energy(kJ/mol)'].unique())
tot = sorted(data[' Total Energy(kJ/mol)'].unique())
temp = sorted(data[' Temperture(K)'].unique())

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
ax.plot(steps, pot, label='Potential Energy')
ax.set_xlabel('Step')
ax.set_ylabel('Energy (kJ/mol)')
ax.set_title('Potential Energy vs Time')
plt.savefig('Potential Energy.png')
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
ax.plot(steps, kin, label='Kinetic Energy')
ax.set_xlabel('Step')
ax.set_ylabel('Energy (kJ/mol)')
ax.set_title('Kinetic Energy vs Time')
plt.savefig('Kinetic Energy.png')
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
ax.plot(steps, tot, label='Total Energy')
ax.set_xlabel('Step')
ax.set_ylabel('Energy (kJ/mol)')
ax.set_title('Total Energy vs Time')
plt.savefig('Total Energy.png')
plt.show()

fig, ax = plt.subplots(1, 1, figsize=(5, 4))
ax.plot(steps, temp, label='Temperature')
ax.set_xlabel('Step')
ax.set_ylabel('Temperature (K)')
ax.set_title('Instantaneous Temperature  vs Time')
plt.savefig('Temperature.png')
plt.show()
