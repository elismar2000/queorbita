import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table as T
from tables import *

#===================================
#Read radial_velocity.h5 files and create a matrix from them
#===================================

path_vel = '/home/elismar/Documentos/Fisica/IC/queorbita/simulations_orbits'
path_params = '/home/elismar/Documentos/Fisica/IC/queorbita/selected_orbits'

parameter = 'q'

list_ij = []

for i in np.arange(1, 9, 1):
    for j in np.arange(0, 5, 1):

        vels = path_vel + str(i) + '/simulation_arp245-' + str(i) + '.' + str(j) + '/radial_velocity.h5'

        with open_file(vels, 'a') as vels_h5file:
            vels_table = vels_h5file.root.velocity.readout
            radial_velocity = np.array([j['radial_vel'] for j in vels_table.iterrows()])

        params = path_params + str(i) + '.txt'
        params_table = T.read(params, format='ascii')

        if (i == 1) & (j == 0):
            radial_vel_matrix = radial_velocity

            params_array = [params_table[parameter][j]]

        else:
            radial_vel_matrix = np.vstack((radial_vel_matrix, radial_velocity))

            params_array.append(params_table[parameter][j])

        list_ij.append(f'{i}, {j}')


#===================================
#Plot results
#===================================

def highlight_cell(x,y, ax=None, **kwargs):
    rect = plt.Rectangle((x-.5, y-.5), 1,1, fill=False, **kwargs)
    ax = ax or plt.gca()
    ax.add_patch(rect)
    return rect

sorted_idx = np.argsort(params_array)
radial_vel_matrix = radial_vel_matrix[sorted_idx]

indices = np.argwhere((radial_vel_matrix < -119+34) & (radial_vel_matrix > -119-34))

fig, ax = plt.subplots()

im = ax.imshow(radial_vel_matrix, cmap='seismic')

for i in range(indices.shape[0]):
     highlight_cell(indices[i][1], indices[i][0], color="limegreen", linewidth=2)

labels = [f"{i} -- {j}" for i, j in zip(np.array(list_ij)[sorted_idx], np.array(params_array)[sorted_idx])]

ax.set_yticks(np.arange(0, 40))
ax.set_yticklabels(labels)

ax.tick_params(labelsize=8, axis='y')

ax.set_ylabel(parameter, fontsize=10)
ax.set_xlabel('Snapshot', fontsize=10)

cbar = plt.colorbar(mappable=im, ax=ax)
cbar.ax.set_ylabel('Radial velocity (km/s)', rotation=90, fontsize=10)

plt.show()
