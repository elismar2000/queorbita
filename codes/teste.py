import numpy as np
import matplotlib.pyplot as plt

from pygadgetreader import *
from tables import *


snap = 1

snapshot = '/home/elismar/Documentos/Fisica/IC/queorbita/orbits_9th_attempt/bigger_disk/snapshot_0030' #+ str(snap)

mins = '/home/elismar/Documentos/Fisica/IC/queorbita/orbits_9th_attempt/bigger_disk/potential_minima.h5'


def _distances(position, point):
    '''
    Evaluates distances of particles from a given point
    '''
    distances = np.sqrt(np.sum(np.square(position[i] - point[i]) for i in range(3)))
    return distances


##########################################

pos = readsnap(snapshot, 'pos', 'disk')

mass = readsnap(snapshot, 'mass', 'disk')

ids = readsnap(snapshot, 'pid', 'disk')


##########################################

mins_h5file = open_file(mins, 'a')
mins_table = mins_h5file.root.potential.readout

xmin1 = np.array([j['xmin1'] for j in mins_table.iterrows()])
ymin1 = np.array([j['ymin1'] for j in mins_table.iterrows()])
zmin1 = np.array([j['zmin1'] for j in mins_table.iterrows()])
coords_min1 = np.vstack((xmin1, ymin1, zmin1))

xmin2 = np.array([j['xmin2'] for j in mins_table.iterrows()])
ymin2 = np.array([j['ymin2'] for j in mins_table.iterrows()])
zmin2 = np.array([j['zmin2'] for j in mins_table.iterrows()])
coords_min2 = np.vstack((xmin2, ymin2, zmin2))

##########################################

mask2992 = (ids >= 34000) & (ids <= 43999)

mask2993 = (ids >= 44000) & (ids <= 53999)


# dists2992 = np.array([_distances(pos[mask2992][i], coords_min2[:, snap]) for i in range(len(pos[mask2992]))])
# dists2993 = np.array([_distances(pos[mask2993][i], coords_min1[:, snap]) for i in range(len(pos[mask2993]))])

disk_mass = np.sum(mass[mask2992])
print(disk_mass)

plt.hist(mass, bins=50)
plt.show()
