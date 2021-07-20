import numpy as np
import matplotlib.pyplot as plt

from pygadgetreader import *
from tables import *

import sys


class Velocity:
    def __init__(self, path):
        self.path = path
        self.snap_path = str()
        self.ids = np.array([], dtype=int)
        self.vel = np.array([], dtype=float)
        self.mask1 = np.array([], dtype=bool)
        self.mask2 = np.array([], dtype=bool)
        self.radial_vel_i = float()


    def velocities(self):
        vel_halo = readsnap(self.snap_path, 'vel', 'dm')
        vel_disk = readsnap(self.snap_path, 'vel', 'disk')
        vel_gas = readsnap(self.snap_path, 'vel', 'gas')
        vel_bulge = readsnap(self.snap_path, 'vel', 'bulge')

        self.vel = np.concatenate((vel_halo, vel_disk, vel_gas, vel_bulge))


    def _ids(self):
        id_halo = readsnap(self.snap_path, 'pid', 'dm')
        id_disk = readsnap(self.snap_path, 'pid', 'disk')
        id_gas = readsnap(self.snap_path, 'pid', 'gas')
        id_bulge = readsnap(self.snap_path, 'pid', 'bulge')

        self.ids = np.concatenate((id_halo, id_disk, id_gas, id_bulge))


    def masks(self):
        #For 120.000 particles:
        #ids galaxy 1: 0 - 9999; 20000 - 39999; 60000 - 83999; 108000 - 113999
        #ids galaxy 2: 10000 - 19999; 40000 - 59999; 84000 - 107999; 114000 - 119999

        #for 60.000 particles:
        self.mask1 = (self.ids <= 4999) | (self.ids >= 10000) & (self.ids <= 19999) | (self.ids >= 30000) \
                & (self.ids <= 41999) | (self.ids >= 54000) & (self.ids <= 56999)


        self.mask2 = (self.ids >= 5000) & (self.ids <= 9999) | (self.ids >= 20000) & (self.ids <= 29999) \
                | (self.ids >= 42000) & (self.ids <= 53999) | (self.ids >= 57000) & (self.ids <= 59999)


    def mean_vel(self):

        self.velocities()
        self._ids()
        self.masks()

        vel1 = self.vel[self.mask1]
        vel2 = self.vel[self.mask2]

        self.radial_vel_i = np.mean(vel2[:, 2]) - np.mean(vel1[:, 2])


    def save_table(self, num_of_snaps, step):
        '''
        perform the iteration over all snapshots and
        save table with the values for radial velocity

        Parameters
        ----------
        num_of_snaps : int
            The significant figure of the last snapshot of the simulation

        step : int
            Optional step value in case one wants to iterate
            just over an even spaced subset of the snapshots
        '''

        class Table(IsDescription):
            snapshot = Int32Col()
            radial_vel = Float32Col()

        h5file = open_file(self.path + 'radial_velocity.h5', mode='w', title='Radial Velocity')
        group = h5file.create_group("/", 'velocity', 'Radial velocity of the galaxies')
        table = h5file.create_table(group, 'readout', Table, 'Readout table')
        vel = table.row

        for i in range(0, num_of_snaps+1, step):
            snapshot = 'snapshot_' + '%04d' % (i,)
            self.snap_path = self.path + snapshot
            self.mean_vel()

            vel['snapshot'] = i
            vel['radial_vel'] = self.radial_vel_i
            vel.append()


        table.flush()


if __name__ == '__main__':

    path         = sys.argv[1]
    num_of_snaps = sys.argv[2]
    step         = sys.argv[3]

    v = Velocity(path)
    v.save_table(num_of_snaps=int(num_of_snaps), step=int(step))
