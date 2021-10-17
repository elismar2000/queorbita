#This code can be used to calculate relative velocities and positions
#of galaxies in a Gadget snapshot, which may be useful, for instance to check if the
#initial condition of the simulation was properly assembled (with the correct values for delta_r and delta_v)
#It's possible that it doesn't work properly if the snapshot contains star particles (it is,
#stars formed during the simulation).

#python relative_pos&vel.py <snapshot>


import numpy as np
import matplotlib.pyplot as plt

from pygadgetreader import *
from tables import *

import sys


class Velocity:
    def __init__(self, snapshot):
        self.snapshot = snapshot
        self.ids = np.array([], dtype=int)
        self.vel = np.array([], dtype=float)
        self.pos = np.array([], dtype=float)
        self.pot = np.array([], dtype=float)
        self.mask1 = np.array([], dtype=bool)
        self.mask2 = np.array([], dtype=bool)
        self.coords_min1 = np.array([], dtype=bool)
        self.coords_min2 = np.array([], dtype=bool)
        self.relative_vel = np.array([], dtype=bool)
        self.relative_pos = np.array([], dtype=bool)


    def positions(self):
        pos_halo = readsnap(self.snapshot, 'pos', 'dm')
        pos_disk = readsnap(self.snapshot, 'pos', 'disk')
        pos_gas = readsnap(self.snapshot, 'pos', 'gas')
        pos_bulge = readsnap(self.snapshot, 'pos', 'bulge')

        self.pos = np.concatenate((pos_halo, pos_disk, pos_gas, pos_bulge))


    def velocities(self):
        vel_halo = readsnap(self.snapshot, 'vel', 'dm')
        vel_disk = readsnap(self.snapshot, 'vel', 'disk')
        vel_gas = readsnap(self.snapshot, 'vel', 'gas')
        vel_bulge = readsnap(self.snapshot, 'vel', 'bulge')

        self.vel = np.concatenate((vel_halo, vel_disk, vel_gas, vel_bulge))


    def _ids(self):
        id_halo = readsnap(self.snapshot, 'pid', 'dm')
        id_disk = readsnap(self.snapshot, 'pid', 'disk')
        id_gas = readsnap(self.snapshot, 'pid', 'gas')
        id_bulge = readsnap(self.snapshot, 'pid', 'bulge')

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



    def _distances(self, position, point):
        '''
        Evaluates distances of particles from a given point
        '''
        distances = np.sqrt(np.sum(np.square(position[i] - point[i]) for i in range(3)))
        return distances


    def relative_positions(self):

        self.positions()
        self._ids()
        self.masks()

        pot_halo = readsnap(self.snapshot, 'pot', 'dm')
        pot_disk = readsnap(self.snapshot, 'pot', 'disk')
        pot_gas = readsnap(self.snapshot, 'pot', 'gas')
        pot_bulge = readsnap(self.snapshot, 'pot', 'bulge')

        self.pot = np.concatenate((pot_halo, pot_disk, pot_gas, pot_bulge))

        pot1 = self.pot.min()
        index_min = np.where(self.pot == pot1)[0][0]
        self.coords_min1 = self.pos[index_min]
        print('id da particula com menor potencial = ', self.ids[index_min])
        dist = np.array([self._distances(self.pos[i], self.coords_min1) for i in range(len(self.pos))])

        radius = 10 #10Kpc more or less the radius of the disk of a galaxy
        mask_out = dist > radius

        pot2 = self.pot[mask_out].min()
        index_min2 = np.where(self.pot == pot2)[0][0]
        self.coords_min2 = self.pos[index_min2]

        self.relative_pos = self.coords_min1 - self.coords_min2


    def relative_velocities(self):

        self.velocities()
        self._ids()
        self.masks()

        vel1 = self.vel[self.mask1]
        vel2 = self.vel[self.mask2]

        self.relative_vel = np.mean(vel2, axis=0) - np.mean(vel1, axis=0)



if __name__ == '__main__':

    snapshot = sys.argv[1]

    v = Velocity(snapshot)
    v.relative_positions()
    v.relative_velocities()

    print('Delta_r = {}'.format(v.relative_pos))
    print('Delta_v = {}'.format(v.relative_vel))
