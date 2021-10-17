from astropy.table import Table

import os
import subprocess
import sys

#============================
#Reading input
#============================

orbits_group = sys.argv[1] #the number of the group of orbits to be tested

orbits = 'selected_orbits' + orbits_group

#============================
#Checking directories and running calcula_orbitas
#============================

cwd = os.getcwd()
print("You start from the directory:", cwd)
print("You have to be in the same directory as calcula_orbitas_queorbita.py\n")

print("Running calcula_orbitas_queorbita.py")
subprocess.Popen(['python', 'calcula_orbitas_queorbita.py', 'orbits'])

#============================
#Renaming orbits_temp directory
#============================

orbits_temp_path = 'orbits_temp' + orbits_group

os.rename('orbits_temp', orbits_temp_path)

print('')

#============================
#Creating directories for the simulations
#============================

os.mkdir('simulation_orbits' + orbits_group)
os.chdir('simulation_orbits' + orbits_group)

for i in range(5):
    os.mkdir('simulation_arp245-' + orbits_group + '.' + i + '/')

print('Joining galaxies\n')

#============================
#Running join_galaxies
#============================

os.chdir('/home/elismar/Documentos/Fisica/IC/queorbita/' + orbits_temp_path)

ngc2992_ic = '/home/elismar/Documentos/Fisica/IC/simulations_ICs/galstep/galstep/ngc2992_rotated.ic'
ngc2993_ic = '/home/elismar/Documentos/Fisica/IC/simulations_ICs/galstep/galstep/ngc2993_rotated.ic'

for orbit in range(5):
    orb = 'orb_0000' + str(orbit) + '.dat'
    t = Table.read(orb, format='ascii')

    print(f'Joining orbit {orbit}\n')

    x_kepl = t['x_kepl'][0]
    y_kepl = t['y_kepl'][0]
    z_kepl = t['z_kepl'][0]
    Vx_kepl = t['Vx_kepl'][0]
    Vy_kepl = t['Vy_kepl'][0]
    Vz_kepl = t['Vz_kepl'][0]

    output = #colocar aqui o diretório do output para a condição inicial

    subprocess.Popen(['python', '/home/elismar/Documentos/Fisica/IC/simulations_ICs/codes/join_galaxies2.py',
                    'x_kepl', 'y_kepl', 'z_kepl', 'Vx_kepl', 'Vy_kepl', 'Vz_kepl'])
