from astropy.table import Table
import matplotlib.pyplot as plt


# orb_type = 'kepl'
orb_type = 'f'
# orb_type = 'nf'

e = 1.0

orbits_dir = 'orbits_temp_e' + str(e)
if e == 1.0: orbits=16
if e == 0.9: orbits=52


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for i in range(0, orbits):
    orb_file = '../' + orbits_dir + '/orb_%05d.dat' %i
    orb_table = Table.read(orb_file, format='ascii')

    # for orb_type in ['kepl', 'f', 'nf']:
    for orb_type in ['f']:
        x = orb_table['x_'+str(orb_type)]
        y = orb_table['y_'+str(orb_type)]
        z = orb_table['z_'+str(orb_type)]

        ax.plot(x[::100], y[::100], z[::100], label='orb_%05d' %i, linewidth=1.0)

        # if i == orbit:
        ax.scatter(x[0], y[0], z[0], marker='X', color='red', label='Initial point')
        ax.scatter(x[-1], y[-1], z[-1], marker='X', color='purple', label='Final point')


ax.scatter(0, 0, 0, marker='o', color='k')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.legend(loc=0)
plt.show()
