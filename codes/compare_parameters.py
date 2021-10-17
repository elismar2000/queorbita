from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

import glob

# sim1 = '8.3'
# sim2 = '4.3'

v = glob.glob('../selected_orbits*.txt')

v.sort()

for i in v:
    t = Table.read(i, format='ascii')
    for sim in range(0,1):

    # t2 = Table.read('../selected_orbits' + str(sim2[0]) + '.txt', format='ascii')

        params = np.delete(np.vstack(np.array(t[sim][str(cols)]) for cols in t.columns).T, [21, 23])
        # params2 = np.delete(np.vstack(np.array(t2[int(sim2[2])][str(cols)]) for cols in t2.columns).T, [21, 23])

        params = np.array([float(params[i]) for i in range(0, len(params))])
        # params2 = np.array([float(params2[i]) for i in range(0, len(params2))])

        columns = np.delete(np.array([cols for cols in t.columns]), [21, 23])

        plt.plot(columns, params, label=i[-5]+str(sim))
        # plt.plot(columns, params2, color='gold', label=sim2)


plt.xticks(rotation=45)
# plt.title(sim1 + ' & ' + sim2)
plt.legend()
plt.show()
