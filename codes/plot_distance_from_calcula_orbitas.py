from astropy import table
import matplotlib.pyplot as plt
import numpy as np


def dist(orbfile: str, orbtype: str):
    '''
    Calculates distances in 3D

    Parameters
    ----------
    orbfile: str
        file containing the orbits

    orbtype: str
        must be "f", "nf" or "kepl"

    Returns
    -------
        R: Distance between the galaxies in a given orbit in a given time
    '''

    x = orbfile["x_" + orbtype]
    y = orbfile["y_" + orbtype]
    z = orbfile["z_" + orbtype]

    return np.sqrt(x**2 + y**2 + z**2)


x_dt001 = np.linspace(0, 1, 199869)
x_dt01 = np.linspace(0, 1, 20003)
fig, axs = plt.subplots(1, 3, figsize=(30, 10))


for i in range(0, 44, 1):
    orb_dt001 = table.Table.read("../orbits_24th_attempt/orbits_temp_e0.9_dt0.01/orb_%05d.dat" %i, format='ascii', data_start=1)
    orb_dt01 = table.Table.read("../orbits_24th_attempt/orbits_temp_e0.9_dt0.1/orb_%05d.dat" %i, format='ascii', data_start=1)

    print("Calculating distances for f, nf and kepl of dt0.01 orbits. Orbit: {:2d}".format(i))
    dist_f_dt001 = dist(orb_dt001, orbtype="f")
    dist_nf_dt001 = dist(orb_dt001, orbtype="nf")
    dist_kepl_dt001 = dist(orb_dt001, orbtype="kepl")

    print("Calculating distances for f, nf and kepl of dt0.1 orbits. Orbit: {:2d}".format(i)")
    dist_f_dt01 = dist(orb_dt01, orbtype="f")
    dist_nf_dt01 = dist(orb_dt01, orbtype="nf")
    dist_kepl_dt01 = dist(orb_dt01, orbtype="kepl")

    axs[0].plot(x_dt001, dist_f_dt001, color="orange", alpha=0.2)
    axs[0].plot(x_dt01, dist_f_dt01, color="indigo", alpha=0.2)

    axs[1].plot(x_dt001, dist_nf_dt001, color="orange", alpha=0.2)
    axs[1].plot(x_dt01, dist_nf_dt01, color="indigo", alpha=0.2)

    axs[2].plot(x_dt001, dist_kepl_dt001, color="orange", alpha=0.2)
    axs[2].plot(x_dt01, dist_kepl_dt01, color="indigo", alpha=0.2)


axs[0].plot(x_dt001, dist_f_dt001, color="orange", alpha=0.2, label="dt=0.01, friction")
axs[0].plot(x_dt01, dist_f_dt01, color="indigo", alpha=0.2, label="dt=0.1, friction")
axs[1].plot(x_dt001, dist_nf_dt001, color="orange", alpha=0.2, label="dt=0.01, without friction")
axs[1].plot(x_dt01, dist_nf_dt01, color="indigo", alpha=0.2, label="dt=0.1, without friction")
axs[2].plot(x_dt001, dist_kepl_dt001, color="orange", alpha=0.2, label="dt=0.01")
axs[2].plot(x_dt01, dist_kepl_dt01, color="indigo", alpha=0.2, label="dt=0.1")
axs[0].set_title("Orbits with friction")
axs[1].set_title("Orbits without friction")
axs[2].set_title("Keplerian orbits")
axs[0].set_ylabel("Distance between galaxies [kpc]")
axs[1].set_xlabel("Normalized time")
plt.legend()
plt.show()
