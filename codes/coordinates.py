import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord


c1x = SkyCoord(ra=146.42*u.deg, dec=-14.32*u.deg)
c2x = SkyCoord(ra=146.42*u.deg, dec=-14.36*u.deg)

dx = c1x.separation(c2x)

#====================================

c1y = SkyCoord(ra=146.42*u.deg, dec=-14.36*u.deg)
c2y = SkyCoord(ra=146.45*u.deg, dec=-14.36*u.deg)

dy = c1x.separation(c2y)

#====================================

distance = 38e+3 * u.kpc

dx = np.deg2rad(dx) * distance
dy = np.deg2rad(dy) * distance
