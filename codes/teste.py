from pygadgetreader import *

snapshot = '/home/elismar/Documentos/Fisica/IC/queorbita/orbits_3rd_attempt/orb39_e0.9-pot/snapshot_0000'
vel = readsnap(snapshot, 'vel')

print(vel)
