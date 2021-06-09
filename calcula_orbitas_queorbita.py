# Programa "calcula_orbitas_queorbita.py#
# Lê arquivo de saída do "queorbita"
# Pra carda órbita, roda o "pot_hbd_din_fric-2" pra frente e pra trás no tempo
#    por 1 milhão de anos. Os arquivos de saída ficam no diretório "orbits_temp"


import numpy as np
from astropy import table
import pandas as pd                 # Para manipular arquivos: read_csv
import plotly.graph_objects as go   # Para plot em 3d
import subprocess                   # Para rodar programas externos
import os                           # Para rodar programas externos



#queorb = pd.read_csv('queorbita_new.out' ,  header=50, sep=",")
queorb = table.Table.read('selected_orbits2.txt', format='ascii')
queorb = queorb.to_pandas()

#queorb.dtypes

#queorb.head()

#queorb.tail()

#print(queorb)

e           = queorb['e']
q           = queorb['q']
rnow        = queorb['rnow']
vq          = queorb['vq']
gx          = queorb['gx']
gy          = queorb['gy']
gz          = queorb['gz']
VX_ceu      = queorb['VX_ceu']
VY_ceu      = queorb['VY_ceu']
VZ_ceu      = queorb['VZ_ceu']
x           = queorb['x']
y           = queorb['y']
vx          = queorb['vx']
vy          = queorb['vy']
spx         = queorb['spx']
spy         = queorb['spy']
spz         = queorb['spz']
qx          = queorb['qx']
qy          = queorb['qy']
qz          = queorb['qz']
vsys        = queorb['vsys']
Dir         = queorb['Dir']
# spin-orb    = queorb['spin-orb']
PERIC       = queorb['PERIC']
# pos-peri    = queorb['pos-peri']

# Crio diretório para arquivos temporários
dirName = "orbits_temp"
if not os.path.exists(dirName):
    os.mkdir(dirName)
    print("Directory " , dirName ,  " Created ")
else:
    print("Directory " , dirName ,  " already exists")


# Parâmetros para o pot_din_fric
M200  = 57.9E10          # MAIN Galaxy M200 (MakeNewDisk output) [Msun]:
R200  = 200.0            # MAIN Galaxy R200 (MakeNewDisk output) [kpc]:
aHalo = 14.92            # HALO scalelength (aHalo, MakeNewDisk output - RH) [kpc]:
dmfrac= 0.04             # DISK mass fraction (dmfrac):
aDisk = 2.1              # DISK scalelength (aDisk, end of MakeNewDisk output) [kpc]:
Z0    = 0.42             # DISK vertical scalelength (DiskHeight, from MakeNewDisk parameters file) [kpc]:
bmfrac= 0.01             # BULGE mass fraction (bmfrac, from MakeNewDisk parameters file - MD):
aBfrac= 0.019            # BULGE scalelength fraction (aBulge, from MakeNewDisk parameters file - BulgeSize):
bhfrac= 5.8e-5           # BLACK HOLE mass fraction (bhfrac, from MakeNewDisk parameters file - MBH):
# Pot 2 Position (x,y,z) [kpc] --> sai da tabela queorbita pra cada órbita
# Pot 2 Velocity (vx,vy,vz) [km/s]: --> sai da tabela queorbita pra cada órbita
m_p2    = 22.0E10
k_p2    = 5.58          # Pot 2 radial scalelength [kpc]:
tf      = 1000.0         # Final time (tf=0 is pericenter) (Myr)
dt      = 0.1            # Timestep

# Varrendo todas as órbitas BACKWARD
for i in queorb.index:
#for i in [2555]:
    fname = "orb_BACKWARD_%05d" % i
    fname = dirName + "/" + fname
    fileIn = open(fname + ".in", 'w')

    qq = np.array([qx[i], qy[i], qz[i]])      # Vetor de pericentro
    ss = np.array([spx[i], spy[i], spz[i]])   # Vetor de spin
    vecvq = np.cross(ss, qq)*vq[i]            # Veloc. no pericentro: produto vetorial spin X q

    texto = [M200,R200,aHalo,dmfrac,aDisk,Z0,bmfrac,aBfrac,bhfrac,
             str(gx[i]) + " " + str(gy[i]) + " " + str(gz[i]),
             str(VX_ceu[i]) + " " + str(VY_ceu[i]) + " " + str(VZ_ceu[i]),
             m_p2,k_p2,-tf,-dt]
    # Mando pro arquivo de entrada
    for tt in texto:
        fileIn.write(str(tt)+"\n")
    fileIn.close()
    # Rodo pot_hbd_din_fric-2
    fnameIn = fname  + ".in"
    fnameOrbit = fname + ".dat"
    ret=subprocess.call("./pot_hbd_din_fric-2 < " + fnameIn, shell=True, stdout=subprocess.DEVNULL)
    ret=subprocess.call("mv orbits_pot_hbd_din_fric.dat " + fnameOrbit, shell=True, stdout=subprocess.DEVNULL)
    # Criei BACKWARD. Ele vai de time=0 a -1000. Tenho que inverter pra ir de time=-1000 a 0
    orb = pd.read_csv(fnameOrbit,  header=2, sep="\t")  # Leio arquivo da órbita
    fnameFINAL = dirName + "/" + "orb_%05d.dat" % i     # Nome do arquivo de saída
    orb=orb[:0:-1]                                      # Inverto a matriz
    orb.to_csv(fnameFINAL, float_format='%.5f', index=False, sep="\t")



# Varrendo todas as órbitas FORWARD
for i in queorb.index:
#for i in [2555]:
    fname = "orb_FORWARD_%05d" % i
    fname = dirName + "/" + fname
    fileIn = open(fname + ".in", 'w')

    qq = np.array([qx[i], qy[i], qz[i]])      # Vetor de pericentro
    ss = np.array([spx[i], spy[i], spz[i]])   # Vetor de spin
    vecvq = np.cross(ss, qq)*vq[i]            # Veloc. no pericentro: produto vetorial spin X q

    texto = [M200,R200,aHalo,dmfrac,aDisk,Z0,bmfrac,aBfrac,bhfrac,
             str(gx[i]) + " " + str(gy[i]) + " " + str(gz[i]),
             str(VX_ceu[i]) + " " + str(VY_ceu[i]) + " " + str(VZ_ceu[i]),
#             str(qx[i]*q[i]) + " " + str(qy[i]*q[i]) + " " + str(qz[i]*q[i]),
#             str(vecvq[0]) + " " + str(vecvq[1]) + " " + str(vecvq[2]),
             m_p2,k_p2,tf,dt]
    # Mando pro arquivo de entrada
    for tt in texto:
        fileIn.write(str(tt)+"\n")
    fileIn.close()
    # Rodo pot_hbd_din_fric-2
    fnameIn = fname  + ".in"
    fnameOrbit = fname + ".dat"
    ret=subprocess.call("./pot_hbd_din_fric-2 < " + fnameIn, shell=True, stdout=subprocess.DEVNULL)
    ret=subprocess.call("mv orbits_pot_hbd_din_fric.dat " + fnameOrbit, shell=True, stdout=subprocess.DEVNULL)
    # Criei FORWARD. Tenho que meter no arquivo que já foi criado pela rotina do BACKWARD
    orb = pd.read_csv(fnameOrbit,  header=2, sep="\t")  # Leio arquivo da órbita
    fnameFINAL = dirName + "/" + "orb_%05d.dat" % i     # Nome do arquivo de saída
    orb.to_csv(fnameFINAL, float_format='%.5f', index=False, sep="\t", mode='a', header=False)


# Apago os arquivos parciais
ret=subprocess.call("rm -f " + dirName + "/" + "orb_BACKWARD* "  + dirName + "/" + "orb_FORWARD* " + dirName + "/" + "*.in", shell=True, stdout=subprocess.DEVNULL)
