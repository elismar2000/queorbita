from astropy import table
import matplotlib.pyplot as plt

#A unica coisa que precisa mudar possivelmente é o caminho/nome do arquivo com a saída do queorbita na linha abaixo
t = table.Table.read('../orbits_3rd_attempt/queorbita_new.out', data_start=50, header_start=49, format='ascii')

pos = t['PERIC'] == 'POS'
# mask = (t['PERIC'] == 'POS') & (t['spin-orb'] < 50) & (t['e'] == 0.9) & (t['VZ_ceu'] > -119 - 2) & (t['VZ_ceu'] < -119 + 2)

t = t[pos]

plt.style.use('Solarize_Light2')

fig, axs = plt.subplots(2, 5, figsize=(12, 7))

axs[0, 0].hist(t['e'], bins=50, color='gold')
axs[0, 0].set_title('excentricity')
axs[0, 0].set_ylabel('number of orbits')
axs[0, 0].set_xlabel('excentricity')
axs[0, 0].tick_params(length=2, tickdir='in', labelsize=8)


axs[0, 1].hist(t['gz'], bins=50, color='gold')
axs[0, 1].set_title('gz')
axs[0, 1].set_xlabel('kpc')
axs[0, 1].tick_params(length=2, tickdir='in', labelsize=8)
axs[0, 1].set_ylim((0, 1100))


axs[0, 2].hist(t['q'], bins=50, color='gold')
axs[0, 2].set_title('distance of pericenter')
axs[0, 2].set_xlabel('kpc')
axs[0, 2].tick_params(length=2, tickdir='in', labelsize=8)


axs[0, 3].hist(t['rnow'], bins=50, color='gold')
axs[0, 3].set_title('r_now')
axs[0, 3].set_xlabel('kpc')
axs[0, 3].tick_params(length=2, tickdir='in', labelsize=8)


axs[0, 4].hist(t['vq'], bins=50, color='gold')
axs[0, 4].set_title('velocity of pericenter')
axs[0, 4].set_xlabel('km/s')
axs[0, 4].tick_params(length=2, tickdir='in', labelsize=8)


axs[1, 0].hist(t['vx'], bins=50, color='gold')
axs[1, 0].set_title('vx')
axs[1, 0].set_ylabel('number of orbits')
axs[1, 0].set_xlabel('km/s')
axs[1, 0].tick_params(length=2, tickdir='in', labelsize=8)


axs[1, 1].hist(t['vy'], bins=50, color='gold')
axs[1, 1].set_title('vy')
axs[1, 1].set_xlabel('km/s')
axs[1, 1].tick_params(length=2, tickdir='in', labelsize=8)


axs[1, 2].hist(t['VX_ceu'], bins=50, color='gold')
axs[1, 2].set_title('VX_ceu')
axs[1, 2].set_xlabel('km/s')
axs[1, 2].tick_params(length=2, tickdir='in', labelsize=8)


axs[1, 3].hist(t['VY_ceu'], bins=50, color='gold')
axs[1, 3].set_title('VY_ceu')
axs[1, 3].set_xlabel('km/s')
axs[1, 3].tick_params(length=2, tickdir='in', labelsize=8)


axs[1, 4].hist(t['VZ_ceu'], bins=50, color='gold')
axs[1, 4].set_title('VZ_ceu')
axs[1, 4].set_xlabel('km/s')
axs[1, 4].tick_params(length=2, tickdir='in', labelsize=8)

plt.subplots_adjust(wspace=0.3, hspace=0.3)
# plt.savefig('histograms.png', dpi=150)
plt.show()
