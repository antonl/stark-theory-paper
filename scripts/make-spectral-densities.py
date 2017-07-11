import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from pyqcfp.runqcfp import render_spectral_density

#plt.rcParams.update(
#    {'figure.figsize': (7, 3.67)}
#)

# units in 1/cm
lambda0 = 35
gamma0 = 40
gammaj = 3.

omegah = 1250
gammah = 500
lambdah = 500

omega = np.arange(1e-4, 10000, 0.5)
Ji = 2*lambda0*omega*gamma0/(omega**2 + gamma0**2)

lbd = lambda0

modes = pd.read_csv('../rawdata/oscillator-modes.csv',
                    names=['wj', 'Sj'],
                    delimiter=' ',
                    skiprows=1)
modes = modes.sort_values(by='wj')
print(modes.tail())
Jii = Ji.copy()
Jiii = Ji.copy()
Jnjp = Ji.copy()

omega_g = 105
lambda_g = 170
Jgauss = 2*lambda_g*omega/(np.sqrt(2*np.pi)*omega_g)*np.exp(
     -omega**2/(2*omega_g**2))

omega_low, dw = np.linspace(1e-8, 100, 1000, retstep=True)
lambda_low = 85
gamma_low = 0.33356
Jlowfreq = 2*lambda_low*omega_low*gamma_low/(omega_low**2 + gamma_low**2)

# add a vibrational mode to Ji
V_mode = 150.
V_strength = 0.05
V = Ji.copy()
V_lambda = V_mode*V_strength
ogj = V_mode*gammaj
V += 2*V_lambda*V_mode**2*(ogj/((V_mode**2 - omega**2)**2 + ogj**2))

# superohmic form
So_w = 1.
So_cutoff = 10.
So_s = 3
So_strength = 0.4
Jso = So_strength*(omega/So_w)**(So_s - 1)*omega*np.exp(-omega/So_cutoff)

for idx, (wj, Sj) in modes.iterrows():
    if wj > 1000:
        continue

    lambdaj = Sj*wj
    lbd += lambdaj
    ogj = omega*gammaj
    Jii += 2*lambdaj*wj**2*(ogj/((wj**2 - omega**2)**2 + ogj**2))

lambda0 = 40
Jnjp = 2*lambda0*omega*gamma0/(omega**2 + gamma0**2)
Jnjp += 2*np.sqrt(2)*lambdah*omegah**2*(gammah*omega/((omegah**2 -
                                                       omega**2)**2 +
                                                      2*gammah**2*omega**2))

print(lbd)
print(lambda0 + lambdah)
for idx, (wj, Sj) in modes.iterrows():
    lambdaj = Sj*wj
    lbd += lambdaj
    ogj = omega*gammaj
    Jiii += 2*lambdaj*wj**2*(ogj/((wj**2 - omega**2)**2 + ogj**2))

fig = plt.figure(dpi=150, figsize=(7, 3.67))

ax = fig.add_subplot(111)

ji_line, = ax.plot(omega, Ji, ls='--')
jii_line, = ax.plot(omega, Jii)
#jiii_line, = ax.plot(omega, Jiii, zorder=-1)
#jnjp_line, = ax.plot(omega, Jnjp)
#jlow_line, = ax.plot(omega_low, Jlowfreq)
#jgauss_line, = ax.plot(omega, Jgauss)
jv_line, = ax.plot(omega, V)
jso_line, = ax.plot(omega, Jso)

print(ji_line.get_color())
#print(jii_line.get_color())
#print(jiii_line.get_color())
ax.semilogy()
ax.set_ylim(10e-1, 1e5)
ax.set_xlim(0, 2500)
ax.set_xlabel(r'$\omega$ ($\mathrm{cm}^{-1}$)')
ax.set_ylabel(r"C($\omega$) ($\mathrm{cm}^{-1}$)")
fig.tight_layout()

plt.show()

#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)

'''
from matplotlib.transforms import TransformedBbox, Affine2D
tight_bbox_raw = ax.get_tightbbox(fig.canvas.get_renderer())
tight_bbox_fraw = fig.get_tightbbox(fig.canvas.get_renderer())
tight_bbox = TransformedBbox(tight_bbox_raw, Affine2D().scale(1./fig.dpi)).get_points()
tight_bbox_f = TransformedBbox(tight_bbox_fraw, Affine2D().scale(1./fig.dpi)).get_points()
extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#ax.patch.set_facecolor('lightslategray')
'''
fig.savefig('sims-spectral-densities.png', bbox_inches='tight')

with open('../spectral-densities/Ji-spd.txt', 'w') as f:
    f.write(render_spectral_density(omega, Ji))
with open('../spectral-densities/Jii-spd.txt', 'w') as f:
    f.write(render_spectral_density(omega, Jii))
with open('../spectral-densities/Jnjp-spd.txt', 'w') as f:
    f.write(render_spectral_density(omega, Jnjp))
with open('../spectral-densities/Jlowfreq-spd.txt', 'w') as f:
    f.write(render_spectral_density(omega_low, Jlowfreq))
with open('../spectral-densities/Jgauss-spd.txt', 'w') as f:
    f.write(render_spectral_density(omega, Jgauss))
with open('../spectral-densities/V-spd.txt', 'w') as f:
    f.write(render_spectral_density(omega, V))
with open('../spectral-densities/Jso-spd.txt', 'w') as f:
    f.write(render_spectral_density(omega, Jso))
