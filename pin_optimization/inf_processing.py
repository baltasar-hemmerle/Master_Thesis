#######################################################################################################################################################
#######################################################################################################################################################
################################################### File for Processing Generated Statepoint Files ####################################################
#######################################################################################################################################################
#######################################################################################################################################################

#######################################################################################################################################################
################################################################## Relevant Libraries #################################################################
#######################################################################################################################################################
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("/home/baltasar/Documents/Master/plot_parameters.mplstyle")
import openmc
from scipy import integrate

#######################################################################################################################################################
####################################################################### Plotting ######################################################################
#######################################################################################################################################################
k           = np.loadtxt("k_and_factor_data/k_data")
dk          = np.loadtxt("k_and_factor_data/dk_data")
######################################################################################################################################################
pitch_UROX  = k[1:,0]
radii_UROX  = k[0,1:]
######################################################################################################################################################
N_colors = 8
cm = plt.get_cmap("Set1") # Set1
colors = [cm(1.*n/N_colors) for n in np.arange(0, 10, 1)]
format = ["o-", "v-", "d-", "s-"]*2
markers = ['o', 'v', 'd', 's']*2
######################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
for n, r in enumerate(radii_UROX):
    plt.errorbar(pitch_UROX, k[1:,n+1], dk[1:,n+1], fmt = format[n], linewidth = 1, ms = 2, color = colors[n], label= f'{r:.2f} cm')
plt.ylabel(r"$k_{eff}$")
plt.xlabel("Fuel-HP Pitch [cm]")
plt.xticks(ticks = np.arange(2.4, 5.0, 0.2))
plt.yticks(ticks = np.arange(1.00, 1.30, 0.05))
plt.ylim(1.00, 1.25)
plt.xlim(2.3, 4.8)
plt.grid()
plt.legend(title = "Fuel Radius")
# plt.savefig(fname = "pin_optimization", bbox_inches = "tight")
# plt.show()

#######################################################################################################################################################
#################################################################### Flux Plotting ####################################################################
#######################################################################################################################################################
energy_range    = np.logspace(-5, np.log10(20e6), num = 1001, base = 10.0)
du              = np.log(energy_range[1:]/energy_range[:-1])[0]
sp_45_24        = openmc.StatePoint("flux_data/r45_p24.h5")
sp_45_48        = openmc.StatePoint("flux_data/r45_p48.h5")
sp_60_24        = openmc.StatePoint("flux_data/r60_p24.h5")
sp_60_48        = openmc.StatePoint("flux_data/r60_p48.h5")
sp_80_24        = openmc.StatePoint("flux_data/r80_p24.h5")
sp_80_48        = openmc.StatePoint("flux_data/r80_p48.h5")
colors45        = ['red', 'red']
colors60        = ['blue', 'blue']
colors80        = ['green', 'green']
linestyles      = ["--", "-"]
r_fuel          = [0.45, 0.60, 0.80]
pitch           = [2.4, 4.8]
#######################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
for n, sp in enumerate([sp_45_24, sp_45_48]):
    flux        = sp.get_tally(name = "Flux").mean.ravel()
    integral    = integrate.trapezoid(flux, dx = du)
    flux        /= integral
    plt.plot(energy_range[:-1]*1e-6, flux, linestyles[n], linewidth = 1, color = colors45[n], label = f'r = {r_fuel[0]:.2f}, p = {pitch[n]:.1f}')
for n, sp in enumerate([sp_60_24, sp_60_48]):
    flux        = sp.get_tally(name = "Flux").mean.ravel()
    integral    = integrate.trapezoid(flux, dx = du)
    flux        /= integral
    plt.plot(energy_range[:-1]*1e-6, flux, linestyles[n], linewidth = 1, color = colors60[n], label = f'r = {r_fuel[1]:.2f}, p = {pitch[n]:.1f}')
for n, sp in enumerate([sp_80_24, sp_80_48]):
    flux        = sp.get_tally(name = "Flux").mean.ravel()
    integral    = integrate.trapezoid(flux, dx = du)
    flux        /= integral
    plt.plot(energy_range[:-1]*1e-6, flux, linestyles[n], linewidth = 1, color = colors80[n], label = f'r = {r_fuel[2]:.2f}, p = {pitch[n]:.1f}')
plt.xscale("log")
plt.xlim(1e-8, 20)
plt.xlabel("Neutron Energy [MeV]")
plt.ylabel("Normalized Neutron Flux") # per unit lethargy
plt.legend(ncols = 1, alignment = 'right', fontsize = 'small')
plt.grid()
# plt.savefig(fname = "pin_spectra", bbox_inches = "tight")
plt.show()