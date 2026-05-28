###########################################################################################################################################################
#################################################################### Relevant Libraries ###################################################################
###########################################################################################################################################################
import numpy as np
import openmc
import matplotlib.pyplot as plt
from scipy import integrate
from uncertainties import ufloat

###########################################################################################################################################################
#################################################################### Constants/Plotting ###################################################################
###########################################################################################################################################################
from reactor_sim import reflectors, t_axial_inf, t_radial_inf
plt.style.use("/home/baltasar/Documents/Master/plot_parameters.mplstyle")
###########################################################################################################################################################
cm                              = plt.get_cmap("Set1")
colors                          = [cm(1.*n/9) for n in range(3)]
for n, item in enumerate(reflectors):
    if n > 0:
        reflectors[n]            = item.capitalize()

###########################################################################################################################################################
############################################################# Comparing BeO and Graphite Flux #############################################################
###########################################################################################################################################################
r_grid                          = np.linspace(0, 91 + 20, 100, endpoint = True)
r_grid_NR                       = np.linspace(0, 91, 100, endpoint = True)

r_matrix                        = np.zeros((3, len(r_grid)))
r_matrix[:-1]                   = r_grid
r_matrix[-1]                    = r_grid_NR

sp_BeO                          = openmc.StatePoint("flux_files/statepoint_BeO.h5")
sp_graphite                     = openmc.StatePoint("flux_files/statepoint_graphite.h5")
sp_no_ref                       = openmc.StatePoint("flux_files/statepoint_NR.h5")

labels                          = ["BeO", "Graphite", "Unreflected"]

reactor_power                   = 5e6
energy_range                    = np.logspace(-5, np.log10(20e6), num = 1001, base = 10.0)
###########################################################################################################################################################
###########################################################################################################################################################
V_array_rad                     = np.ones((3, len(r_grid)-1))
V_array_rad[:-1]                *= (175) * np.pi * (r_grid[1:]**2 - r_grid[:-1]**2)
V_array_rad[-1]                 *= (175) * np.pi * (r_grid_NR[1:]**2 - r_grid_NR[:-1]**2)
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
for n, sp in enumerate([sp_BeO, sp_graphite, sp_no_ref]):
    heating                     = sp.get_tally(name = "Heating").mean.ravel()
    normalization               = reactor_power/(1.602e-19 * heating * V_array_rad[n])
    radial_flux                 = sp.get_tally(name = "Radial Flux").get_pandas_dataframe()
    flux_mean                   = radial_flux["mean"].to_numpy()
    plt.plot(r_matrix[n,:-1], flux_mean * normalization, 
             markersize = 1, linewidth = 1, color = colors[n], label = labels[n])
plt.xlabel("Radial Position [cm]", size = 14)
plt.ylabel(r"Neutron Flux [$ns^{-1}cm^{-2}$]", size = 14)
plt.legend()
plt.grid()
# plt.savefig(fname = "radial_flux", bbox_inches = "tight")
# plt.show()
###########################################################################################################################################################
###########################################################################################################################################################
V_reactor                       = np.ones(3)
V_reactor[:-1]                  *= np.pi * (91 + 20)**2 * (175)
V_reactor[-1]                   *= np.pi * (91)**2 * (175)
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
for n, sp in enumerate([sp_BeO, sp_graphite, sp_no_ref]):
    heating                     = sp.get_tally(name = "Heating").mean.ravel()
    normalization               = reactor_power/(1.602e-19 * heating * V_reactor[n])
    flux                        = sp.get_tally(name = "Flux").mean.ravel()
    plt.plot(energy_range[:-1]*1e-6, flux * normalization, 
             markersize = 1, linewidth = 1, color = colors[n], label = labels[n])
plt.xscale("log")
plt.xlim(1e-9, 20)
plt.xlabel("Neutron Energy [MeV]", size = 14)
plt.ylabel(r"Neutron Flux [$ns^{-1}cm^{-2}$]", size = 14)
plt.legend()
plt.grid()
# plt.savefig(fname = "energy_spectra", bbox_inches = "tight")
# plt.show()

###########################################################################################################################################################
################################################################# Axial Reflector Analysis ################################################################
###########################################################################################################################################################
k_data                          = np.loadtxt("k_data/axial_k_data")
sigma_k_data                    = np.loadtxt("k_data/axial_dk_data")
k_inf_data                      = np.loadtxt("k_data/inf_axial_k_data")
sigma_k_inf_data                = np.loadtxt("k_data/inf_axial_dk_data")
###########################################################################################################################################################
t_axial_array                   = k_data[1:,0]
rho_data                        = 1 - 1/k_data[1:,1:]
drho_data                       = sigma_k_data[1:,1:]/k_data[1:,1:]**2
rho_data                        *= 10**5
drho_data                       *= 10**5
###########################################################################################################################################################
rho_inf_data                    = 1 - 1/k_inf_data
drho_inf_data                   = sigma_k_inf_data/k_inf_data**2
rho_inf_data                    *= 10**5
drho_inf_data                   *= 10**5
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
for r, reflector in enumerate(reflectors):
    rho_worth                   = np.zeros(len(t_axial_array))
    drho_worth                  = np.zeros(len(t_axial_array))
    for n, rho in enumerate(rho_data[:, r]):
        if n > 0:
            rho_worth[n]        = rho - rho_data[0, r]
            drho_worth[n]       = np.sqrt(drho_data[n, r]**2 + drho_data[0, r]**2)
    plt.errorbar(t_axial_array, rho_worth, drho_worth, 
                 fmt = ".-", linewidth = 1, markersize = 4, color = colors[r], label = str(reflector))
for r, reflector in enumerate(reflectors):
    inf_worth                   = rho_inf_data[r] - rho_data[0, r]
    dinf_worth                  = np.sqrt(drho_inf_data[r]**2 + drho_data[0, r]**2)
    plt.hlines(inf_worth, t_axial_array[0], t_axial_array[-1], 
               linewidth = 1, colors = colors[r], linestyles = '--')
plt.xticks(ticks = np.arange(0, t_axial_array[-1]+5, 5))
plt.ylim(0, 6100)
plt.xlim(-2.5, 57.5)
plt.xlabel("Axial Reflector Thickness [cm]")
plt.ylabel("Reactivity Gain [pcm]")
plt.grid()
plt.legend(title = "Ref. Material")
# plt.savefig(fname = "axial_reflector", bbox_inches = "tight")
# plt.show()

###########################################################################################################################################################
################################################################ Radial Reflector Analysis ################################################################
###########################################################################################################################################################
k_data                          = np.loadtxt("k_data/radial_k_data")
sigma_k_data                    = np.loadtxt("k_data/radial_dk_data")
k_inf_data                      = np.loadtxt("k_data/inf_radial_k_data")
sigma_k_inf_data                = np.loadtxt("k_data/inf_radial_dk_data")
###########################################################################################################################################################
t_radial_array                  = k_data[1:,0]
rho_data                        = 1 - 1/k_data[1:,1:]
drho_data                       = sigma_k_data[1:,1:]/k_data[1:,1:]**2
rho_data                        *= 10**5
drho_data                       *= 10**5
###########################################################################################################################################################
rho_inf_data                    = 1 - 1/k_inf_data
drho_inf_data                   = sigma_k_inf_data/k_inf_data**2
rho_inf_data                    *= 10**5
drho_inf_data                   *= 10**5
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
for r, reflector in enumerate(reflectors):
    rho_worth                   = np.zeros(len(t_radial_array))
    drho_worth                  = np.zeros(len(t_radial_array))
    for n, rho in enumerate(rho_data[:, r]):
        if n > 0:
            rho_worth[n]        = rho - rho_data[0, r]
            drho_worth[n]       = np.sqrt(drho_data[n, r]**2 + drho_data[0, r]**2)
    plt.errorbar(t_radial_array, rho_worth, drho_worth, 
                 fmt = ".-", linewidth = 1, markersize = 4, color = colors[r], label = str(reflector))
for r, reflector in enumerate(reflectors):
    inf_worth                   = rho_inf_data[r] - rho_data[0, r]
    dinf_worth                  = np.sqrt(drho_inf_data[r]**2 + drho_data[0, r]**2)
    plt.hlines(inf_worth, t_radial_array[0], t_radial_array[-1], 
               linewidth = 1, colors = colors[r], linestyles = '--')
plt.xticks(ticks = np.arange(0, t_radial_array[-1]+5, 5))
plt.xlabel("Radial Reflector Thickness [cm]")
plt.ylabel("Reactivity Gain [pcm]")
plt.xlim(-2.5, 57.5)
plt.ylim(bottom = 0, top = 13000)
plt.grid()
plt.legend(title = "Ref. Material")
# plt.savefig(fname = "radial_reflector", bbox_inches = "tight")
plt.show()