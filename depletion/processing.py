###########################################################################################################################################################
#################################################################### Relevant Libraries ###################################################################
###########################################################################################################################################################
import openmc
import openmc.deplete
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colorsf
from scipy import integrate
plt.style.use("/home/baltasar/Documents/Master/plot_parameters.mplstyle")

###########################################################################################################################################################
#################################################################### Calculating Burn-Up ##################################################################
###########################################################################################################################################################
UROX_density                = np.round(0.95 * 10960, 0) # in kg/m3
UROX_volume                 = (np.pi * (0.60/100)**2 * 1.75) * 19 * 144 # in m3
UROX_mass                   = UROX_density * UROX_volume
UROX_molar_mass             = 0.2 * 235.04392 + 0.8 * 238.05078 + 2 * 15.999
uranium_percentage          = (0.2 * 235.04392 + 0.8 * 238.05078)/UROX_molar_mass
UROX_mass                   *= uranium_percentage
###########################################################################################################################################################
print(UROX_density, 'kg/m3')
print(UROX_volume, 'm3')
print(UROX_mass, 'kg')
print(UROX_molar_mass, 'g/mol')
print(uranium_percentage, 'unitless')
print(UROX_mass, 'kg')
###########################################################################################################################################################
power = [5, 10, 15]
energy_extracted = np.array([5, 10, 15])*2000 # in MWd
burn_up = energy_extracted/np.round(UROX_mass,0)
print(np.round(UROX_mass, 0))
print("Power [MWth]".rjust(12), "Burn-up [MWd/kg]".rjust(18))
for n, P in enumerate(power):
    print(f"{P}".rjust(12), f"{burn_up[n]:.2f}".rjust(18))

###########################################################################################################################################################
################################################################## Loading Depletion Files ################################################################
###########################################################################################################################################################
dp_5                        = openmc.deplete.Results("5MW/depletion_results.h5")
dp_10                       = openmc.deplete.Results("10MW/depletion_results.h5")
dp_15                       = openmc.deplete.Results("15MW/depletion_results_old.h5")

###########################################################################################################################################################
################################################################ Plotting k-effective vs time #############################################################
###########################################################################################################################################################
t_5, k_5                    = dp_5.get_keff(time_units = "d")
t_10, k_10                  = dp_10.get_keff(time_units = "d")
t_15, k_15                  = dp_15.get_keff(time_units = "d")
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
# plt.errorbar(t_5, k_5[:,0], k_5[:,1], fmt = ".-", markersize = 2, linewidth = 0.8, color = "blue", label = "5 MWth")
# plt.errorbar(t_10, k_10[:,0], k_10[:,1], fmt = ".-", markersize = 2, linewidth = 0.8, color = "purple", label = "10 MWth")
plt.errorbar(t_15, k_15[:,0], k_15[:,1], fmt = ".-", markersize = 2, linewidth = 0.8, color = "red", label = "15 MWth")
plt.xlabel("Time [days]")
plt.ylabel(r"$k_{eff}$")
# plt.legend()
plt.grid()
# plt.savefig(fname = "15MW_k_vs_time", bbox_inches = "tight")
# plt.show()
###########################################################################################################################################################
for n, t in enumerate(t_15):
    print(f"{n}".ljust(3), f"{t}".ljust(6), f"{k_15[n,0]:.5f}", f"{k_15[n,1]:.5f}", f"{1e5*(1-1/k_15[n,0]):.5f}")

###########################################################################################################################################################
########################################################### Fissile Material Development over Time ########################################################
###########################################################################################################################################################
t_U235, U235            = dp_15.get_atoms("1", nuc = "U235", nuc_units = "atoms", time_units = "d")
U235                    += dp_15.get_atoms("2", nuc = "U235", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_U238, U238            = dp_15.get_atoms("1", nuc = "U238", nuc_units = "atoms", time_units = "d")
U238                    += dp_15.get_atoms("2", nuc = "U238", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_Pu237, Pu237          = dp_15.get_atoms("1", nuc = "Pu237", nuc_units = "atoms", time_units = "d")
Pu237                   += dp_15.get_atoms("2", nuc = "Pu237", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_Pu238, Pu238          = dp_15.get_atoms("1", nuc = "Pu238", nuc_units = "atoms", time_units = "d")
Pu238                   += dp_15.get_atoms("2", nuc = "Pu238", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_Pu239, Pu239          = dp_15.get_atoms("1", nuc = "Pu239", nuc_units = "atoms", time_units = "d")
Pu239                   += dp_15.get_atoms("2", nuc = "Pu239", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_Pu240, Pu240          = dp_15.get_atoms("1", nuc = "Pu240", nuc_units = "atoms", time_units = "d")
Pu240                   += dp_15.get_atoms("2", nuc = "Pu240", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_Pu241, Pu241          = dp_15.get_atoms("1", nuc = "Pu241", nuc_units = "atoms", time_units = "d")
Pu241                   += dp_15.get_atoms("2", nuc = "Pu241", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_Pu242, Pu242          = dp_15.get_atoms("1", nuc = "Pu242", nuc_units = "atoms", time_units = "d")
Pu242                   += dp_15.get_atoms("2", nuc = "Pu242", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_Np237, Np237          = dp_15.get_atoms("1", nuc = "Np237", nuc_units = "atoms", time_units = "d")
Np237                   += dp_15.get_atoms("2", nuc = "Np237", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
t_Pu242, Pu242          = dp_15.get_atoms("1", nuc = "Pu242", nuc_units = "atoms", time_units = "d")
Pu242                   += dp_15.get_atoms("2", nuc = "Pu242", nuc_units = "atoms", time_units = "d")[1]
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
plt.plot(t_U235, U235, "o-", color = "red", markersize = 1.5, linewidth = 1, label = "U235")
plt.xlabel("Time [days]")
plt.ylabel("Number of U235 Atoms")
plt.grid()
# plt.savefig(fname = "U235_vs_time", bbox_inches = "tight")
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
plt.plot(t_Pu239, Pu239, "o-", color = "blue", markersize = 1.5, linewidth = 1, label = "Pu239")
plt.xlabel("Time [days]")
plt.ylabel("Number of Pu239 Atoms")
plt.grid()
# plt.savefig(fname = "Pu239_vs_time", bbox_inches = "tight")
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
plt.plot(t_Pu238, Pu238, "o-", color = "purple", markersize = 1.5, linewidth = 1, label = "Pu238")
plt.plot(t_Pu240, Pu240, "o-", color = "green", markersize = 1.5, linewidth = 1, label = "Pu240")
plt.plot(t_Pu241, Pu241, "o-", color = "magenta", markersize = 1.5, linewidth = 1, label = "Pu241")
plt.xlabel("Time [days]")
plt.ylabel("Number of Atoms")
plt.grid()
plt.legend()
# plt.savefig(fname = "Pu_vs_time", bbox_inches = "tight")
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
plt.plot(t_Np237, Np237, "o-", color = "green", markersize = 1.5, linewidth = 1, label = "Np237")
plt.xlabel("Time [days]")
plt.ylabel("Number of Atoms")
plt.grid()
plt.legend()
# plt.savefig(fname = "MA_vs_time", bbox_inches = "tight")
# plt.show()

###########################################################################################################################################################
############################################################### Plotting neutron flux over time ###########################################################
###########################################################################################################################################################
L                       = 175
t_axial_reflector       = 15
t_radial_reflector      = 20
r_grid, dr              = np.linspace(0, 91 + t_radial_reflector, 100, endpoint = True, retstep = True)
z_grid, dz              = np.linspace(-(L/2 + t_axial_reflector), L/2 + t_axial_reflector, 200, endpoint = True, retstep = True)
energy_range            = np.logspace(-5, np.log10(20e6), num = 251, base = 10.0)
reactor_power           = 15e6
###########################################################################################################################################################
sp_start                = openmc.StatePoint("15MW/openmc_simulation_n0.h5")
sp_end                  = openmc.StatePoint("15MW/openmc_simulation_n24.h5")
###########################################################################################################################################################
colors                  = ["blue", "red"]
labels                  = ["0 days", "2000 days"]
###########################################################################################################################################################
V_reactor = np.pi * (91 + t_radial_reflector)**2 * (175 + t_axial_reflector * 2)
###########################################################################################################################################################
plt.figure(figsize = (5.78, 3.85))
for n, sp in enumerate([sp_start, sp_end]):
    heating             = sp.get_tally(name = "Heating").mean.ravel()
    normalization       = reactor_power/(1.602e-19 * heating * V_reactor)
    flux                = sp.get_tally(name = "Flux").mean.ravel()
    dflux               = sp.get_tally(name = "Flux").std_dev.ravel()
    du                  = np.log(energy_range[1:]/energy_range[:-1])[0]
    integral            = integrate.trapezoid(flux, dx = du)
    # flux                *= normalization
    # dflux               *= normalization
    flux                /= integral
    dflux               /= integral
    plt.errorbar(energy_range[:-1]*1e-6, flux, dflux, fmt = "-", markersize = 1, linewidth = 0.5, color = colors[n], label = labels[n])
plt.xscale("log")
plt.xlim(1e-9, 20)
plt.xlabel("Neutron Energy [MeV]")
# plt.ylabel(r"Neutron Flux [$ns^{-1}cm^{-2}$]")
plt.ylabel("Normalized Neutron Flux")
plt.legend()
plt.grid()
# plt.savefig(fname = "energy_spectra", bbox_inches = "tight")
plt.show()