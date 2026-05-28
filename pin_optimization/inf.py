###########################################################################################################################
#################################################### Relevant Libraries ###################################################
###########################################################################################################################
import numpy as np
import openmc
import matplotlib.pyplot as plt
import os
plt.style.use("/home/baltasar/Documents/Master/plot_parameters.mplstyle")

###########################################################################################################################
###################################################### Model Function #####################################################
###########################################################################################################################
def make_model(r_fuel : float, t_fuel_gap : float, t_fuel_clad : float, 
               r_HP : float, t_HP_clad : float, L_HP : float, 
               pitch : float, 
               T_array : np.ndarray[any, float], enrichment : float, shannon_entropy : bool = False):
    '''
    Function for creating xml input files for OpenMC k-eigenvalue calculation. \n 
    Model created by this function corresponds to 2 uranium-dioxide fuel pins in Zircaloy-4 cladding and 1 sodium heat pipe with SS316 cladding 
    arranged in a hexagonal unit cell with reflective boundary conditions. \n 
    Can be used to investigate the relationship between fuel pin paramters, pitch, and the neutron multiplication factor

    Parameters
    ----------
    r_fuel : radius of the UROX fuel pin [cm]
    t_fuel_gap : thickness of the fuel-cladding helium gap [cm]
    t_fuel_clad : thickness of the zircaloy fuel cladding [cm]
    r_HP : radius of the heat pipe vapor core [cm]
    t_HP_clad : thickness of the heat pipe wick/cladding structure [cm]
    L_HP : length of the heat pipe evaporator section/effective core height [cm]
    pitch : fuel/heat pipe pitch [cm]
    T_array : list/array containing the material temperatures (in Kelvin) [fuel, fuel-gap, fuel-cladding, block, HP]
    enrichment : float between 0 and 1 specifying the U-235 enrichment (atomic percent)
    shannon_entropy : boolean indiciating whether to calculate the shannon entropy or not. Default is false

    Returns
    -------
    Exports the materials.xml, geometry.xml, settings.xml, and tallies.xml files representing the model to the directory in which the script is run
    '''
    #######################################################################################################################
    ###################################################### Materials ######################################################
    #######################################################################################################################
    materials_list                  = []
    #######################################################################################################################
    UROX                            = openmc.Material(name = "UROX")
    UROX.add_element("O", 2)
    UROX.add_nuclide("U238", 1 - enrichment)
    UROX.add_nuclide("U235", enrichment)
    UROX.add_s_alpha_beta("c_U_in_UO2")
    UROX.set_density("kg/m3", np.round(0.95 * 10960, 0))
    UROX.temperature                = T_array[0]
    materials_list.append(UROX)
    #######################################################################################################################
    helium                          = openmc.Material(name = "Gap Helium")
    helium.add_element("He", 1)
    helium.temperature              = T_array[1]
    helium.set_density("kg/m3", 0.973)
    materials_list.append(helium)
    #######################################################################################################################
    zircaloy                        = openmc.Material(name = "Zircaloy-4")
    zircaloy.add_element("Sn", 1.45/100, "wo")
    zircaloy.add_element("Fe", 0.21/100, "wo")
    zircaloy.add_element("Cr", 0.10/100, "wo")
    zircaloy.add_element("O", 0.125/100, "wo")
    zircaloy.add_element("Zr", (100-1.45-0.21-0.10-0.125)/100, "wo")
    zircaloy.set_density("kg/m3", 6501)
    zircaloy.temperature            = T_array[2]
    materials_list.append(zircaloy)
    #######################################################################################################################
    graphite                        = openmc.Material(name = "G-348 Graphite")
    graphite.add_element("C", 1)
    graphite.set_density("kg/m3", 1800)
    graphite.add_s_alpha_beta("c_Graphite")
    graphite.temperature            = T_array[3]
    materials_list.append(graphite)
    #######################################################################################################################
    SS316                           = openmc.Material(name = "Stainless Steel 316")
    SS316.add_element('Fe', 0.61575)
    SS316.add_element('C', 0.0008)
    SS316.add_element('Cr', 0.18)
    SS316.add_element('Ni', 0.14)
    SS316.add_element('Mn', 0.02)
    SS316.add_element('Mo', 0.03)
    SS316.add_element('Si', 0.01)
    SS316.add_element('P', 0.00045)
    SS316.add_element('S', 0.0003)
    SS316.set_density('kg/m3', np.round(0.80*7954,0))
    SS316.temperature               = T_array[4]
    materials_list.append(SS316)
    #######################################################################################################################
    sodium                          = openmc.Material(name = "Sodium")
    sodium.add_element("Na", 1)
    sodium.set_density("kg/m3", np.round(0.80 * 1.70 * 10**-2 + 0.20*805, 2))
    sodium.temperature              = T_array[4]
    materials_list.append(sodium)
    #######################################################################################################################
    materials                       = openmc.Materials(materials_list)
    #######################################################################################################################
    ###################################################### Geometry #######################################################
    #######################################################################################################################
    z_bottom                        = openmc.ZPlane(z0 = -L_HP/2,
                                                    boundary_type = "vacuum")
    z_top                           = openmc.ZPlane(z0 = L_HP/2,
                                                    boundary_type = "vacuum")
    #######################################################################################################################
    fuel_cylinder                   = openmc.ZCylinder(r = r_fuel,
                                                       boundary_type = "transmission")
    gap_cylinder                    = openmc.ZCylinder(r = r_fuel + t_fuel_gap,
                                                       boundary_type = "transmission")
    clad_cylinder                   = openmc.ZCylinder(r = r_fuel + t_fuel_gap + t_fuel_clad, 
                                                       boundary_type = "transmission")
    #######################################################################################################################
    HP_cylinder                     = openmc.ZCylinder(r = r_HP,
                                                       boundary_type = "transmission")
    HP_clad_cylinder                = openmc.ZCylinder(r = r_HP + t_HP_clad,
                                                       boundary_type = "transmission")
    #######################################################################################################################
    fuel_cell                       = openmc.Cell(fill = UROX,
                                                  region = - fuel_cylinder &+ z_bottom &- z_top)
    gap_cell                        = openmc.Cell(fill = helium,
                                                  region = - gap_cylinder &+ fuel_cylinder &+ z_bottom &- z_top)
    clad_cell                       = openmc.Cell(fill = zircaloy,
                                                  region = - clad_cylinder &+ gap_cylinder &+ z_bottom &- z_top)
    fuel_outer_cell                 = openmc.Cell(fill = graphite,
                                                  region = + clad_cylinder &+ z_bottom &- z_top)
    HP_cell                         = openmc.Cell(fill = sodium,
                                                  region = - HP_cylinder &+ z_bottom &- z_top)
    HP_clad_cell                    = openmc.Cell(fill = SS316,
                                                  region = - HP_clad_cylinder &+ HP_cylinder &+ z_bottom &- z_top)
    HP_outer_cell                   = openmc.Cell(fill = graphite,
                                                  region = + HP_clad_cylinder &+ z_bottom &- z_top)
    #######################################################################################################################
    fuel_pin                        = openmc.Universe(name = "Fuel Rod", 
                                                      cells = [fuel_cell, gap_cell, clad_cell, fuel_outer_cell])
    HP                              = openmc.Universe(name = "Heat Pipe",  
                                                      cells = [HP_cell, HP_clad_cell, HP_outer_cell])
    #######################################################################################################################
    outer_cell                      = openmc.Cell(fill = graphite)
    outer                           = openmc.Universe(cells = [outer_cell])
    #######################################################################################################################
    lat                             = openmc.HexLattice()
    lat.center                      = (0, 0)
    lat.pitch                       = (pitch, )
    lat.orientation                 = 'y'
    lat.outer                       = outer
    lat.universes                   = [[fuel_pin, HP]*3, [fuel_pin]]
    #######################################################################################################################
    hex                             = openmc.model.HexagonalPrism(pitch, 
                                                                  orientation = 'y', 
                                                                  boundary_type = "reflective")
    lat_region                      = - hex &+ z_bottom &- z_top
    #######################################################################################################################
    reactor_cell                    = openmc.Cell(fill = lat, 
                                                  region = lat_region)
    reactor_universe                = openmc.Universe(cells = [reactor_cell, outer_cell])
    #######################################################################################################################
    geometry                        = openmc.Geometry()
    geometry.root_universe          = reactor_universe

    #######################################################################################################################
    ###################################################### Settings #######################################################
    #######################################################################################################################
    settings                        = openmc.Settings()
    settings.photon_transport       = True
    settings.batches                = 200
    settings.inactive               = 50
    settings.particles              = 10000
    settings.generations_per_batch  = 3
    settings.output                 = {'tallies': False}
    settings.temperature            = {'method': 'interpolation', 'range': (250, 1200)}
    bounds                          = [-pitch/2, -pitch/2, -L_HP/2, 
                                       pitch/2, pitch/2, L_HP/2] 
    uniform_dist                    = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source                 = openmc.IndependentSource(space = uniform_dist, 
                                                               constraints = {'domains': [reactor_cell]})
    if shannon_entropy == False:
        pass
    else:
        entropy_mesh                = openmc.RegularMesh()
        entropy_mesh.lower_left     = [-pitch/2, -pitch/2, -L_HP/2]    
        entropy_mesh.upper_right    = [pitch/2, pitch/2, L_HP/2]
        entropy_mesh.dimension      = (25, 25, 25)
        settings.entropy_mesh       = entropy_mesh

    #######################################################################################################################
    ####################################################### Tallies #######################################################
    #######################################################################################################################
    tallies                         = openmc.Tallies() 
    #######################################################################################################################
    particle_filter                 = openmc.ParticleFilter('neutron')
    energy_range                    = np.logspace(-5, np.log10(20e6), num = 1001, base = 10.0)
    energy_filter                   = openmc.EnergyFilter(energy_range)
    ########################################################################################################################
    energy_spectra                  = openmc.Tally(name = "Flux")
    energy_spectra.filters          = [energy_filter, particle_filter]
    energy_spectra.scores           = ["flux"]
    tallies.append(energy_spectra)

    #######################################################################################################################
    #################################################### Making Model #####################################################
    #######################################################################################################################
    model                           = openmc.model.Model(geometry = geometry, materials = materials, 
                                                         settings = settings, tallies = tallies)
    model.export_to_xml()

    #######################################################################################################################
    #################################################### Making Plots #####################################################
    #######################################################################################################################
    # plt.figure(figsize = (5.78, 5.78))
    # geometry.plot(basis = "xy", pixels = [1200, 1200], width = (5.2, 5.2), 
    #                 origin = (0, 0, 0), 
    #                 color_by = "material", 
    #                 colors = {UROX: "saddlebrown", helium: "white", zircaloy: "silver", 
    #                           graphite: "black", SS316: "darkorange", sodium: "mediumblue" })
    # plt.xlabel('x [cm]')
    # plt.ylabel('y [cm]')
    # plt.savefig(fname = "xy_pin_optimization", bbox_inches = "tight")

###########################################################################################################################
################################################# General Model Parameters ################################################
###########################################################################################################################
enrichment                          = 0.20
###########################################################################################################################
t_fuel_gap                          = 0.01
t_fuel_clad                         = 0.09
###########################################################################################################################
r_HP                                = 0.85
t_HP_clad                           = 0.35
L_HP                                = 175
###########################################################################################################################
T_fuel                              = 1000
T_gap                               = 1000
T_cladding                          = 1000
T_graphite                          = 1000
T_HP                                = 900
T_array                             = [T_fuel, T_gap, T_cladding, T_graphite, T_HP]
###########################################################################################################################
OMC_threads                         = 14 # number of OpenMP threads to use

###########################################################################################################################
################################################## k-effective Calculation ################################################
###########################################################################################################################
fuel_radii                          = np.arange(0.45, 0.85, 0.05)
pitch_array                         = np.arange(2*(r_HP+t_HP_clad), 4*(r_HP+t_HP_clad), 0.1)
###########################################################################################################################
k_data                              = np.zeros((len(pitch_array)+1, len(fuel_radii)+1))
dk_data                             = np.zeros((len(pitch_array)+1, len(fuel_radii)+1))
###########################################################################################################################
k_data[0,1:]                        = fuel_radii
dk_data[0,1:]                       = fuel_radii
k_data[1:, 0]                       = pitch_array
dk_data[1:, 0]                      = pitch_array
###########################################################################################################################
if __name__ == "__main__":
    for j, r in enumerate(fuel_radii):
        for i, pitch in enumerate(pitch_array):
            print(f'Radius = {r}, Pitch = {pitch}')
            make_model(r, t_fuel_gap, t_fuel_clad, r_HP, t_HP_clad, L_HP, pitch, T_array, enrichment)
            openmc.run(threads = OMC_threads)
            with openmc.StatePoint('statepoint.200.h5') as sp:
                k                   = sp.keff
                k_data[i+1, j+1]    = k.nominal_value
                dk_data[i+1, j+1]   = k.std_dev
    #######################################################################################################################
    np.savetxt('k_data', k_data)
    np.savetxt('dk_data', dk_data)

###########################################################################################################################
#################################################### Spectra Calculation ##################################################
###########################################################################################################################
# if __name__ == "__main__":
#     for r in [0.45, 0.60, 0.80]:
#         for p in [2.4, 3.6, 4.8]:
#                 print(f"r = {r}", f"p = {p}")
#                 make_model(r, t_fuel_gap, t_fuel_clad, r_HP, t_HP_clad, L_HP, p, T_array, enrichment)
#                 openmc.run(threads = 16)
#                 fname = 'r'+str(int(r*100))+'_p'+str(int(p*10))+".h5"
#                 os.rename('statepoint.200.h5', fname)