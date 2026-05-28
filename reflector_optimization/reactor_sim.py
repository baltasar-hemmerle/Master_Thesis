###########################################################################################################################################################
#################################################################### Relevant Libraries ###################################################################
###########################################################################################################################################################
import numpy as np
import openmc
import matplotlib.pyplot as plt
from typing import Literal
plt.style.use("/home/baltasar/Documents/Master/plot_parameters.mplstyle")

###########################################################################################################################################################
###################################################################### Model Function #####################################################################
###########################################################################################################################################################
def make_model(r_fuel: float, t_fuel_gap : float, t_fuel_clad : float, 
               r_HP : float, t_HP_clad : float, 
               L_HP : float, pitch : float,
               T_array : np.ndarray[any, float],
               reflector : Literal["BeO", "graphite"],
               enrichment : float,
               t_axial_reflector : float, t_radial_reflector : float,
               shannon_entropy : bool = False):
    '''
    Function for creating xml input files for OpenMC k-eigenvalue calculation. \n 
    Model created by this function consists of UROX fuel and sodium heat pipes arranged in hexagonal graphite blocks to form assemblies. 19 of these assemblies
    are arranged to form a core. This core can be reflected or unreflected depending on the paramters provided to the function.

    Parameters
    ----------
    r_fuel : radius of the UROX fuel pin [cm]
    t_fuel_gap : thickness of the fuel-cladding helium gap [cm]
    t_fuel_clad : thickness of the zircaloy fuel cladding [cm]
    r_HP : radius of the heat pipe vapor core [cm]
    t_HP_clad : thickness of the heat pipe wick/cladding structure [cm]
    L_HP : length of the heat pipe evaporator section/effective core height [cm]
    pitch : fuel/heat pipe pitch [cm]
    T_array : list/array containing the material temperatures (in Kelvin) [fuel, fuel-gap, fuel-cladding, block, HP, reflector]
    reflector : string specifying reflector material. Must be "BeO" for beryllium-oxide or "graphite" for graphite reflector
    enrichment : float between 0 and 1 specifying the U-235 enrichment (atomic percent)
    t_axial_reflector : axial reflector thickness [cm]
    t_radial_reflector : radial reflector thickness [cm]
    shannon_entropy : boolean indiciating whether to calculate the shannon entropy or not. Default is false

    Returns
    -------
    Exports the materials.xml, geometry.xml, settings.xml, and tallies.xml files representing the model to the directory in which the script is run
    '''
    #####################################################################################################################################################
    ##################################################################### Materials #####################################################################
    #####################################################################################################################################################
    materials_list = []
    #####################################################################################################################################################
    UROX = openmc.Material(name = "UROX")
    UROX.add_element("O", 2)
    UROX.add_nuclide("U238", 1 - enrichment)
    UROX.add_nuclide("U235", enrichment)
    UROX.add_s_alpha_beta("c_U_in_UO2")
    UROX.set_density("kg/m3", np.round(0.95 * 10960, 0))
    UROX.temperature = T_array[0]
    materials_list.append(UROX)
    #####################################################################################################################################################
    helium = openmc.Material(name = "Gap Helium")
    helium.add_element("He", 1)
    helium.set_density("kg/m3", 0.973)
    helium.temperature = T_array[1]
    materials_list.append(helium)
    #####################################################################################################################################################
    zircaloy = openmc.Material(name = "Zircaloy-4")
    zircaloy.add_element("Sn", 1.45/100, "wo")
    zircaloy.add_element("Fe", 0.21/100, "wo")
    zircaloy.add_element("Cr", 0.10/100, "wo")
    zircaloy.add_element("O", 0.125/100, "wo")
    zircaloy.add_element("Zr", (100-1.45-0.21-0.10-0.125)/100, "wo")
    zircaloy.set_density("kg/m3", 6501)
    zircaloy.temperature = T_array[2]
    materials_list.append(zircaloy)
    #####################################################################################################################################################
    graphite = openmc.Material(name = "G-348 Graphite")
    graphite.add_element("C", 1)
    graphite.set_density("kg/m3", 1800)
    graphite.add_s_alpha_beta("c_Graphite")
    graphite.temperature = T_array[3]
    materials_list.append(graphite)
    #####################################################################################################################################################
    SS316 = openmc.Material(name = "Stainless Steel 316")
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
    SS316.temperature = T_array[4]
    materials_list.append(SS316)
    #####################################################################################################################################################
    sodium = openmc.Material(name = "Sodium")
    sodium.add_element("Na", 1)
    sodium.set_density("kg/m3", np.round(0.80 * 1.70 * 10**-2 + 0.20*805, 2))
    sodium.temperature = T_array[4]
    materials_list.append(sodium)
    #####################################################################################################################################################
    BeO = openmc.Material(name = "BeO")
    BeO.add_element("O", 1)
    BeO.add_element("Be", 1)
    BeO.set_density("kg/m3", 3005)
    BeO.add_s_alpha_beta("c_Be_in_BeO")
    materials_list.append(BeO)
    #####################################################################################################################################################
    materials = openmc.Materials(materials_list)
    #####################################################################################################################################################
    if reflector == "graphite":
        ref_mat = graphite
    elif reflector == "BeO":
        ref_mat = BeO

    #####################################################################################################################################################
    ################################################################### Pin Geometry ####################################################################
    #####################################################################################################################################################
    if t_axial_reflector == 0:
        pin_boundary = "vacuum"
    else:
        pin_boundary = "transmission"
    #####################################################################################################################################################
    z_core_bottom                       = openmc.ZPlane(z0 = -L_HP/2 - t_axial_reflector,   boundary_type = "vacuum")
    z_core_top                          = openmc.ZPlane(z0 = L_HP/2 + t_axial_reflector,    boundary_type = "vacuum")
    z_fuel_bottom                       = openmc.ZPlane(z0 = -L_HP/2,                       boundary_type = pin_boundary)
    z_fuel_top                          = openmc.ZPlane(z0 = L_HP/2,                        boundary_type = pin_boundary)
    #####################################################################################################################################################
    fuel_cylinder                       = openmc.ZCylinder(r = r_fuel)
    gap_cylinder                        = openmc.ZCylinder(r = r_fuel + t_fuel_gap)
    clad_cylinder                       = openmc.ZCylinder(r = r_fuel + t_fuel_gap + t_fuel_clad)
    #####################################################################################################################################################
    fuel_cell                           = openmc.Cell(fill = UROX,          region = - fuel_cylinder &+ z_fuel_bottom &- z_fuel_top)
    gap_cell                            = openmc.Cell(fill = helium,        region = + fuel_cylinder &- gap_cylinder &+ z_fuel_bottom &- z_fuel_top)
    clad_cell                           = openmc.Cell(fill = zircaloy,      region = + gap_cylinder &- clad_cylinder &+ z_fuel_bottom &- z_fuel_top)
    outer_cell                          = openmc.Cell(fill = graphite,      region = + clad_cylinder &+ z_fuel_bottom &- z_fuel_top)
    #####################################################################################################################################################
    HP_cylinder                         = openmc.ZCylinder(r = r_HP)
    HP_clad_cylinder                    = openmc.ZCylinder(r = r_HP + t_HP_clad)
    #####################################################################################################################################################
    HP_cell                             = openmc.Cell(fill = sodium,    region = - HP_cylinder &+ z_fuel_bottom &- z_fuel_top)
    HP_clad_cell                        = openmc.Cell(fill = SS316,     region = + HP_cylinder &- HP_clad_cylinder &+ z_fuel_bottom &- z_fuel_top)
    HP_outer_cell                       = openmc.Cell(fill = graphite,  region = + HP_clad_cylinder &+ z_fuel_bottom &- z_fuel_top)
    #####################################################################################################################################################
    HP_cell_t                           = openmc.Cell(fill = sodium,    region = - HP_cylinder &+ z_fuel_top &- z_core_top)
    HP_clad_cell_t                      = openmc.Cell(fill = SS316,     region = + HP_cylinder &- HP_clad_cylinder &+ z_fuel_top &- z_core_top)
    HP_top_reflector                    = openmc.Cell(fill = ref_mat,   region = + HP_clad_cylinder &+ z_fuel_top &- z_core_top)
    HP_top_reflector.temperature        = T_array[5]
    #####################################################################################################################################################
    top_reflector_cell                  = openmc.Cell(fill = ref_mat,   region = + z_fuel_top &- z_core_top)
    top_reflector_cell.temperature      = T_array[5]
    #####################################################################################################################################################
    fuel_pin                            = openmc.Universe(name = "Fuel Rod",           cells = [fuel_cell, gap_cell, clad_cell, outer_cell])
    HP                                  = openmc.Universe(name = "Heat Pipe",          cells = [HP_cell, HP_clad_cell, HP_outer_cell])
    HP_t                                = openmc.Universe(name = "Upper Heat Pipe",    cells = [HP_cell_t, HP_clad_cell_t, HP_top_reflector])
    ref_pin                             = openmc.Universe(name = "Top Reflector Pin",  cells = [top_reflector_cell])

    #####################################################################################################################################################
    ################################################################# Assembly Geometry #################################################################
    #####################################################################################################################################################
    main_outer_lattice_cell             = openmc.Cell(fill = graphite)
    top_outer_lattice_cell              = openmc.Cell(fill = ref_mat)
    top_outer_lattice_cell.temperature  = T_array[5]
    main_outer_lattice                  = openmc.Universe(cells = [main_outer_lattice_cell])
    top_outer_lattice                   = openmc.Universe(cells = [top_outer_lattice_cell])
    ################################################################# Hexagonal Lattice #################################################################
    main_lat                            = openmc.HexLattice()
    main_lat.center                     = (0, 0)
    main_lat.pitch                      = (pitch, )
    main_lat.orientation                = 'x'
    main_lat.outer                      = main_outer_lattice
    #####################################################################################################################################################
    top_lat                             = openmc.HexLattice()
    top_lat.center                      = (0, 0)
    top_lat.pitch                       = (pitch, )
    top_lat.orientation                 = 'x'
    top_lat.outer                       = top_outer_lattice
    #####################################################################################################################################################
    main_ring_list                      = [[HP], [fuel_pin]*6]
    top_ring_list                       = [[HP_t], [ref_pin]*6]
    #####################################################################################################################################################
    ring3_main                          = [HP if n in range(2, 14, 2) else fuel_pin for n in range(1,13)]
    ring3_top                           = [HP_t if n in range(2, 14, 2) else ref_pin for n in range(1,13)]
    main_ring_list.append(ring3_main)
    top_ring_list.append(ring3_top)
    #####################################################################################################################################################
    ring4_main                          = [HP if n in range(1, 19, 3) else fuel_pin for n in range(1,19)]
    ring4_top                           = [HP_t if n in range(1, 19, 3) else ref_pin for n in range(1,19)]
    main_ring_list.append(ring4_main)
    top_ring_list.append(ring4_top)
    #####################################################################################################################################################
    ring5_main                          = [HP if n in range(3, 27, 4) else fuel_pin for n in range(1,25)]
    ring5_top                           = [HP_t if n in range(3, 27, 4) else ref_pin for n in range(1,25)]
    main_ring_list.append(ring5_main)
    top_ring_list.append(ring5_top)
    #####################################################################################################################################################
    ring6_main                          = [HP if n in [2, 5, 7, 10, 12, 15, 17, 20, 22, 25, 27, 30] else fuel_pin for n in range(1,31)]
    ring6_top                           = [HP_t if n in [2, 5, 7, 10, 12, 15, 17, 20, 22, 25, 27, 30] else ref_pin for n in range(1,31)]
    main_ring_list.append(ring6_main)
    top_ring_list.append(ring6_top)
    #####################################################################################################################################################
    ring7_main                          = [HP if n in range(1, 37, 3) else fuel_pin for n in range(1,37)]
    ring7_top                           = [HP_t if n in range(1, 37, 3) else ref_pin for n in range(1,37)]
    main_ring_list.append(ring7_main)
    top_ring_list.append(ring7_top)
    #####################################################################################################################################################
    ring8_main                          = [HP if n in [3, 6, 10, 13, 17, 20, 24, 27, 31, 34, 38, 41] else fuel_pin for n in range(1,43)]
    ring8_top                           = [HP_t if n in [3, 6, 10, 13, 17, 20, 24, 27, 31, 34, 38, 41] else ref_pin for n in range(1,43)]
    main_ring_list.append(ring8_main)
    top_ring_list.append(ring8_top)
    #####################################################################################################################################################
    ring9_main                          = [HP if n in [2, 5, 8, 10, 13, 16, 18, 21, 24, 26, 29, 32, 34, 37, 40, 42, 45, 48] else fuel_pin for n in range(1,49)]
    ring9_top                           = [HP_t if n in [2, 5, 8, 10, 13, 16, 18, 21, 24, 26, 29, 32, 34, 37, 40, 42, 45, 48] else ref_pin for n in range(1,49)]
    main_ring_list.append(ring9_main)
    top_ring_list.append(ring9_top)
    #####################################################################################################################################################
    main_ring_list.reverse()
    main_lat.universes                  = main_ring_list
    top_ring_list.reverse()
    top_lat.universes                   = top_ring_list
    ############################################################### Surfaces for Assembly ###############################################################
    hex                                 = openmc.model.HexagonalPrism(pitch*8.7, orientation = 'x', boundary_type = "transmission")
    ############################################################### Regions for Assembly ################################################################
    main_core_region                    = - hex &+ z_fuel_bottom &- z_fuel_top
    top_core_region                     = - hex &+ z_fuel_top &- z_core_top
    main_outer_region                   = + hex &+ z_fuel_bottom &- z_fuel_top
    top_outer_region                    = + hex &+ z_fuel_top &- z_core_top
    lower_ref_region                    = + z_core_bottom &- z_fuel_bottom
    ########################################################## Making cells and Returning Univere ########################################################
    main_assembly_cell                  = openmc.Cell(fill = main_lat, region = main_core_region)
    main_outer_cell                     = openmc.Cell(fill = graphite, region = main_outer_region)
    top_assembly_cell                   = openmc.Cell(fill = top_lat,  region = top_core_region)
    top_outer_cell                      = openmc.Cell(fill = ref_mat, region = top_outer_region)
    top_outer_cell.temperature          = T_array[5]
    lower_ref_cell                      = openmc.Cell(fill = ref_mat, region = lower_ref_region)
    lower_ref_cell.temperature          = T_array[5]
    assembly                            = openmc.Universe(cells = [main_assembly_cell, top_assembly_cell, main_outer_cell, top_outer_cell, lower_ref_cell])

    #####################################################################################################################################################
    ################################################################## Reactor Geometry #################################################################
    #####################################################################################################################################################
    lr_outer_cell                       = openmc.Cell(fill = ref_mat, region = - z_fuel_bottom &+ z_core_bottom)
    lr_outer_cell.temperature           = T_array[5]
    reactor_outer_cell                  = openmc.Cell(fill = graphite, region = + z_fuel_bottom &- z_fuel_top)
    tr_outer_cell                       = openmc.Cell(fill = ref_mat, region = + z_fuel_top &- z_core_top)
    tr_outer_cell.temperature           = T_array[5]
    outer_reactor_lattice               = openmc.Universe(cells = [lr_outer_cell, reactor_outer_cell, tr_outer_cell])
    ################################################################## Hexagonal Lattice ##################################################################
    lat                                 = openmc.HexLattice()
    lat.center                          = (0, 0)
    lat.pitch                           = (pitch*8.75*np.sqrt(3), )
    lat.orientation                     = 'y'
    lat.outer                           = outer_reactor_lattice
    #######################################################################################################################################################
    ring_list                           = [[assembly]*12, [assembly]*6, [assembly]]
    lat.universes                       = ring_list
    ################################################################## Surfaces for Core ##################################################################
    if t_radial_reflector == 0:
        core_boundary = "vacuum"
    else:
        core_boundary = "transmission"
    #######################################################################################################################################################
    core_cylinder                       = openmc.ZCylinder(r = 91,                         boundary_type = core_boundary)
    reflector_cylinder                  = openmc.ZCylinder(r = 91 + t_radial_reflector,    boundary_type = "vacuum") 
    ################################################################## Regions for Core ###################################################################
    lat_region                          = - core_cylinder &+ z_core_bottom &- z_core_top
    reflector_region                    = + core_cylinder &- reflector_cylinder &+ z_core_bottom &- z_core_top
    ########################################################## Making cells and Returning Univere #########################################################
    core_cell                           = openmc.Cell(fill = lat, region = lat_region)
    reflector_cell                      = openmc.Cell(fill = ref_mat, region = reflector_region)
    reflector_cell.temperature          = T_array[5]
    reactor_cells                       = [core_cell, reflector_cell]
    reactor_universe                    = openmc.Universe(name = "Reactor Core", cells = reactor_cells)

    #######################################################################################################################################################
    ################################################################## Geomtry XML File ###################################################################
    #######################################################################################################################################################
    geometry                            = openmc.Geometry()
    geometry.merge_surfaces             = True
    geometry.root_universe              = reactor_universe

    #######################################################################################################################################################
    ################################################################## Settings XML File ##################################################################
    #######################################################################################################################################################
    settings                            = openmc.Settings()
    settings.photon_transport           = True
    settings.batches                    = 200
    settings.inactive                   = 50
    settings.particles                  = 10000
    settings.generations_per_batch      = 3
    settings.output                     = {'tallies': False}
    settings.temperature                = {'method': 'interpolation', 'range': (250, 2500)}
    bounds                              = [-(91 + t_radial_reflector), -(91 + t_radial_reflector), -(L_HP/2 + t_axial_reflector), 
                                           91 + t_radial_reflector, 91 + t_radial_reflector, L_HP/2 + t_axial_reflector]
    uniform_dist                        = openmc.stats.Box(bounds[:3], bounds[3:])
    settings.source                     = openmc.IndependentSource(space = uniform_dist, constraints = {'domains': reactor_cells})
    if shannon_entropy == False:
        pass
    else:
        entropy_mesh                    = openmc.RegularMesh()
        entropy_mesh.lower_left         = [-(91 + t_radial_reflector), -(91 + t_radial_reflector), -(L_HP/2 + t_axial_reflector)]    
        entropy_mesh.upper_right        = [ 91 + t_radial_reflector, 91 + t_radial_reflector, L_HP/2 + t_axial_reflector]
        entropy_mesh.dimension          = (50, 50, 75)
        settings.entropy_mesh           = entropy_mesh
    #####################################################################################################################################################
    ################################################################### Making Plots ####################################################################
    #####################################################################################################################################################
    # plt.figure(figsize = (5.78, 5.78))
    # geometry.plot(basis = "xy", pixels = [1900, 1900], width = (182, 182), 
    #                 origin = (0, 0, 0), 
    #                 color_by = "material", 
    #                 colors = {UROX: "saddlebrown", helium: "white", zircaloy: "silver", 
    #                           graphite: "black", SS316: "darkorange", sodium: "mediumblue",
    #                           BeO: "lime"})
    # plt.xlabel('x [cm]')
    # plt.ylabel('y [cm]')
    # plt.savefig(fname = "unreflected_core", bbox_inches = "tight")

    #####################################################################################################################################################
    ##################################################################### Tallies #######################################################################
    #####################################################################################################################################################
    tallies                             = openmc.Tallies() 
    #####################################################################################################################################################
    particle_filter                     = openmc.ParticleFilter('neutron')
    #####################################################################################################################################################
    r_grid                              = np.linspace(0, 91 + t_radial_reflector, 100, endpoint = True)
    z_grid                              = np.linspace(-(L_HP/2 + t_axial_reflector), L_HP/2 + t_axial_reflector, 200, endpoint = True)
    cylindrical_mesh_r                  = openmc.CylindricalMesh(r_grid = r_grid, 
                                                                 z_grid = np.array([-(L_HP/2 + t_axial_reflector), L_HP/2 + t_axial_reflector]))
    cylindrical_mesh_z                  = openmc.CylindricalMesh(r_grid = np.array([0, 91 + t_radial_reflector]),
                                                                 z_grid = z_grid)
    mesh_filter_r                       = openmc.MeshFilter(cylindrical_mesh_r)
    mesh_filter_z                       = openmc.MeshFilter(cylindrical_mesh_z)
    #####################################################################################################################################################
    r_flux_tally                        = openmc.Tally(name = "Radial Flux")
    r_flux_tally.filters                = [mesh_filter_r, particle_filter]
    r_flux_tally.scores                 = ["flux"]
    tallies.append(r_flux_tally)
    #####################################################################################################################################################
    z_flux_tally                        = openmc.Tally(name = "Axial Flux")
    z_flux_tally.filters                = [mesh_filter_z, particle_filter]
    z_flux_tally.scores                 = ["flux"]
    tallies.append(z_flux_tally)
    #####################################################################################################################################################
    tally_heating                       = openmc.Tally(name = "Heating")
    tally_heating.scores                = ["heating"]
    tallies.append(tally_heating)
    #####################################################################################################################################################
    energy_range                        = np.logspace(-5, np.log10(20e6), num = 1001, base = 10.0)
    energy_filter                       = openmc.EnergyFilter(energy_range)
    #####################################################################################################################################################
    energy_spectra                      = openmc.Tally(name = "Flux")
    energy_spectra.filters              = [energy_filter, particle_filter]
    energy_spectra.scores               = ["flux"]
    tallies.append(energy_spectra)

    #####################################################################################################################################################
    ################################################################### Making Model ####################################################################
    #####################################################################################################################################################
    model = openmc.model.Model(geometry = geometry, materials = materials, settings = settings, tallies = tallies)
    model.export_to_xml()


###########################################################################################################################################################
################################################################# General Model Parameters ################################################################
###########################################################################################################################################################
enrichment                              = 0.20
###########################################################################################################################################################
r_HP                                    = 0.85
t_HP_clad                               = 0.35
L_HP                                    = 175
###########################################################################################################################################################
r_fuel                                  = 0.60
t_fuel_gap                              = 0.01
t_fuel_clad                             = 0.09
###########################################################################################################################################################
pitch                                   = 2.4
###########################################################################################################################################################
T_fuel                                  = 1000
T_gap                                   = 1000
T_cladding                              = 1000
T_graphite                              = 1000
T_HP                                    = 900
T_reflector                             = 1000
T_array                                 = [T_fuel, T_gap, T_cladding, T_graphite, T_HP, T_reflector]
###########################################################################################################################################################
reflectors                              = ["BeO", "graphite"]
axial_array                             = np.arange(0, 57.5, 2.5)
radial_array                            = np.arange(0, 57.5, 2.5)
t_axial_inf                             = 750
t_radial_inf                            = 750
###########################################################################################################################################################
OMC_threads                             = 14 # number of OpenMP threads to use for OpenMC calculation


###########################################################################################################################################################
################################################################## Spectra Calculation ####################################################################
###########################################################################################################################################################
if __name__ == "__main__":
    "StatePoint files were renamed manually since there were only three simulations"
    ref_mat = "BeO"
    t_axial = 0
    t_radial = 20 # changed to 0 for NR statepoint file
    make_model(r_fuel, t_fuel_gap, t_fuel_clad, r_HP, t_HP_clad, L_HP, pitch, T_array, ref_mat, enrichment, t_axial, t_radial)
    openmc.run(threads = OMC_threads)

###########################################################################################################################################################
############################################################### Axial Reflector Simulation ################################################################
###########################################################################################################################################################
if __name__ == "__main__":
    data                                = np.zeros((len(axial_array) + 1, len(reflectors) + 1))
    sigma_data                          = np.zeros((len(axial_array) + 1, len(reflectors) + 1))
    data[1:, 0]                         = axial_array
    sigma_data[1:, 0]                   = axial_array
    #######################################################################################################################################################
    for j, ref_mat in enumerate(reflectors):
        for i, t_axial in enumerate(axial_array):
            make_model(r_fuel, t_fuel_gap, t_fuel_clad, r_HP, t_HP_clad, L_HP, pitch, T_array, ref_mat, enrichment, t_axial, 0)
            print("\n"+f'{ref_mat} axial reflector, {t_axial} cm'+"\n")
            openmc.run(threads = OMC_threads)
            with openmc.StatePoint('statepoint.200.h5') as sp:
                k                       = sp.keff
                data[i+1, j+1]          = k.nominal_value
                sigma_data[i+1, j+1]    = k.std_dev
    #######################################################################################################################################################
    np.savetxt('axial'+'_k_data', data)
    np.savetxt('axial'+'_dk_data', sigma_data)

###########################################################################################################################################################
############################################################### Radial Reflector Simulation ###############################################################
###########################################################################################################################################################
if __name__ == "__main__":
    data                                = np.zeros((len(radial_array) + 1, len(reflectors) + 1))
    sigma_data                          = np.zeros((len(radial_array) + 1, len(reflectors) + 1))
    data[1:, 0]                         = radial_array
    sigma_data[1:, 0]                   = radial_array
    #######################################################################################################################################################
    for j, ref_mat in enumerate(reflectors):
        for i, t_rad in enumerate(radial_array):
            make_model(r_fuel, t_fuel_gap, t_fuel_clad, r_HP, t_HP_clad, L_HP, pitch, T_array, ref_mat, enrichment, 0, t_rad)
            print("\n"+f'{ref_mat} radial reflector, {t_rad} cm'+"\n")
            openmc.run(threads = OMC_threads)
            with openmc.StatePoint('statepoint.200.h5') as sp:
                k                       = sp.keff
                data[i+1, j+1]          = k.nominal_value
                sigma_data[i+1, j+1]    = k.std_dev
    #######################################################################################################################################################
    np.savetxt('radial'+'_k_data', data)
    np.savetxt('radial'+'_dk_data', sigma_data)

###########################################################################################################################################################
############################################################### "Inf" Reflector Simulation ################################################################
###########################################################################################################################################################
if __name__ == "__main__":
    data                                = np.zeros(len(reflectors))
    sigma_data                          = np.zeros(len(reflectors))
    ######################################################################################################################################################
    for j, ref_mat in enumerate(reflectors):
        make_model(r_fuel, t_fuel_gap, t_fuel_clad, r_HP, t_HP_clad, L_HP, pitch, T_array, ref_mat, enrichment, t_axial_inf, 0)
        print("\n"+f'{ref_mat} reflector'+"\n")
        openmc.run(threads = OMC_threads)
        with openmc.StatePoint('statepoint.200.h5') as sp:
            k                           = sp.keff
            data[j]                     = k.nominal_value
            sigma_data[j]               = k.std_dev
    #######################################################################################################################################################
    np.savetxt('inf_axial'+'k_data', data)
    np.savetxt('inf_axial'+'dk_data', sigma_data)
    #######################################################################################################################################################
    data                                = np.zeros(len(reflectors))
    sigma_data                          = np.zeros(len(reflectors))
    ######################################################################################################################################################
    for j, ref_mat in enumerate(reflectors):
        make_model(r_fuel, t_fuel_gap, t_fuel_clad, r_HP, t_HP_clad, L_HP, pitch, T_array, ref_mat, enrichment, 0, t_radial_inf)
        print("\n"+f'{ref_mat} reflector'+"\n")
        openmc.run(threads = OMC_threads)
        with openmc.StatePoint('statepoint.200.h5') as sp:
            k                           = sp.keff
            data[j]                     = k.nominal_value
            sigma_data[j]               = k.std_dev
    #######################################################################################################################################################
    np.savetxt('inf_radial'+'_k_data', data)
    np.savetxt('inf_radial'+'_dk_data', sigma_data)