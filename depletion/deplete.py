###########################################################################################################################################################
#################################################################### Relevant Libraries ###################################################################
###########################################################################################################################################################
import openmc
import openmc.deplete
import numpy as np
import json

###########################################################################################################################################################
############################################################# Chain-file and XS-file locations ############################################################
###########################################################################################################################################################
'''
HepLab
'''
# material_XS_path                    = '/uio/hume/student-u63/baltasah/Data_Libraries/endfb8.0/cross_sections.xml'
# chain_file_path                     = '/uio/hume/student-u63/baltasah/Documents/chain_endfb80_pwr.xml'
# fissionq_path                       = '/uio/hume/student-u63/baltasah/Documents/giroscope/serpent_fissq.json'
'''
Laptop
'''
material_XS_path                    = '/home/baltasar/Documents/endfb-vii.1-hdf5/cross_sections.xml'
chain_file_path                     = '/home/baltasar/Documents/depletion_chain_files/chain_endfb71_pwr.xml'
fissionq_path                       = '/home/baltasar/Documents/Giroscope/serpent_fissq.json'
'''
Fox
'''
# material_XS_path                    = "/fp/projects01/ec12/ec-baltasah/cross_sections/endfb-viii.0-hdf5/cross_sections.xml"
# fissionq_path                       = "/fp/homes01/u01/ec-baltasah/scripts/serpent_fissq.json"
# chain_file_path                     = "/fp/projects01/ec12/ec-baltasah/chain_files/chain_endfb80_pwr.xml"

###########################################################################################################################################################
#################################################################### Reactor Parameters ###################################################################
###########################################################################################################################################################
r_fuel                              = 0.60
t_fuel_gap                          = 0.01
t_fuel_clad                         = 0.09
###########################################################################################################################################################
r_HP                                = 0.85
t_HP_clad                           = 0.35
L                                   = 175
pitch                               = 2.4
t_axial_reflector                   = 15
t_radial_reflector                  = 20
###########################################################################################################################################################
enrichment                          = 0.20
###########################################################################################################################################################
T_HP                                = 900
T_reflector                         = 900
###########################################################################################################################################################
temp_dictionary                     = {
                                        "5 MWth": {"Power": 5e6,
                                                    "Fuel": 924,
                                                    "Gap": 908,
                                                    "Clad": 902,
                                                    "Block": 901},
                                        "10 MWth": {"Power": 10e6,
                                                    "Fuel": 949,
                                                    "Gap": 917,
                                                    "Clad": 905,
                                                    "Block": 902},
                                        "15 MWth": {"Power": 15e6,
                                                    "Fuel": 974,
                                                    "Gap": 925,
                                                    "Clad": 907,
                                                    "Block": 903}
                                                    }
###########################################################################################################################################################
reactor_power                       = "15 MWth"
timesteps                           = [5, 10, 15, 20, 50, 
                                       100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 
                                       100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                                       100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
                                       100, 100, 100, 100, 100]

###########################################################################################################################################################
######################################################################## Materials ########################################################################
###########################################################################################################################################################
materials_list                      = []
###########################################################################################################################################################
UROX_inner  = openmc.Material(name = "UROX-inner")
UROX_inner.add_element("O", 2)
UROX_inner.add_nuclide("U238", 1 - enrichment)
UROX_inner.add_nuclide("U235", enrichment)
UROX_inner.add_s_alpha_beta("c_U_in_UO2")
UROX_inner.set_density("kg/m3", np.round(0.95 * 10960, 0))
UROX_inner.temperature = temp_dictionary[reactor_power]["Fuel"]
UROX_inner.volume = (np.pi * r_fuel**2 * L) * 7 * 144
materials_list.append(UROX_inner)
###########################################################################################################################################################
UROX_outer = openmc.Material(name = "UROX-outer")
UROX_outer.add_element("O", 2)
UROX_outer.add_nuclide("U238", 1 - enrichment)
UROX_outer.add_nuclide("U235", enrichment)
UROX_outer.add_s_alpha_beta("c_U_in_UO2")
UROX_outer.set_density("kg/m3", np.round(0.95 * 10960, 0))
UROX_outer.temperature = temp_dictionary[reactor_power]["Fuel"]
UROX_outer.volume = (np.pi * r_fuel**2 * L) * 12 * 144
materials_list.append(UROX_outer)
###########################################################################################################################################################
helium = openmc.Material(name = "Gap Helium")
helium.add_element("He", 1)
helium.set_density("kg/m3", 0.973)
helium.temperature = temp_dictionary[reactor_power]["Gap"]
materials_list.append(helium)
###########################################################################################################################################################
zircaloy = openmc.Material(name = "Zircaloy-4")
zircaloy.add_element("Sn", 1.45/100, "wo")
zircaloy.add_element("Fe", 0.21/100, "wo")
zircaloy.add_element("Cr", 0.10/100, "wo")
zircaloy.add_element("O", 0.125/100, "wo")
zircaloy.add_element("Zr", (100-1.45-0.21-0.10-0.125)/100, "wo")
zircaloy.set_density("kg/m3", 6501)
zircaloy.temperature = temp_dictionary[reactor_power]["Clad"]
materials_list.append(zircaloy)
###########################################################################################################################################################
graphite = openmc.Material(name = "G-348 Graphite")
graphite.add_element("C", 1)
graphite.set_density("kg/m3", 1800)
graphite.add_s_alpha_beta("c_Graphite")
graphite.temperature = temp_dictionary[reactor_power]["Block"]
materials_list.append(graphite)
###########################################################################################################################################################
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
SS316.temperature = T_HP
materials_list.append(SS316)
###########################################################################################################################################################
sodium = openmc.Material(name = "Sodium")
sodium.add_element("Na", 1)
sodium.set_density("kg/m3", np.round(0.80 * 1.70 * 10**-2 + 0.20*805, 2))
sodium.temperature = T_HP
materials_list.append(sodium)
###########################################################################################################################################################
BeO = openmc.Material(name = "BeO")
BeO.add_element("O", 1)
BeO.add_element("Be", 1)
BeO.set_density("kg/m3", 3005)
BeO.add_s_alpha_beta("c_Be_in_BeO")
BeO.temperature = T_reflector
materials_list.append(BeO)
###########################################################################################################################################################
materials                           = openmc.Materials(materials_list)

###########################################################################################################################################################
###################################################################### Make Functions #####################################################################
###########################################################################################################################################################
def make_fuel_pin(fuel_mat : openmc.Material):
    '''
    Function for creating openmc.Universe corresponding to fuel rod. \n
    The fuel pin consists of a fuel pin in Zircaloy 4 cladding, with a betryllium-oxide reflector pin at the top

    Parameters
    ----------
    fuel_mat : openmc.Material used to fill the fuel pin region, e.g. UROX

    Returns
    -------
    openmc.Universe corresponding to fuel rod with top reflector
    '''
    z_fuel                          = [openmc.ZPlane(z0 = z) for z in np.linspace(-L/2, L/2, 6, endpoint = True)]
    z_core_top                      = openmc.ZPlane(z0 = L/2 + t_axial_reflector,    boundary_type = "vacuum")
    fuel_cylinder                   = openmc.ZCylinder(r = r_fuel)
    gap_cylinder                    = openmc.ZCylinder(r = r_fuel + t_fuel_gap)
    clad_cylinder                   = openmc.ZCylinder(r = r_fuel + t_fuel_gap + t_fuel_clad)
    #######################################################################################################################################################
    gap_cell                        = openmc.Cell(fill = helium, 
                                              region = + fuel_cylinder &- gap_cylinder &+ z_fuel[0] &- z_fuel[-1])
    clad_cell                       = openmc.Cell(fill = zircaloy,
                                                region = + gap_cylinder &- clad_cylinder &+ z_fuel[0] &- z_fuel[-1])
    fuel_cells                      = [openmc.Cell(fill = fuel_mat,
                                                region = - fuel_cylinder &+ z_fuel[n-1] &- z_fuel[n]) 
                                                for n in range(1, len(z_fuel))]
    fuel_ref_cell                   = openmc.Cell(fill = BeO,
                                                region = + z_fuel[-1] &- z_core_top)
    fuel_outer_cell                 = openmc.Cell(fill = graphite,
                                                region = + clad_cylinder &+ z_fuel[0] &- z_fuel[-1])
    #######################################################################################################################################################
    return openmc.Universe(cells = fuel_cells + [gap_cell, clad_cell, fuel_outer_cell, fuel_ref_cell])

def make_HP():
    '''
    Function for creating openmc.Universe corresponding to a heat pipe \n
    Heat pipe consists of a sodium vapor core and stainless steel wick/cladding structure.

    Returns
    -------
    openmc.Universe corresponding to a heat pipe (evaporator section and adiabatic section in top axial reflector)

    '''
    z_core_top                      = openmc.ZPlane(z0 = L/2 + t_axial_reflector,    boundary_type = "vacuum")
    z_fuel                          = [openmc.ZPlane(z0 = z) for z in np.linspace(-L/2, L/2, 6, endpoint = True)]
    HP_cylinder                     = openmc.ZCylinder(r = r_HP)
    HP_clad_cylinder                = openmc.ZCylinder(r = r_HP + t_HP_clad)
    #######################################################################################################################################################
    HP_core_cell                    = openmc.Cell(fill = sodium,
                                              region = - HP_cylinder &+ z_fuel[0] &- z_core_top)
    HP_clad_cell                    = openmc.Cell(fill = SS316,
                                                region = + HP_cylinder &- HP_clad_cylinder &+ z_fuel[0] &- z_core_top)
    HP_outer_mod_cell               = openmc.Cell(fill = graphite,
                                                region = + HP_clad_cylinder &+ z_fuel[0] &- z_fuel[-1])
    HP_outer_ref_cell               = openmc.Cell(fill = BeO,
                                                region = + HP_clad_cylinder &+ z_fuel[-1] &- z_core_top)
    #######################################################################################################################################################
    return openmc.Universe(cells = [HP_core_cell, HP_clad_cell, HP_outer_mod_cell, HP_outer_ref_cell])

def make_assembly(fuel_mat : openmc.Material):
    '''
    Function for creating openmc.Universe representing a single fuel assembly. \n
    One fuel assembly consists of 144 fuel rods and 73 heat pipes in a graphite block.

    Parameters
    ----------
    fuel_mat : openmc.Material used for filling the fuel pin region
    Returns
    -------
    openmc.Universe corresponding to a single hexagonal fuel assembly
    '''
    HP                              = make_HP()
    fuel_pin                        = make_fuel_pin(fuel_mat)
    #######################################################################################################################################################
    z_fuel                          = [openmc.ZPlane(z0 = z) for z in np.linspace(-L/2, L/2, 6, endpoint = True)]
    z_core_top                      = openmc.ZPlane(z0 = L/2 + t_axial_reflector,    boundary_type = "vacuum")
    assembly_hex                    = openmc.model.HexagonalPrism(pitch*8.7, orientation = 'x')
    #######################################################################################################################################################
    assembly_outer_mod_cell         = openmc.Cell(fill = graphite,
                                              region = + z_fuel[0] &- z_fuel[-1])
    assembly_outer_ref_cell         = openmc.Cell(fill = BeO,
                                                region = + z_fuel[-1] &- z_core_top)
    assembly_outer_universe         = openmc.Universe(cells = [assembly_outer_mod_cell, assembly_outer_ref_cell])
    #######################################################################################################################################################
    assembly_lattice                = openmc.HexLattice()
    assembly_lattice.center         = (0, 0)
    assembly_lattice.pitch          = (pitch, )
    assembly_lattice.orientation    = 'x'
    assembly_lattice.outer          = assembly_outer_universe
    #######################################################################################################################################################
    ring_list                       = [[HP], [fuel_pin]*6]
    ring3                           = [HP if i%2 == 0 else fuel_pin for i in range(1,13)]
    ring_list.append(ring3)
    #######################################################################################################################################################
    ring4                           = [HP if i in range(1, 17, 3) else fuel_pin for i in range(1, 19)]
    ring_list.append(ring4)
    #######################################################################################################################################################
    ring5                           = [HP if i in range(3, 24, 4) else fuel_pin for i in range(1, 25)]
    ring_list.append(ring5)
    #######################################################################################################################################################
    ring6                           = [HP if i in [2, 5, 7, 10, 12, 15, 17, 20, 22, 25, 27, 30]
                                    else fuel_pin for i in range(1,31)]
    ring_list.append(ring6)
    #######################################################################################################################################################
    ring7                           = [HP if i in range(1, 35, 3) else fuel_pin for i in range(1, 37)]
    ring_list.append(ring7)
    #######################################################################################################################################################
    ring8                           = [HP if i in [3, 6, 10, 13, 17, 20, 24, 27, 31, 34, 38, 41]
                                    else fuel_pin for i in range(1, 43)]
    ring_list.append(ring8)
    #######################################################################################################################################################
    ring9                           = [HP if i in [2, 5, 8, 10, 13, 16, 18, 21, 24, 26, 29, 32, 34, 37, 40, 42, 45, 48]
                                    else fuel_pin for i in range(1, 49)]
    ring_list.append(ring9)
    #######################################################################################################################################################
    ring_list.reverse()
    assembly_lattice.universes      = ring_list
    #######################################################################################################################################################
    main_assembly_cell              = openmc.Cell(fill = assembly_lattice,
                                                region = - assembly_hex &+ z_fuel[0] &- z_core_top)
    outer_hex_ref_cell              = openmc.Cell(fill = BeO,
                                                region = + assembly_hex &+ z_fuel[-1] &- z_core_top)
    outer_hex_mod_cell              = openmc.Cell(fill = graphite,
                                                region = + assembly_hex &+ z_fuel[0] &- z_fuel[-1])
    #######################################################################################################################################################
    return openmc.Universe(cells = [main_assembly_cell, outer_hex_mod_cell, outer_hex_ref_cell])

def make_reactor():
    '''
    Function for creating the openmc.Universe corresponding to the full HPR core. \n
    The core consists of 19 fuel assemblies arranged in a hexagonal lattice and a beryllium-oxide reflector (radial and axial)

    Returns
    -------
    universe : openmc.Universe corresponding to the full reactor core
    reactor_cells : list containing several instances openmc.Cell used to create the reactor universe
    '''
    z_fuel                          = [openmc.ZPlane(z0 = z) for z in np.linspace(-L/2, L/2, 6, endpoint = True)]
    z_core_bottom                   = openmc.ZPlane(z0 = -L/2 - t_axial_reflector,   boundary_type = "vacuum")
    z_core_top                      = openmc.ZPlane(z0 = L/2 + t_axial_reflector,    boundary_type = "vacuum")
    core_cylinder                   = openmc.ZCylinder(r = 91)
    reflector_cylinder              = openmc.ZCylinder(r = 91 + t_radial_reflector,    boundary_type = "vacuum") 
    #######################################################################################################################################################
    assembly_1                      = make_assembly(UROX_inner)
    assembly_2                      = make_assembly(UROX_outer)
    #######################################################################################################################################################
    reactor_lattice_mod_cell        = openmc.Cell(fill = graphite,
                                                region = + z_fuel[0] &- z_fuel[-1])
    reactor_lattice_top_ref_cell    = openmc.Cell(fill = BeO,
                                                region = + z_fuel[-1] &- z_core_top)
    reactor_outer_universe          = openmc.Universe(cells = [reactor_lattice_mod_cell, reactor_lattice_top_ref_cell])
    #######################################################################################################################################################
    reactor_lattice                 = openmc.HexLattice()
    reactor_lattice.center          = (0, 0)
    reactor_lattice.pitch           = (pitch*8.75*np.sqrt(3), )
    reactor_lattice.orientation     = 'y'
    reactor_lattice.outer           = reactor_outer_universe
    reactor_lattice.universes       = [[assembly_2]*12, [assembly_1]*6, [assembly_1]]
    #######################################################################################################################################################
    core_cell                       = openmc.Cell(fill = reactor_lattice,
                                                region = - core_cylinder &+ z_fuel[0] &- z_core_top)
    core_lower_reflector_cell       = openmc.Cell(fill = BeO,
                                                region = - reflector_cylinder &+ z_core_bottom &- z_fuel[0])
    core_reflector_cell             = openmc.Cell(fill = BeO,
                                                region = + core_cylinder &- reflector_cylinder &+ z_fuel[0] &- z_core_top)
    #######################################################################################################################################################
    reactor_cells                   = [core_cell, core_lower_reflector_cell, core_reflector_cell]
    return openmc.Universe(name = "Reactor Core", cells = reactor_cells), reactor_cells


###########################################################################################################################################################
######################################################################## Geometry #########################################################################
###########################################################################################################################################################
reactor_universe, cells         = make_reactor()
geometry                        = openmc.Geometry()
geometry.merge_surfaces         = True
geometry.root_universe          = reactor_universe

###########################################################################################################################################################
######################################################################## Settings #########################################################################
###########################################################################################################################################################
settings                        = openmc.Settings()
###########################################################################################################################################################
settings.photon_transport       = True
settings.batches                = 250
settings.inactive               = 50
settings.particles              = 25000
settings.generations_per_batch  = 4
###########################################################################################################################################################
settings.output                 = {'tallies': False}
settings.temperature            = {'method': 'interpolation', 'range': (250, 2500)}
###########################################################################################################################################################
bounds                          = [-(91 + t_radial_reflector), -(91 + t_radial_reflector), -(L/2 + t_axial_reflector), 
                                    91 + t_radial_reflector, 91 + t_radial_reflector, L/2 + t_axial_reflector]
uniform_dist                    = openmc.stats.Box(bounds[:3], bounds[3:])
settings.source                 = openmc.IndependentSource(space = uniform_dist, constraints = {'domains': cells, 'fissionable': True})

###########################################################################################################################################################
######################################################################## Settings #########################################################################
###########################################################################################################################################################
tallies                         = openmc.Tallies() 
#####################################################################################################################################################
particle_filter                 = openmc.ParticleFilter('neutron')
#####################################################################################################################################################
r_grid                          = np.linspace(0, 91 + t_radial_reflector, 100, endpoint = True)
z_grid                          = np.linspace(-(L/2 + t_axial_reflector), L/2 + t_axial_reflector, 200, endpoint = True)
cylindrical_mesh_r              = openmc.CylindricalMesh(r_grid = r_grid, 
                                                        z_grid = np.array([-(L/2 + t_axial_reflector), L/2 + t_axial_reflector]))
cylindrical_mesh_z              = openmc.CylindricalMesh(r_grid = np.array([0, 91 + t_radial_reflector]), 
                                                        z_grid = z_grid)
mesh_filter_r                   = openmc.MeshFilter(cylindrical_mesh_r)
mesh_filter_z                   = openmc.MeshFilter(cylindrical_mesh_z)
#####################################################################################################################################################
r_flux_tally                    = openmc.Tally(name = "Radial Flux")
r_flux_tally.filters            = [mesh_filter_r, particle_filter]
r_flux_tally.scores             = ["flux"]
tallies.append(r_flux_tally)
#####################################################################################################################################################
z_flux_tally                    = openmc.Tally(name = "Axial Flux")
z_flux_tally.filters            = [mesh_filter_z, particle_filter]
z_flux_tally.scores             = ["flux"]
tallies.append(z_flux_tally)
#####################################################################################################################################################
tally_heating                   = openmc.Tally(name = "Heating")
tally_heating.scores            = ["heating"]
tallies.append(tally_heating)
#####################################################################################################################################################
energy_range                    = np.logspace(-5, np.log10(20e6), num = 1001, base = 10.0)
energy_filter                   = openmc.EnergyFilter(energy_range)
#####################################################################################################################################################
energy_spectra                  = openmc.Tally(name = "Flux")
energy_spectra.filters          = [energy_filter, particle_filter]
energy_spectra.scores           = ["flux"]
tallies.append(energy_spectra)


###########################################################################################################################################################
###################################################################### Making Model #######################################################################
###########################################################################################################################################################
model                           = openmc.model.Model(geometry = geometry, materials = materials, settings = settings, tallies = tallies)
model.export_to_model_xml()

###########################################################################################################################################################
#################################################################### Depletion Settings ###################################################################
###########################################################################################################################################################
###########################################################################################################################################################
with open(fissionq_path, 'rb') as fq:
    fission_q                   = json.load(fq)
###########################################################################################################################################################
operator                        = openmc.deplete.CoupledOperator(model, 
                                                                 chain_file = chain_file_path,
                                                                 normalization_mode = 'fission-q',
                                                                 fission_yield_mode = 'average',
                                                                 fission_q = fission_q)
###########################################################################################################################################################
integrator                      = openmc.deplete.CECMIntegrator(operator, 
                                                                timesteps,
                                                                timestep_units = 'd',
                                                                power = temp_dictionary[reactor_power]["Power"],
                                                                solver = 'cram48')
###########################################################################################################################################################
if __name__ == "__main__":
    integrator.integrate()