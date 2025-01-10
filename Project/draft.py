from copy import deepcopy
from ase.build import add_adsorbate, fcc111, molecule 
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
from ase.visualize import view
import numpy as np
import matplotlib.pyplot as plt



CH4_ads = deepcopy(slab)
add_adsorbate(slab=CH4_ads, adsorbate=CH4_mol, height= 4.0, position=(3.82, 2.21))
constraint_CH4 = FixAtoms(mask=[atom.symbol == "Cu" for atom in CH4_ads])
CH4_ads.set_constraint(constraint_CH4)
dyn = QuasiNewton(CH4_ads, trajectory="CH4_Cu(111).traj")
dyn.run(fmax=0.05)
energy_CH4_ads = CH4_ads.get_potential_energy()

adsorption_energy_CH4 = energy_CH4_ads - (energy_slab + energy_CH4_gas)
print("Adsorption energy", adsorption_energy_CH4)


vib_energy_CH4_gas = [0.3843, 0.3840, 0.3840, 0.3685, 0.1881, 0.1879,
                    0.1595, 0.1593, 0.1592]
vib_energy_CH4_ads = [0.3815, 0.3758, 0.3758, 0.3625, 0.1850, 0.1848,
                    0.1589, 0.1584, 0.1559, 0.0161, 0.0161, 0.0112, 0.0061, 0.0061, 0.0061]
thermo_CH4_gas = IdealGasThermo(vib_energies=vib_energy_CH4_gas,
                        geometry="nonlinear",
                        potentialenergy=energy_CH4_gas,
                        atoms=CH4_mol,
                        symmetrynumber=12,
                        spin=0)
thermo_CH4_ads = HarmonicThermo(vib_energies=vib_energy_CH4_ads, potentialenergy=energy_CH4_ads)
temp = 300
pressure = 1.0e+5

g_CH4_gas = thermo_CH4_gas.get_gibbs_energy(temperature=temp, pressure=pressure, verbose=False)
g_CH4_ads = thermo_CH4_ads.get_helmholtz_energy(temperature=temp, verbose=False)
g_slab = energy_slab
Pa_to_bar = 1.0e-5
adsorption_free_energy_CH4 = g_CH4_ads - (g_slab + g_CH4_gas)

print("Adsorption free energy of CH4:", adsorption_free_energy_CH4)