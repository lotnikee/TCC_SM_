from copy import deepcopy
from ase.build import add_adsorbate, fcc111, molecule 
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
from ase.visualize import view
import numpy as np
import matplotlib.pyplot as plt

# Define constants 
temp = 300                                 
pressure = 1.0e+5    
kB = 8.61733326e-05                      
Pa_to_bar = 1.0e-5
T = np.linspace(100,1000)               # Temperature range needed for selectivity plot 

# Create a Cu(111) slab and get energy
slab = fcc111("Cu", size=(4,4,2), vacuum=10.0)
slab.calc = EMT()
energy_slab = slab.get_potential_energy()

# Get energy of CO molecule 
CO_mol = molecule("CO")
CO_mol.calc = EMT()
energy_CO_gas = CO_mol.get_potential_energy()

# Get energy of a CH4 molecule 
CH4_mol = molecule("CH4")
CH4_mol.calc = EMT()
energy_CH4_gas = CH4_mol.get_potential_energy()

# Run geometry optimisation of CO and CH4 on Cu(111) slab
CO_ads = deepcopy(slab)
add_adsorbate(slab=CO_ads, adsorbate=CO_mol, height=3.0, position=(3.82, 2.21))
constraint_CO = FixAtoms(mask=[atom.symbol == "Cu" for atom in CO_ads])
CO_ads.set_constraint(constraint_CO)
dyn = QuasiNewton(CO_ads, trajectory="CO_Cu(111).traj")
dyn.run(fmax=0.05)
energy_CO_ads = CO_ads.get_potential_energy()

CH4_ads = deepcopy(slab)
add_adsorbate(slab=CH4_ads, adsorbate=CH4_mol, height= 4.0, position=(3.82, 2.21))
constraint_CH4 = FixAtoms(mask=[atom.symbol == "Cu" for atom in CH4_ads])
CH4_ads.set_constraint(constraint_CH4)
dyn = QuasiNewton(CH4_ads, trajectory="CH4_Cu(111).traj")
dyn.run(fmax=0.05)
energy_CH4_ads = CH4_ads.get_potential_energy()

# Calculate and print adsoroption energies of CO and CH4
adsorption_energy_CO = energy_CO_ads - (energy_slab + energy_CO_gas)
adsorption_energy_CH4 = energy_CH4_ads - (energy_slab + energy_CH4_gas)
print(f"Adsorption energy of CO on Cu(111): {adsorption_energy_CO: .3f} eV")
print(f"Adsorption energy of CH4 on Cu(111): {adsorption_energy_CH4: .3f} eV")

# Calculate the adsorption free energy of CO 
vib_energy_CO_gas = [0.2634]
vib_energy_CO_ads = [0.2404, 0.0827, 0.0601, 0.0600, 0.0072, 0.0065]
thermo_CO_gas = IdealGasThermo(vib_energies=vib_energy_CO_gas,
                        geometry="linear",
                        potentialenergy=energy_CO_gas,
                        atoms=CO_mol,
                        symmetrynumber=1,
                        spin=0)
thermo_CO_ads = HarmonicThermo(vib_energies=vib_energy_CO_ads, potentialenergy=energy_CO_ads)

g_CO_gas = thermo_CO_gas.get_gibbs_energy(temperature=temp, pressure=pressure, verbose=False)
g_CO_ads = thermo_CO_ads.get_helmholtz_energy(temperature=temp, verbose=False)
g_slab_CO = energy_slab
adsorption_free_energy_CO = g_CO_ads - (g_slab_CO + g_CO_gas)

# Calculate the adsorption free energy of CH4
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

g_CH4_gas = thermo_CH4_gas.get_gibbs_energy(temperature=temp, pressure=pressure, verbose=False)
g_CH4_ads = thermo_CH4_ads.get_helmholtz_energy(temperature=temp, verbose=False)
g_slab = energy_slab
adsorption_free_energy_CH4 = g_CH4_ads - (g_slab + g_CH4_gas)

# Print the adsorption free energies of CO and Ch4
print(f"Adsorption free energy of CO on Cu(111) at {temp}K and {pressure*Pa_to_bar} bar: {adsorption_free_energy_CO: .3f} eV")
print(f"Adsorption free energy of CH4 on Cu(111) at {temp}K and {pressure*Pa_to_bar} bar: {adsorption_free_energy_CH4: .3f} eV")


# Calculate the equilibrium constant for adsorption K_CO
beta = kB * T 
exp_CO = -adsorption_free_energy_CO / beta
K_CO = np.exp(exp_CO)

# Calculate the equilibrium constant for adsorption K_CH4
exp_CH4 = -adsorption_free_energy_CH4 / beta
K_CH4 = np.exp(exp_CH4)

# Determine the selectivity of CO over CH4
S_CO_CH4 = K_CO / K_CH4  

# Determine the selectivity at x = 300 K
x_target = 300 
y_target = np.interp(x_target, T, np.log10(S_CO_CH4))

# Print the selectivity of CO over CH4 at 300 K
print(f"Selectivity (linear scale) at {x_target} K: {10**y_target:.4f}")

# Plot the results side by side 
fig, axs = plt.subplots(1, 2, figsize=(12, 5))  # figsize=(width, height) in inches
plt.title("Selectivity of CO over CH$_4$ on a Cu(111) surface as a function of T", weight='bold')
plt.xlim(0,1000)
plt.grid(linestyle = "--", linewidth=0.5, alpha=0.7)

# First subplot
axs[0].plot(T, np.log10(S_CO_CH4))
axs[0].set_title("Selectivity CO over CH$_4$ (log scale)", weight='bold')
axs[0].set_xlabel("Temperature (K)", weight='bold')
axs[0].set_ylabel("$S_{CO/CH4}$ (log scale)", weight='bold')
axs[0].grid(True)

# Second subplot
axs[1].plot(T, S_CO_CH4)
axs[1].set_title("Selectivity CO over CH$_4$ (linear)", weight='bold')
axs[1].set_xlabel("Temperature (K)", weight='bold')
axs[1].set_ylabel("$S_{CO/CH4}$ (linear scale)", weight='bold')
axs[1].grid(True)

plt.tight_layout()
plt.show()