import numpy as np 
import matplotlib.pyplot as plt 
import random 

# Set random seeds for reproducibility reasons
random.seed(13)
np.random.seed(13)

# Define constants 
T = 300                     # temperature in K
G_ads_CO = 0.106            # CO adsorption free energy
G_ads_CH4 = 0.061           # CH4 adsorption free energy
P_CO = 5.0                  # Partial pressure CO in bar
P_CH4 = 5.0                 # Partial pressure CH4 in bar 
P0 = 1.0                    # Reference pressure in bar
kB = 8.617333262145e-05     # Boltzmann constant (ev / K)
kBT = kB * T                # Thermal energy

# Determine adsorption probabilities 
P_ads_CO = min(1, (P_CO / P0) * np.exp(-G_ads_CO / kBT))
P_ads_CH4 = min(1, (P_CH4 / P0) * np.exp(-G_ads_CH4 / kBT))

# Determine the desorption probabilities 
P_des_CO = min(1, np.exp(G_ads_CO / kBT))
P_des_CH4 = min(1, np.exp(G_ads_CH4 / kBT))

# Determine the exchange probabilities 
# Exchange from CO to CH4
delta_G_CO_to_CH4 = G_ads_CH4 - G_ads_CO
P_exchange_CO_CH4 = min(1, np.exp(-delta_G_CO_to_CH4 / kBT))

# Exchange from CH4 to CO 
delta_G_CH4_to_CO = G_ads_CO - G_ads_CH4
P_exchange_CH4_CO = min(1, np.exp(-delta_G_CH4_to_CO / kBT))

# Print all the probabilities for debugging purposes 
print(f"Adsorption of CO (P_adsorb_CO): {P_ads_CO:.4f}")
print(f"Adsorption of CH4 (P_adsorb_CH4): {P_ads_CH4:.4f}", "\n")
print(f"Desorption of CO (P_desorb_CO): {P_des_CO:.4f}")
print(f"Desorption of CH4 (P_desorb_CH4): {P_des_CH4:.4f}", "\n")
print(f"Exchange from CO to CH4 (P_exchange_CO_to_CH4): {P_exchange_CO_CH4:.4f}")
print(f"Exchange from CH4 to CO (P_exchange_CH4_to_CO): {P_exchange_CH4_CO:.4f}", "\n")

# Build the lattice 
L = 10
lattice = np.zeros((L, L), dtype=int)

# Set up steps and storage of data 
N_steps = 10**6
N_snapshot = 10**3
equilibration_steps = N_steps // 2      # Equilibration will take place during the first half of the total steps 
production_steps = N_steps // 2         # Production of the actual lattice will be last half of the total steps

snapshots = []                  # Store snapshots of the lattice 
CO_count = []                   # Count the CO molecules during the production stage
CH4_count = []                  # Count the CH4 molecules during the production stage
selectivity_list = []

for step in range(1, N_steps + 1):
    # first select a random site to work on 
    random_row = np.random.randint(0, L)
    random_column = np.random.randint(0, L)
    current_state = lattice[random_row, random_column]

    # determine the possible actions based on the current state 
    # empty site; decide between CO or CH4 adsorption
    if current_state == 0:                                      # if the lattice site is empty, randomly generate between CO or CH4 adsorption
        action = random.choices(['adsorb_CO', 'adsorb_CH4'], weights=[0.5, 0.5], k=1)[0]

        if action == 'adsorb_CO':
            rand_val = random.random()                         
            if rand_val < P_ads_CO:                             # if random generated number is less than CO adsorption probability, change lattice site to 1
                lattice[random_row, random_column] = 1
        elif action == 'adsorb_CH4':
            rand_val = random.random()
            if rand_val < P_ads_CH4:
                lattice[random_row, random_column] = 2          # if random generated number is less than CH4 adsorption probability, change lattice site to 2

    # site occupied by CO; decide between CO desorption or exchange with CH4
    elif current_state == 1:
        action = random.choices(['desorb_CO', 'exchange_to_CH4'], weights=[0.5, 0.5], k=1)[0]

        if action == 'desorb_CO':
            rand_val = random.random()
            if rand_val < P_des_CO:
                lattice[random_row, random_column] = 0          # if random generated number is less than CO desorption probability, go back to an empty lattice
        elif action == 'exchange_to_CH4': 
            rand_val = random.random()
            if rand_val < P_exchange_CO_CH4:
                lattice[random_row, random_column] = 2          # if random generated number is less than probability of exchanging CO with CH4, CH4 is adsorped
    
    # site occupied by CH4; decide between CH4 desorption or exchange with CO
    elif current_state == 2:
        action = random.choices(['desorb_CH4', 'exchange_to_CO'], weights=[0.5, 0.5], k=1)[0]

        if action == 'desorb_CH4':
            rand_val = random.random()
            if rand_val < P_des_CH4:
                lattice[random_row, random_column] = 0          # if random generated number is less than CH4 desorption probability, go back to an empty lattice
        elif action == 'exchange_to_CO':
            rand_val = random.random()
            if rand_val < P_exchange_CH4_CO:
                lattice[random_row, random_column] = 1          # if random generated number is less than probability of exchanging CH4 with CO, CO is adsorped


    # take snapshots at the specified intervals 
    if step % N_snapshot == 0:
        snapshots.append(lattice.copy())

    # if entering the second half of the Monte Carlo steps, store the collected data from the production phase
    if step > equilibration_steps:
        CO = np.sum(lattice == 1)       # sum all the CO molecules on the lattice 
        CH4 = np.sum(lattice == 2)      # sum all the CH4 molecules on the lattice
        CO_count.append(CO)
        CH4_count.append(CH4)

        if CH4 > 0:
            selectivity = CO / CH4
        else:
            selectivity = np.inf  # Avoid division by zero
        
        selectivity_list.append(selectivity)

# Now calculate the average surface coverages 
avg_CO = np.mean(CO_count)
avg_CH4 = np.mean(CH4_count) 

# Calculate selectivity based on average surface coverage 
finite_selectivity = [s for s in selectivity_list if np.isfinite(s)]
if finite_selectivity: 
    avg_selectivity = np.mean(finite_selectivity)
else: 
    avg_selectivity = np.inf

# print the CO and CH4 surface coverage and the selectivity
print(f"Average CO surface coverage: {avg_CO:.4f}")
print(f"Average CH4 surface coverage: {avg_CH4:.4f}")
print(f"Selectivity (CO/CH4): {avg_selectivity:.4f}")