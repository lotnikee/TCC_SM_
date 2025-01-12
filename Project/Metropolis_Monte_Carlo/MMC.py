import numpy as np 
import random

# Define some constants 
t = 300     # K
G_ads_CO = -0.434   # eV
G_ads_CH4 = -0.239  # eV
P_CO = 5    # bar
P_CH4 = 5   # bar
P_o = 1     # bar
kB = 8.61733326e-05                      

# Calculate the probability of adsorption for CO 
constant_CO = - G_ads_CO / (kB * t)
exp_CO = np.exp(constant_CO)
pres_constant_CO = (P_CO/P_o) * exp_CO
prob_ads_CO = min(1, pres_constant_CO)

# Calculate the probability of adsorption for CH4
constant_CH4 = - G_ads_CH4 / (kB * t)
exp_CH4 = np.exp(constant_CH4)
pres_constant_CH4 = (P_CH4/P_o) * exp_CH4
prob_ads_CH4 = min(1, pres_constant_CH4)

# Print adsorption expontentials for debugging purposes
print(f"The adsorption exponential constant for CO is: {pres_constant_CO: .3f}")
print(f"The adsorption exponential constant for CH4 is: {pres_constant_CH4: .3f}")

# Calculate the probability of desorption for CO and CH4
prob_des_CO = min(1, exp_CO)
prob_des_CH4 = min(1, exp_CH4)

# Print desorption expontentials for debugging purposes
print(f"The desorption expontential for CO is: {exp_CO: .3f}")
print(f"The desorption expontential for CH4 is: {exp_CH4: .3f}")


# Calculate the exchange probability 
delta_G = - G_ads_CH4 - G_ads_CO
constant = delta_G / (kB * t)
exp = np.exp(constant)
pres_constant = (P_CH4 / P_CO) * exp
prob_exchange = min(1, pres_constant)

# Print exchange probability 
print(f"The exchange expontential is: {pres_constant: .3f}")

# Build an LxL lattice 
L = np.zeros((3, 3), dtype=int)
print(L)

number = random.randint(1, 1e+8)

if 1<= number <= prob_ads_CH4:
    zero_indices = np.argwhere(L==0)
    
    if zero_indices > 0: 
        selected_index = tuple(zero_indices[random.randint(0, len(zero_indices)-1)])

        L[selected_index] = 2

elif prob_ads_CH4 <= number <= prob_ads_CO:
    zero_indices = np.argwhere(L==0)
    
    if zero_indices > 0: 
        selected_index = tuple(zero_indices[random.randint(0, len(zero_indices)-1)])

        L[selected_index] = 1

print(number)
print("Updated matrix:")
print(L)
