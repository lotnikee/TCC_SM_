import matplotlib.pyplot as plt 
import numpy as np 
import pandas as pd

# Define constants 
R = 8.314                                # J/mol
conversion_factor = 96485                # Add a conversion factor to easily convert to eV/K units 
h = 6.626e-34                            # J/s
c = 2.998e10                             # cm/s
v = 2160                                 # cm^-1
k = 1.381e-23                            # J/K

# Calculate the vibrational temperature vib_temp
vib_temp = (h * c * v) / k

# Define the translational and rotational contribution to contant volume heat capacity 
Cv_trans = (3/2) * R
Cv_rot = R 

# Calculating the theoretical values for the plot 
theory_trans_rot = 10000 * (5/2) * R / conversion_factor         # Multiply by E04 to get all values in order E-04 for a neater plot
theory_trans_rot_vib = 10000 * (7/2) * R / conversion_factor     # Multiply by E04 to get all values in order E-04 for a neater plot

# Add a range for the temperature (in K)
T = np.linspace(300, 6000, 1000)

# Calculating the vibrational contributions to C_v
vib_temp_T = vib_temp / T
vib_temp_T_square = vib_temp_T ** 2

exp_vib_temp_T = np.exp(vib_temp_T)
exp_vib_temp_T_adjusted = exp_vib_temp_T - 1
exp_vib_temp_T_squared = exp_vib_temp_T_adjusted ** 2
total_exp_vib_temp_T = exp_vib_temp_T / exp_vib_temp_T_squared

# Calculating the total C_v vibrational contribution
Cv_vib_total = R * vib_temp_T_square * total_exp_vib_temp_T

# Calculating total C_v constant 
Cv_total = 10000 * (Cv_trans + Cv_rot + Cv_vib_total) / conversion_factor

# Load the 'NO' datasheet and extract relevant columns 
data = pd.ExcelFile('~/Desktop/TCC_SM_/CW_2/CO/SM_CW2_CO.xlsx')
CO_data_sheet = data.parse('Sheet1') 
temperature_column = CO_data_sheet.columns[0]
cv_nist_column = CO_data_sheet.columns[2]

filtered_CO_data = CO_data_sheet[[temperature_column, cv_nist_column]].iloc[2:60]
filtered_CO_data.columns = ['T', 'C_v (NIST)']

# Plotting 
plt.figure(figsize=(10, 8))
# Plot theoretical horizontal lines
plt.axhline(theory_trans_rot, label='Theory (Trans + Rot: 5/2 * $k_B$)', color='black', linestyle='--')
plt.axhline(theory_trans_rot_vib, label='Theory (Trans + Rot + Vib: 7/2 * $k_B$)', color='blue', linestyle='--')
# Plot theoretical Cv_total vs Temperature
plt.plot(T, Cv_total, label='Theory (Trans + Rot + Vib)', color='blue')
# Plot NO data from the NIST Webbook 
plt.plot(filtered_CO_data['T'], filtered_CO_data['C_v (NIST)'], label='Data (NIST Webbook; Chase, 1998)', color='red')
# Plot design 
plt.xlabel('Temperature (K)', fontsize=14, weight='bold')
plt.ylabel('$C_v$ (eV/K) * E-04', fontsize=14, weight='bold')
plt.title('Constant Volume Heat Capacity $C_v$ of CO', fontsize=16, weight='bold')
plt.legend(loc='upper left', fontsize=8)
plt.grid(True)
plt.xlim(0, 6000)
plt.ylim(2, 3.2)
plt.tight_layout()
plt.show()
