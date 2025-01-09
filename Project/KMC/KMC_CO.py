import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np 

data = pd.ExcelFile('~/Desktop/TCC_SM_/Project/KMC/KMC_CO.xlsx')
KMC_CO_data_sheet = data.parse('Sheet1') 
pressure_column = KMC_CO_data_sheet.columns[0]
coverage_column = KMC_CO_data_sheet.columns[3]

# Selecting rows 6-14 (inclusive)
filtered_KMC_CO_data = KMC_CO_data_sheet.iloc[4:15]  # Python uses zero-based indexing
filtered_KMC_CO_data = filtered_KMC_CO_data[[pressure_column, coverage_column]]

# Rename columns
filtered_KMC_CO_data.columns = ['p', 'Coverage per site']

# Convert columns to numeric
filtered_KMC_CO_data['p'] = pd.to_numeric(filtered_KMC_CO_data['p'], errors='coerce')
filtered_KMC_CO_data['Coverage per site'] = pd.to_numeric(filtered_KMC_CO_data['Coverage per site'], errors='coerce')

# Drop rows with NaN values
filtered_KMC_CO_data = filtered_KMC_CO_data.dropna()

delta_g = 0.106  # eV 
T = 300  # K
boltzmann = 8.6173e-05  # eV / K
P = np.logspace(-2, 6, num =500)

constant = -delta_g / (T * boltzmann)
K = np.exp(constant)
theta = (K * P) / (1 + (K * P))

plt.figure(figsize=(10, 8))
plt.plot(P, theta, label='Langmuir isotherm', color='red')
plt.plot(filtered_KMC_CO_data['p'], filtered_KMC_CO_data['Coverage per site'], label='KMC Data', color='red', marker='o', linestyle="")

plt.xlabel('Pressure (bar)', fontsize=14, weight='bold')
plt.ylabel('Surface coverage', fontsize=14, weight='bold')
plt.title('Surface coverage vs pressure', fontsize=16, weight='bold')

# Custom x-axis scale and ticks
plt.xscale("log")

# Manually specify major ticks
major_ticks = [1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6]
plt.xticks(major_ticks, [f"$10^{{{int(np.log10(tick))}}}$" for tick in major_ticks])  # Custom tick labels

# Grid and layout
plt.grid(True, which="both", linestyle="--", linewidth=0.5)  # Grid for both major and minor ticks
plt.legend(loc='upper left', fontsize=8)
plt.tight_layout()
plt.show()


