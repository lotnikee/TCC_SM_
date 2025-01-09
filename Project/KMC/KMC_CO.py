import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 

# Reading the Excel sheet with all the adsorption isotherms
data = pd.ExcelFile('~/Desktop/TCC_SM_/Project/KMC/KMC_CO.xlsx')
KMC_CO_data_sheet = data.parse('Sheet1') 
pressure_column = KMC_CO_data_sheet.columns[0]
coverage_column = KMC_CO_data_sheet.columns[3]

filtered_KMC_CO_data = KMC_CO_data_sheet[[pressure_column, coverage_column]].iloc[15:27]
filtered_KMC_CO_data.columns = ['p', 'Coverage per site']

# Plotting
plt.figure(figsize=(10, 8))
plt.plot(filtered_KMC_CO_data['p'], filtered_KMC_CO_data['Coverage per site'], label='KMC Data', color='red', marker="o", linestyle="")
plt.xlabel('Pressure (Pa)', fontsize=14, weight='bold')
plt.ylabel('Surface coverage', fontsize=14, weight='bold')
plt.title('Surface coverage vs pressure', fontsize=16, weight='bold')

# Custom x-axis scale and ticks
plt.xscale("log")
plt.xlim(0.5e3, 1.5e11)  # Set x-axis limits

# Manually specify major ticks
major_ticks = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11]
plt.xticks(major_ticks, [f"$10^{{{int(np.log10(tick))}}}$" for tick in major_ticks])  # Custom tick labels

# Grid and layout
plt.grid(True, which="both", linestyle="--", linewidth=0.5)  # Grid for both major and minor ticks
plt.legend(loc='upper left', fontsize=8)
plt.ylim(-0.05, 1.05)
plt.tight_layout()
plt.show()
