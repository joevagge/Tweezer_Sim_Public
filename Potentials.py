#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 11:43:06 2024

@author: joe
"""

import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from Trap_Simulation import *

# Sample R values
R = np.logspace(-9,  -6, 5000)

# Add blue lines for each detuning
colors = ['darkgreen', 'darkgreen', 'darkgreen']
linestyles = ['solid', 'dashed', 'dotted']

Detuning_1 = 50e6
Detuning_2 = 50e6
Detuning_3 = 50e6

Laser_Freq_1 = D1_Freq + Detuning_1
Laser_Freq_2 = D1_Freq + Detuning_2
Laser_Freq_3 = D1_Freq + Detuning_3

# Laser frequencies
Laser_Frequencies = [Laser_Freq1, Laser_Freq2, Laser_Freq3]

# Calculate results
results = [Collision_Separation_Calculator(R, freq, Delta) for freq in Laser_Frequencies]

# Calculate potentials
SP_potential = SPPotential(Delta, CouplingPotential(R))
PS_potential = PSPotential(Delta, CouplingPotential(R))
SS_potential = SSPotential(Delta, CouplingPotential(R))

# Add a vertical line at the flattening point
flattening_point = find_flattening_point(R, SP_potential, 0.001 * Detuning1)
# print(f"Asymptote at {flattening_point}")


####### Plotting PECs with broken y-axis ###########
chosen_linewidth = 3
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(10, 11), dpi=500) 
                                # ,gridspec_kw={'height_ratios': [2, 1]})  # Adjust the height ratios here

# Upper subplot
ax1.plot(R * 1e6, SP_potential * 1e-12, color='blue',linewidth = chosen_linewidth)
ax1.plot(R * 1e6, PS_potential * 1e-12, color='red',linewidth = chosen_linewidth)
# ax1.set_ylabel('Potential (THz)', fontsize=20)
# ax1.legend(loc='upper right')

# Lower subplot
# ax2.plot(R * 1e6, SP_potential * 1e-12, label=f"SP", color='blue',linewidth = chosen_linewidth)
# ax2.plot(R * 1e6, PS_potential * 1e-12, label=f"PS", color='red',linewidth = chosen_linewidth)
# ax2.plot(R * 1e6, SS_potential * 1e-12, label=r"SS", color='orange',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, SP_potential * 1e-12, label=f"SP $= \Delta + U$", color='blue',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, PS_potential * 1e-12, label=f"PS $= \Delta - U$", color='red',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, SS_potential * 1e-12, label=r"SS $= \Delta - \sqrt{\Delta^2 + U^2}$", color='orange',linewidth = chosen_linewidth)
ax2.set_xlabel('Distance (µm)', fontsize=20)
ax2.set_ylabel('Potential (THz)', fontsize=20, loc="top")
ax2.yaxis.set_label_coords(-0.16, 1.23)
# Set limits
ax1.set_ylim([(1 - 3e-8) * (SP_potential[-1] * 1e-12), (1 + 3e-8) * (SP_potential[-1] * 1e-12)])  # Set limits for the upper part
ax1.set_xlim([6e-2, 1])  # Set x-axis limits in micrometers
ax2.set_ylim([-0.01, 0.04])  # Set limits for the lower part

# Add lines for each detuning
for i, (min_R, sp_potential_at_min_R, ps_potential_at_min_R, ss_potential_at_min_R) in enumerate(results):
    detuning = Laser_Frequencies[i] - Delta
    color = colors[i % len(colors)]  # Cycle through colors
    linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
    
    add_line_to_plot(ax1, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, chosen_linewidth)
    add_line_to_plot(ax2, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, chosen_linewidth)
    
# Add a vertical line at the flattening point
# ax1.axvline(flattening_point * 1e6, color='purple', linestyle='--', label=f'Flattening Point at R = {flattening_point*1e6:.2} µm')
# ax2.axvline(flattening_point * 1e6, color='purple', linestyle='--', label=f'Flattening Point at R = {flattening_point*1e6:.2} µm')

ax1.tick_params(axis='both', labelsize=16)
ax2.tick_params(axis='both', labelsize=16)
ax2.legend(loc='upper right', fontsize=16)

ax1.get_yaxis().get_major_formatter().set_useOffset(False)

ax1.set_xscale('log')
ax2.set_xscale('log')

# plt.suptitle('Potentials as a Function of Distance (R)', fontsize=20)
plt.show()

# Upper subplot
ax1.plot(R * 1e6, SP_potential * 1e-12, color='blue',linewidth = chosen_linewidth)
ax1.plot(R * 1e6, PS_potential * 1e-12, color='red',linewidth = chosen_linewidth)
# ax1.set_ylabel('Potential (THz)', fontsize=20)
# ax1.legend(loc='upper right')

# Lower subplot
ax2.plot(R * 1e6, SP_potential * 1e-12, label=f"SP", color='blue',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, PS_potential * 1e-12, label=f"PS", color='red',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, SS_potential * 1e-12, label=r"SS", color='orange',linewidth = chosen_linewidth)
# ax2.plot(R * 1e6, SP_potential * 1e-12, label=f"SP $= \Delta + U$", color='blue',linewidth = chosen_linewidth)
# ax2.plot(R * 1e6, PS_potential * 1e-12, label=f"PS $= \Delta - U$", color='red',linewidth = chosen_linewidth)
# ax2.plot(R * 1e6, SS_potential * 1e-12, label=r"SS $= \Delta - \sqrt{\Delta^2 + U^2}$", color='orange',linewidth = chosen_linewidth)
ax2.set_xlabel('Distance (µm)', fontsize=20)
ax2.set_ylabel('Potential (THz)', fontsize=20, loc="top")
ax2.yaxis.set_label_coords(-0.16, 1.23)
# Set limits
ax1.set_ylim([(1 - 3e-8) * (SP_potential[-1] * 1e-12), (1 + 3e-8) * (SP_potential[-1] * 1e-12)])  # Set limits for the upper part
ax1.set_xlim([6e-2, 1])  # Set x-axis limits in micrometers
ax2.set_ylim([-0.01, 0.04])  # Set limits for the lower part

# # Add lines for each detuning
# for i, (min_R, sp_potential_at_min_R, ps_potential_at_min_R, ss_potential_at_min_R) in enumerate(results):
#     detuning = Laser_Frequencies[i] - Delta
#     color = colors[i % len(colors)]  # Cycle through colors
#     linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
    
#     add_line_to_plot(ax1, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, white, linestyle, chosen_linewidth)
#     add_line_to_plot(ax2, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, chosen_linewidth)
    
# Add a vertical line at the flattening point
# ax1.axvline(flattening_point * 1e6, color='purple', linestyle='--', label=f'Flattening Point at R = {flattening_point*1e6:.2} µm')
# ax2.axvline(flattening_point * 1e6, color='purple', linestyle='--', label=f'Flattening Point at R = {flattening_point*1e6:.2} µm')

ax1.tick_params(axis='both', labelsize=16)
ax2.tick_params(axis='both', labelsize=16)
ax2.legend(loc='upper right', fontsize=16)

ax1.get_yaxis().get_major_formatter().set_useOffset(False)

ax1.set_xscale('log')
ax2.set_xscale('log')

# plt.suptitle('Potentials as a Function of Distance (R)', fontsize=20)
plt.show()
# Upper subplot
ax1.plot(R * 1e6, SP_potential * 1e-12, color='blue',linewidth = chosen_linewidth)
ax1.plot(R * 1e6, PS_potential * 1e-12, color='red',linewidth = chosen_linewidth)
# ax1.set_ylabel('Potential (THz)', fontsize=20)
# ax1.legend(loc='upper right')

# Lower subplot
ax2.plot(R * 1e6, SP_potential * 1e-12, label=f"SP", color='blue',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, PS_potential * 1e-12, label=f"PS", color='red',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, SS_potential * 1e-12, label=r"SS", color='orange',linewidth = chosen_linewidth)
# ax2.plot(R * 1e6, SP_potential * 1e-12, label=f"SP $= \Delta + U$", color='blue',linewidth = chosen_linewidth)
# ax2.plot(R * 1e6, PS_potential * 1e-12, label=f"PS $= \Delta - U$", color='red',linewidth = chosen_linewidth)
# ax2.plot(R * 1e6, SS_potential * 1e-12, label=r"SS $= \Delta - \sqrt{\Delta^2 + U^2}$", color='orange',linewidth = chosen_linewidth)
ax2.set_xlabel('Distance (µm)', fontsize=20)
ax2.set_ylabel('Potential (THz)', fontsize=20, loc="top")
ax2.yaxis.set_label_coords(-0.16, 1.23)
# Set limits
ax1.set_ylim([(1 - 3e-8) * (SP_potential[-1] * 1e-12), (1 + 3e-8) * (SP_potential[-1] * 1e-12)])  # Set limits for the upper part
ax1.set_xlim([6e-2, 1])  # Set x-axis limits in micrometers
ax2.set_ylim([-0.01, 0.04])  # Set limits for the lower part

# # Add lines for each detuning
# for i, (min_R, sp_potential_at_min_R, ps_potential_at_min_R, ss_potential_at_min_R) in enumerate(results):
#     detuning = Laser_Frequencies[i] - Delta
#     color = colors[i % len(colors)]  # Cycle through colors
#     linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
    
#     add_line_to_plot(ax1, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, white, linestyle, chosen_linewidth)
#     add_line_to_plot(ax2, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, chosen_linewidth)
    
# Add a vertical line at the flattening point
# ax1.axvline(flattening_point * 1e6, color='purple', linestyle='--', label=f'Flattening Point at R = {flattening_point*1e6:.2} µm')
# ax2.axvline(flattening_point * 1e6, color='purple', linestyle='--', label=f'Flattening Point at R = {flattening_point*1e6:.2} µm')

ax1.tick_params(axis='both', labelsize=16)
ax2.tick_params(axis='both', labelsize=16)
ax2.legend(loc='upper right', fontsize=16)

ax1.get_yaxis().get_major_formatter().set_useOffset(False)

ax1.set_xscale('log')
ax2.set_xscale('log')

# plt.suptitle('Potentials as a Function of Distance (R)', fontsize=20)
plt.show()






# ###### Plot the potentials without break ######
# fig, ax = plt.subplots(figsize=(10, 10))
# ax.plot(R * 1e6, SS_potential * 1e-12, label=r"SS $= \Delta - \sqrt{\Delta^2 + U^2}$", color='orange')
# ax.plot(R * 1e6, SP_potential * 1e-12, label=f"SP $= \Delta + U$", color='blue')
# ax.plot(R * 1e6, PS_potential * 1e-12, label=f"PS $= \Delta - U$", color='red')

# # Add a vertical line at the flattening point
# # ax.axvline(flattening_point * 1e6, color='purple', linestyle='--', label=f'Flattening Point at R = {flattening_point*1e6:.2} µm')


# # Add wavy transition lines
# for i, (min_R, sp_potential_at_min_R, ps_potential_at_min_R, ss_potential_at_min_R) in enumerate(results):
#     detuning = Laser_Frequencies[i] - Delta
#     color = colors[i % len(colors)]  # Cycle through colors
#     linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
    
#     add_line_to_plot(ax, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, chosen_linewidth)

# ax1.set_xlim([6e-2, 1])  # Set x-axis limits in micrometers
# # plt.ylim(-1e-7,1e-7)
# # plt.ylim(SPPotential(Delta, CouplingPotential(1e-7))* 1e-12-1e-7,SPPotential(Delta, CouplingPotential(1e-7))* 1e-12+1e-7)
# plt.xlabel('Distance (µm)', fontsize=16)
# plt.ylabel('Potential (THz)', fontsize=16)
# plt.title('Potentials as a Function of Distance (R)', fontsize=18)
# plt.legend(fontsize=14)
# plt.xscale('log')
# plt.show()




# ###### Collision Range Plot #########

# #Generate Laser Frequencies
# # Detuning_values = np.logspace(4, 8, 100)

# start = 0.01e6  # 0.05 MHz
# end = 30e6  # 30 MHz
# start_log = np.log10(start)
# end_log = np.log10(end)
# Detuning_values = np.logspace(start_log, end_log, 100)

# Laser_Frequency_Range = Delta + Detuning_values

# # Calculate Separations for each Laser Frequency
# Separations = np.zeros_like(Laser_Frequency_Range)

# for i, freq in enumerate(Laser_Frequency_Range):
#     Separations[i] = Collision_Separation_Calculator(R, freq, Delta)[0]

# # Convert Laser Frequencies to MHz and Separations to um
# x_data = (Laser_Frequency_Range - Delta) * 1e-6
# y_data = Separations * 1e6

# # Log-log transformation
# log_x_data = np.log(x_data)
# log_y_data = np.log(y_data)

# # Perform linear fit on the log-log data
# linear_fit_params = np.polyfit(log_x_data, log_y_data, 1)
# a, b = np.exp(linear_fit_params[1]), linear_fit_params[0]

# # Generate fitted data on the original scale
# fitted_separations_power_law = a * x_data**b


# # Plotting on regular axes with clearer colors and styles
# plt.figure(figsize=(11, 6))

# # Fitted Power-law Data
# plt.plot(x_data, fitted_separations_power_law, linewidth=8, label=f'Fit: $R = {a:.2f} \cdot \delta^{{{b:.2f}}}$')
# plt.plot(x_data, y_data, label='Original Data', color='red', linestyle='dotted',linewidth=8)  # blue with circle markers


# # Labels and Legends
# plt.xlabel('Detuning (MHz)', fontsize=20)
# plt.ylabel('Collision Separation (µm)', fontsize=20)
# plt.legend(fontsize=20)
# plt.xscale('log')
# plt.xticks(fontsize=18)
# plt.yticks(fontsize=18)

# # Grid for readability
# plt.grid(True)

# plt.show()


# # ###### Red Collision Range Plot #########

# # Generate Laser Frequencies
# Detuning_values = np.logspace(4, 8, 100)
# Laser_Frequency_Range = Delta - Detuning_values

# # Calculate Separations for each Laser Frequency
# Separations = np.zeros_like(Laser_Frequency_Range)

# for i, freq in enumerate(Laser_Frequency_Range):
#     Separations[i] = Red_Collision_Separation_Calculator(R, freq, Delta)[0]

# # Convert Laser Frequencies to MHz and Separations to um
# x_data = -(Laser_Frequency_Range - Delta) * 1e-6
# y_data = Separations * 1e6

# # Log-log transformation
# log_x_data = np.log(x_data)
# log_y_data = np.log(y_data)

# # Perform linear fit on the log-log data
# linear_fit_params = np.polyfit(log_x_data, log_y_data, 1)
# a, b = np.exp(linear_fit_params[1]), linear_fit_params[0]

# # Generate fitted data on the original scale
# fitted_separations_power_law = a * x_data**b

# # Plotting
# plt.figure(figsize=(10, 6))
# plt.plot(x_data, y_data, label='Original Data', color='b')
# plt.plot(x_data, fitted_separations_power_law, label=f'Fit: $y = {a:.3f} \cdot x^{b:.3f}$', color='r', linestyle='--')
# plt.xlabel('Detuning (MHz)', fontsize=14)
# plt.ylabel('Collision Separation (µm)', fontsize=14)
# plt.title('Collision Separation as a Function of Red Detuning', fontsize=16)
# plt.xscale('log')
# plt.yscale('log')
# plt.legend()
# plt.grid(True)

# # Setting more y-ticks
# y_ticks = np.logspace(np.log10(min(y_data)), np.log10(max(y_data)), num=10)
# plt.yticks(y_ticks)
# plt.gca().get_yaxis().set_major_formatter(FuncFormatter(format_func))

# plt.show()



# ###### Timescale Plot #######
# # Detuning_values = np.logspace(4, 8, 100)
# Laser_Frequency_Range = Delta + Detuning_values
# times = np.zeros_like(Laser_Frequency_Range)
# Laser_Energies = Laser_Frequency_Range * sc.h

# # Calculate flattening points for each detuning value
# flattening_points = np.zeros_like(Detuning_values)
# for i, detuning in enumerate(Detuning_values):
#     flattening_points[i] = find_flattening_point(R, SP_potential, 0.001 * detuning)

# # Calculate times using the flattening points
# for i, (Energy, sep) in enumerate(zip(Laser_Energies, Separations)):
#     try:
#         U_f = Delta * sc.h
#         U_i = Energy
#         X = flattening_points[i] - sep  # Corrected this line to use the i-th flattening point
#         times[i] = potential_hill_time(U_f, U_i, m, sep, flattening_points[i])  # Pass the correct i-th flattening point
#     except Exception as e:
#         times[i] = np.nan
#         print(f"Error calculating time for frequency {Laser_Frequency_Range[i]}: {e}")

# # Plotting times
# plt.figure(figsize=(10, 6))
# plt.plot(Detuning_values/ (U0 / sc.h), times * 1e6, linewidth=4, label=r"$T = \left( \frac{9 \cdot X^3\cdot m}{8\cdot \left( U(0) - U(t) \right)} \right)^{1/3}$", color='g')
# plt.xlabel('Detuning (Units of Trap Depth)', fontsize=20)
# plt.ylabel(r"Time ( $\mu$s )", fontsize=20)
# # plt.xscale('log')
# plt.legend(fontsize=20)
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)
# plt.gca().set_facecolor('none')

# plt.show()

