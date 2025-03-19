#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 16:18:40 2024

@author: joe v
"""

import numpy as np 
import matplotlib.pyplot as plt
from Trap_Simulation import *
import time
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap

Detuning = -150e6 # [Hz]

num_runs = 200

# List to store results from each run
results = []
position_histories = []
velocity_histories = []
total_trapped_list = []
timelength_list = []
collision_times = []
collision_coords_list = []
collision_distance_list = []
involved_in_collisions_list = []
number_of_collisions = []
Time_Stamps = []
SeparationHistory_list = []
collision_colour_list = []
ApproachVelocities_list = []
scattering_coords_list = []
scattering_times_list = []

iteration_times = []
start_time = time.time()


for run_number in range(num_runs):
    start_run_time = time.time()
    result, PositionHistory, TotalTrapped, timelength, collision_time, collision_coords, involved_in_collisions, collision_count, beta, TimeStamps, VelocityHistory, collision_distance, SeparationHistory, collision_colour, ApproachVelocities, scattering_coords, scattering_times = run_simulation(Detuning)
    run_time = time.time() - start_run_time
    results.append(result)
    position_histories.append(PositionHistory)
    velocity_histories.append(VelocityHistory)
    total_trapped_list.append(TotalTrapped)
    timelength_list.append(timelength)
    collision_times.append(collision_time)
    collision_coords_list.append(collision_coords)
    collision_distance_list.append(collision_distance)
    involved_in_collisions_list.append(involved_in_collisions)
    number_of_collisions.append(collision_count)
    Time_Stamps.append(TimeStamps)
    SeparationHistory_list.append(SeparationHistory)
    collision_colour_list.append(collision_colour)
    ApproachVelocities_list.append(ApproachVelocities)
    scattering_coords_list.append(scattering_coords)
    scattering_times_list.append(scattering_times)

    iteration_times.append(run_time)

    
    # Update time estimation every 5 iterations
    if (run_number + 1) % 1 == 0:
        print(f"Completed run {run_number + 1}/{num_runs}")

        # Calculate average time per iteration
        if len(iteration_times) > 1:
            average_time_per_iteration = np.mean(iteration_times)
        else:
            average_time_per_iteration = iteration_times[-1]  # Fallback to last iteration time

        # Estimate total remaining time
        completed_iterations = run_number + 1
        remaining_iterations = num_runs - completed_iterations
        estimated_time_remaining = average_time_per_iteration * remaining_iterations

        # Calculate the estimated finish time
        finish_time = datetime.now() + timedelta(seconds=estimated_time_remaining)
        finish_time_str = finish_time.strftime("%Y-%m-%d %H:%M:%S")

        # Print the estimated finish time and remaining time
        hours, remainder = divmod(estimated_time_remaining, 3600)
        minutes, seconds = divmod(remainder, 60)
        print(f"Time Remaining: {int(hours)}:{int(minutes)}:{int(seconds)}")
        print(f"Estimated Finish Time: {finish_time_str}")
        

average_atoms_left = np.mean(results)
average_number_of_collisions = np.mean(number_of_collisions)
one_atom_likelihood = (results.count(1) / num_runs) * 100

# Filter results based on even and odd TotalTrapped
even_trapped_results = [result for result, total_trapped in zip(results, total_trapped_list) if total_trapped % 2 == 0]
odd_trapped_results = [result for result, total_trapped in zip(results, total_trapped_list) if total_trapped % 2 != 0]

# Calculate average particles left for even and odd TotalTrapped
average_particles_left_even = np.mean(even_trapped_results)
average_particles_left_odd = np.mean(odd_trapped_results)





totaltime = time.time() - start_time
hours, remainder = divmod(totaltime, 3600)
minutes, seconds = divmod(remainder, 60)
print(f"Time Taken: {int(hours)}:{int(minutes)}:{int(seconds)}")

param_names = [
    "Single atom occupation probability",
    "",
    "Mean number of collisions",
    "",
    "Detuning (Trap center)",
    "Initial Trapped Particle Temperature",
    "Trap Depth",
    "Trap Width",
    "Initial Number of Particles",
    "Time span",
    "Initial Position Standard Deviation",
    "Damping Coefficient"
]

param_values = [
    f"{one_atom_likelihood:.1f}%",
    "",
    f"{average_number_of_collisions:.1f}",
    "",
    f"{Detuning * 1e-6:.2f} MHz",
    f"{Temp * 1e6:.2f} $\mu$K",
    f"{np.round(Trap_Depth_mK, 3)} mK / {np.round(Trap_Depth_MHz, 3)} MHz",
    f"{w0_tweezer * 1e6:.2f} $\mu$m",
    f"{np.average(total_trapped_list):.0f} $\pm$1",
    f"{ActualTime*1e3} ms",
    f"{PosStdDev * 1e6:.2f} $\mu$m",
    f"{beta:.2e}Kg/s ({(m/(beta+0.0001))*1e6:.0f} $\mu$s)"
]

max_name_length = max(len(name) for name in param_names)

parameters = [
    f"{name.ljust(max_name_length)} : {value}" if name else ""
    for name, value in zip(param_names, param_values)
]

param_text = "\n".join(parameters)


# Create the first figure for histogram and text plot
fig1, axs = plt.subplots(2, 1, figsize=(14, 10))

# Histogram plot for results
axs[0].hist(results, edgecolor='black', weights=np.ones_like(results) * 100.0 / num_runs)
axs[0].set_ylabel('Frequency (%)', fontsize=18)
axs[0].set_title(f"Trap Occupation Rate After {num_runs} Iterations", fontsize=18)

# Text plot for key parameters
axs[1].text(0.1, 0.5, param_text, fontsize=18, verticalalignment='center', fontfamily='monospace')
axs[1].axis('off')

plt.tight_layout()
plt.show()







# # Odd vs Even Histogram

# all_data = [results, even_trapped_results, odd_trapped_results]
# labels = ['All Results', 'Even TotalTrapped', 'Odd TotalTrapped']

# # Define bin edges
# bins = np.linspace(min(min(results), min(even_trapped_results), min(odd_trapped_results)),
#                     max(max(results), max(even_trapped_results), max(odd_trapped_results)), 21)

# fig, ax = plt.subplots(figsize=(14, 8))
# spacing = (bins[1] - bins[0]) * 0.3  # Adjust this value to increase or decrease spacing

# for i, (data, label) in enumerate(zip(all_data, labels)):
#     counts, _ = np.histogram(data, bins=bins)
#     # Convert counts to percentages
#     percentages = (counts / np.sum(counts)) * 100
#     offset = i * (bins[1] - bins[0] + spacing)  # Calculate offset for spacing
#     ax.bar(bins[:-1] + offset, percentages, width=bins[1] - bins[0], alpha=0.5, label=label, edgecolor='black')

# # Set labels and title
# ax.set_xlabel('Number of Particles Left', fontsize=18)
# ax.set_ylabel('Frequency (%)', fontsize=18)
# ax.set_title(f'Trap Occupation Rate After {num_runs} Iterations', fontsize=18)

# ax.legend(loc='upper right', fontsize=14)
# plt.tight_layout()
# plt.show()




######## Separation vs Detuning Plot ########


# # Detuning values (in Hz)
# detuning_values = [2e6, 5e6, 10e6, 20e6]

# # Number of samples
# num_samples = 10000

# # Generate x components and magnitudes for each detuning value
# x_components_energy = []
# magnitudes_energy = []

# for detuning in detuning_values:
#     x_vals = []
#     magnitudes = []
#     for _ in range(num_samples):
#         velocity = DetuningKEVelocity(detuning, m)
#         x_vals.append(velocity[0])
#         magnitudes.append(np.linalg.norm(velocity))
    
#     # Convert x components and magnitudes to energy in MHz
#     x_energy_mhz = [0.5 * m * x_val**2 / h * 1e-6 for x_val in x_vals]
#     magnitudes_energy_mhz = [0.5 * m * mag**2 / h * 1e-6 for mag in magnitudes]
    
#     x_components_energy.append(x_energy_mhz)
#     magnitudes_energy.append(magnitudes_energy_mhz)
    
#     # Debug: Print min, max, and sample values
#     print(f"Detuning = {detuning / 1e6:.1f} MHz:")
#     print(f"  Sample x component energies (MHz): {x_energy_mhz[:10]}")
#     print(f"  X component energies range from {min(x_energy_mhz):.2e} to {max(x_energy_mhz):.2e}")
#     print(f"  Sample magnitudes energies (MHz): {magnitudes_energy_mhz[:10]}")
#     print(f"  Magnitudes energies range from {min(magnitudes_energy_mhz):.2e} to {max(magnitudes_energy_mhz):.2e}")

# # Plotting
# fig, axs = plt.subplots(4, 1, figsize=(16, 16))

# # Flatten the axs array for easy iteration
# axs = axs.flatten()

# for i, (x_energy_mhz, magnitudes_energy_mhz, detuning) in enumerate(zip(x_components_energy, magnitudes_energy, detuning_values)):
#     # axs[i].hist(x_energy_mhz, bins=30, density=True, edgecolor='black', alpha=0.1, label='X Component Energy')
#     axs[i].hist(magnitudes_energy_mhz, bins=30, density=True, edgecolor='blue', alpha=1, label='Magnitude Energy')
#     axs[i].set_xlabel('Energy (MHz)')
#     axs[i].set_ylabel('Frequency')
#     axs[i].set_title(f'Detuning = {detuning / 1e6:.1f} MHz')
#     # axs[i].set_xlim(left=0, right=max(max(x_components_energy[-1]), max(magnitudes_energy[-1])))  # Adjust x-axis limits
#     axs[i].legend()

# plt.tight_layout()
# plt.show()



# ###### COM Distribution ##########
# fig, ax = plt.subplots(figsize=(10, 6))
# bin_range = (-1, 1)
# counts, bin_edges = np.histogram(COMvelocities, bins=40, range=bin_range)
# percentages = (counts / counts.sum()) * 100
# ax.bar(bin_edges[:-1], percentages, width=np.diff(bin_edges), edgecolor='black')
# ax.set_xlabel('Center of Mass Velocity (m/s)', fontsize=14)
# ax.set_ylabel('Frequency (%)', fontsize=14)
# ax.set_title('COM Velocities', fontsize=16)
# plt.show()


# ####### Detuning Distribution ##########
# fig, ax = plt.subplots(figsize=(11, 7))

# # Generate the histogram data
# counts, bin_edges = np.histogram(1e-6 * np.array(broadened_detunings), bins=30)
# percentages = (counts / counts.sum()) * 100

# # Split the bins into negative and non-negative detunings
# negative_indices = bin_edges[:-1] < 0
# positive_indices = ~negative_indices

# # Calculate the total percentage for each group
# negative_percentage = percentages[negative_indices].sum()
# positive_percentage = percentages[positive_indices].sum()

# # Plot negative detuning bars in lightcoral
# ax.bar(bin_edges[:-1][negative_indices], percentages[negative_indices], 
#         width=np.diff(bin_edges)[negative_indices], color='lightcoral', edgecolor='black', 
#         label=f"Red Collisions: {negative_percentage:.0f}%")

# # Plot non-negative detuning bars in the default color (blue)
# ax.bar(bin_edges[:-1][positive_indices], percentages[positive_indices], 
#         width=np.diff(bin_edges)[positive_indices], edgecolor='black',
#         label=f"Blue Collisions: {positive_percentage:.0f}%")

# # Add a vertical line at the Detuning value
# ax.axvline(1e-6 * Detuning, linestyle='--', linewidth=4, color='orange', 
#             label=rf"$\delta$ = {1e-6 * Detuning:.2f} MHz")

# # Set labels and title
# ax.set_xlabel('Detuning (MHz)', fontsize=20)
# ax.set_ylabel('Frequency (%)', fontsize=20)
# ax.set_title(f"New Model", fontsize=16)

# # Add a legend with the correct labels and Detuning line
# ax.legend(fontsize=20, facecolor='none')

# # Show the plot
# plt.show()

 

##### Pretty 3D Potential ##########


# # Define the grid of positions
# x = np.linspace(-5e-6, 5e-6, 1000)
# z = np.linspace(-15e-6, 15e-6, 1000)

# X, Z = np.meshgrid(x, z)
# Y = np.zeros_like(X)  # Fix y at 0 for simplicity

# # Initialize an array to store the potential
# Potential = np.zeros_like(X)

# # Calculate the potential at each grid point in MHz
# for i in range(len(x)):
#     for j in range(len(z)):
#         Potential[i, j] = potential(X[i, j], Y[i, j], Z[i, j], w0, zR, U0) / sc.h * 1e-6  # Convert potential to MHz

# # Create a custom colormap
# min_val, max_val = 0.0, 0.7
# n = 256
# orig_cmap = plt.cm.Blues_r
# colors = orig_cmap(np.linspace(min_val, max_val, n))
# cmap = LinearSegmentedColormap.from_list("mycmap", colors)

# # Plot the potential in 3D
# fig = plt.figure(figsize=(20,15), dpi=500)
# ax = fig.add_subplot(111, projection='3d')

# # Create the surface plot
# surf = ax.plot_surface(X, Z, Potential, cmap=cmap, linewidth=0.1, edgecolor='white')  # Add edgecolor for grid mesh

# # Customize the plot

# ax.set_box_aspect([4, 12, 6])  # [X, Y, Z] -> Make X twice as wide as Y and Z

# # Remove the background
# # ax.xaxis.pane.fill = False
# # ax.yaxis.pane.fill = False
# # ax.zaxis.pane.fill = False
# ax.grid(False)

# # Set the viewing angle
# ax.view_init(elev=20)  # Adjust the elevation and azimuth angle

# # Set axis labels with padding
# ax.set_xlabel(f"X Position (μm)", fontsize=20, labelpad=10)
# ax.set_ylabel(f"Z Position (μm)", fontsize=20, labelpad=15)
# ax.set_zlabel('Potential (MHz)', fontsize=20, labelpad=25)

# # Set X and Z ticks to whole numbers by scaling
# ax.set_xticks(np.arange(X.min(), X.max() + 1e-6, 5e-6))
# ax.set_xticklabels([f'{int(x * 1e6)}' for x in np.arange(X.min(), X.max() + 1e-6, 5e-6)], fontsize=15)

# ax.set_yticks(np.arange(z.min(), z.max() + 5e-6, 5e-6))
# ax.set_yticklabels([f'{int(z * 1e6)}' for z in np.arange(z.min(), z.max() + 5e-6, 5e-6)], fontsize=15)

# # Add a color bar
# cbar = fig.colorbar(surf, ax=ax, label='Potential (MHz)', fraction=0.026, pad=0.04, shrink=0.6)
# cbar.ax.tick_params(labelsize=14)  # Increase the font size of the colorbar ticks
# cbar.set_label('Potential (MHz)', fontsize=20) 
# # Show the plot
# plt.show()


##### Stark Shift Display (Blue/Red) #######

# from Trap_Functions import *
# from Trap_Simulation import *
# Detuning = 16e6 # [Hz]


# xvalues = np.linspace(-2.5e-6, 2.5e-6, 500)
# potential_at_center = potential(0, 0, 0, w0, zR, U0) / sc.h
# potential_at_pos1 = potential(-0.9e-6, 0, 0, w0, zR, U0) / sc.h
# potential_at_pos2 = potential(-0.75e-6, 0, 0, w0, zR, U0) / sc.h
# actual_detuning = []
# potential_values = []
# for x in xvalues:
#     pot_at_x = potential(x, 0, 0, w0, zR, U0) / sc.h
#     potential_values.append(pot_at_x*1e-6)
#     actual_detuning.append((Detuning + potential_at_center)*1e-6)


# pale_red = (1, 0, 0, 0.1)
# pale_blue = (0, 0, 1, 0.1)
# plt.figure(figsize=(10, 6),dpi=500)
# # plt.plot(xvalues * 1e6, actual_detuning, color = 'black', linestyle = 'dotted')
# plt.plot(xvalues * 1e6, potential_values, color = 'black')
# # plt.fill_between(xvalues * 1e6, actual_detuning, potential_values, where=(np.array(actual_detuning) >= np.array(potential_values)), facecolor=pale_blue, alpha=0.3, label = 'Blue-Detuned')
# # plt.fill_between(xvalues * 1e6, actual_detuning, potential_values, where=(np.array(actual_detuning) < np.array(potential_values)), facecolor=pale_red, alpha=0.3, label = 'Red-Detuned')
# plt.fill_between(xvalues * 1e6, potential_values, -13, facecolor='silver', alpha=0.3)
# # plt.vlines(0, ymin=potential_at_center*1e-6, ymax=(Detuning + potential_at_center)*1e-6, color='black', linestyle='--', linewidth=2)

# # Plot vertical lines and arrows
# plt.vlines(0, ymin=potential_at_center * 1e-6, ymax=(Detuning + potential_at_center) * 1e-6, 
#             color='black', linestyle='-', linewidth=2, label = f"{Detuning * 1e-6} MHz")
# plt.vlines(-0.9, ymin=potential_at_pos1 * 1e-6, ymax=(Detuning + potential_at_center) * 1e-6, 
#             color='blue', linestyle='-', linewidth=2, label = f"{(Detuning - (-potential_at_center + potential_at_pos1))*1e-6:.0f} MHz")
# plt.vlines(-0.75, ymin=potential_at_pos2 * 1e-6, ymax=(Detuning + potential_at_center) * 1e-6, 
#             color='blue', linestyle=(0,(5,1)), linewidth=2, label = f"{(Detuning - (-potential_at_center + potential_at_pos2))*1e-6:.0f} MHz")

# # Add arrows at the top of the blue lines
# plt.annotate('', xy=(-0.9, (Detuning + potential_at_center) * 1e-6), 
#               xytext=(-0.9, (Detuning + potential_at_center) * 1e-6 - 0.1), 
#               arrowprops=dict(facecolor='blue', arrowstyle='->', mutation_scale=20, linewidth=2))
# plt.annotate('', xy=(-0.75, (Detuning + potential_at_center) * 1e-6), 
#               xytext=(-0.75, (Detuning + potential_at_center) * 1e-6 - 0.1), 
#               arrowprops=dict(facecolor='blue', arrowstyle='->', mutation_scale=20, linewidth=2,color ='blue'))

# # Add small dots at the bottom of the blue lines
# plt.scatter([-0.9, -0.75], [potential_at_pos1 * 1e-6, potential_at_pos2 * 1e-6], color='red', s=80, zorder=5)


# # Add black arrows for the potential center
# plt.annotate('', xy=(0, (Detuning + potential_at_center) * 1e-6), 
#               xytext=(0, (Detuning + potential_at_center) * 1e-6 - 0.1), 
#               arrowprops=dict(facecolor='black', arrowstyle='->', mutation_scale=30, linewidth=2))
# plt.annotate('', xy=(0, (potential_at_center) * 1e-6), 
#               xytext=(0, (potential_at_center) * 1e-6 + 0.1), 
#               arrowprops=dict(facecolor='black', arrowstyle='->', mutation_scale=30, linewidth=2))


# # plt.text(0.25, (0.1 *Detuning + potential_at_center) * 1e-6, f"{Detuning * 1e-6} MHz", verticalalignment='center', fontsize=18, color='black')
# plt.text(0.03, (0.6 *Detuning + potential_at_center) * 1e-6, f"{Detuning * 1e-6} MHz", verticalalignment='center', fontsize=15, color='black')
# plt.text(-1.4, (0.8 *Detuning + potential_at_center) * 1e-6, f"{(Detuning - (-potential_at_center + potential_at_pos1))*1e-6:.0f} MHz", verticalalignment='center', fontsize=15, color='black')
# plt.text(-0.72, (0.8 *Detuning + potential_at_center) * 1e-6, f"{(Detuning - (-potential_at_center + potential_at_pos2))*1e-6:.0f} MHz", verticalalignment='center', fontsize=15, color='black')
# plt.xlabel('x Position (µm)', fontsize=20)
# plt.ylabel('Energy (MHz)', fontsize=20)
# # plt.title(f"Detuning as a Function of Position for {Detuning * 1e-6} MHz Detuned Light", fontsize=16)
# # plt.legend(fontsize=18, loc='lower right')
# plt.gca().set_facecolor('none')
# plt.show()


# ##### Initial Conditions Plot #######

# # Ensure all records have the same shape, for example (5, 3)
# consistent_length = 2  # Set this to the expected number of particles

# # Pad or truncate records to have consistent shape
# def adjust_shape(records, length):
#     adjusted_records = []
#     for record in records:
#         arr = np.array(record)
#         if arr.shape[0] < length:
#             # Pad with zeros if shorter
#             padded = np.pad(arr, ((0, length - arr.shape[0]), (0, 0)), mode='constant')
#         else:
#             # Truncate if longer
#             padded = arr[:length]
#         adjusted_records.append(padded)
#     return np.array(adjusted_records)

# positions_array = adjust_shape(Initial_Positions_Record, consistent_length)
# velocities_array = adjust_shape(Initial_Velocities_Record, consistent_length)
# trapped_array = Total_Trapped_Record


# # Function to plot histograms as frequency percentages
# def plot_histogram(data, bins, xlabel, ylabel, custom_ticks=None):
#     # Create the figure and axis with specified size
#     plt.figure(figsize=(11, 7))
    
#     # Calculate histogram data
#     counts, bin_edges = np.histogram(data, bins=bins)
#     percentages = (counts / counts.sum()) * 100  # Convert to percentages
    
#     # Calculate bin centers for accurate bar placement
#     bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    
#     # Plot the histogram
#     plt.bar(bin_centers, percentages, width=np.diff(bin_edges),edgecolor='black', align='center')
    
#     # Set labels and title with updated fontsize
#     plt.xlabel(xlabel, fontsize=20)
#     plt.ylabel(ylabel, fontsize=20)
#     # plt.title(title, fontsize=20)
    
#     plt.yticks(fontsize=16)
#     # Set custom ticks if provided
#     if custom_ticks:
#         plt.xticks(custom_ticks, fontsize=16)  # Set custom ticks
#     else:
#         plt.xticks(fontsize=16)  # Use default ticks
    
#     plt.show()

# # Plot histogram of x positions
# bins = np.array(np.linspace(-0.075, 0.075, 14))
# plot_histogram(positions_array[:, 0].flatten()*1e6, bins=30,
#                 xlabel=rf'X Position ($\mu$m)', ylabel='Frequency (%)')
#                 # title='Histogram of Initial X Positions')

# # Plot histogram of velocities
# bins = np.array(np.linspace(-5, 5, 10))
# plot_histogram(velocities_array.flatten(), bins=30,
#                 xlabel=rf'Velocity ($\mu$m/s)', ylabel='Frequency (%)')
#                 # title='Histogram of Initial Velocities')

# # Histogram of the total trapped particles
# bins = np.array([3.5, 4.5, 5.5, 6.5])  # Bins for 4, 5, 6 with edges around integers
# # Define custom ticks for the x-axis
# custom_ticks = [4, 5, 6]
# # plot_histogram(trapped_array, bins=bins,
# #                 xlabel='Total Trapped Particles', ylabel='Frequency (%)',
# #                 # title='Histogram of Total Trapped Particles',
# #                 custom_ticks=custom_ticks) 


######### New Model Picture #########


# from Trap_Functions import *
# from Trap_Simulation import *
# from matplotlib.ticker import ScalarFormatter
# from matplotlib.ticker import MaxNLocator
# D1_Freq = 377107.463380e9
# Detuning = 5e6 # [Hz]


# xvalues = np.linspace(-2.5e-6, 2.5e-6, 500)
# potential_at_center = potential(0, 0, 0, w0, zR, U0) / sc.h
# potential_at_pos1 = potential(particles_xpos*1e-6, 0, 0, w0, zR, U0) / sc.h
# actual_detuning = []
# potential_values = []
# for x in xvalues:
#     pot_at_x = potential(x, 0, 0, w0, zR, U0) / sc.h
#     potential_values.append(pot_at_x*1e-6)
#     actual_detuning.append((Detuning + potential_at_center)*1e-6)

# particles_xpos = -0.75
# Actual_Laser = D1_Freq - 2* potential_at_center
# upper_curve = np.zeros_like(potential_values) + Actual_Laser*1e-6 - potential_values
# # ####### broken y-axis ###########

# fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(9, 13), dpi=200) 

# # Upper subplot
# ax1.plot(xvalues * 1e6, upper_curve, color = 'black')
# # ax1.scatter([-0.9], [(Detuning - 2 * potential_at_pos1 + D1_Freq) * 1e-6], color='red', s=80, zorder=5)

# # # Blue Line
# # ax1.vlines(particles_xpos, ymin=potential_at_pos1 * 1e-6, ymax=(Actual_Laser - potential_at_pos1 +  Detuning + 2 * (potential_at_pos1 - potential_at_center)) * 1e-6, 
# #             color='blue', linestyle='-', linewidth=2, label = f"{(Detuning - (-potential_at_center + potential_at_pos1))*1e-6:.0f} MHz")
# # Black Line
# ax1.vlines(0, ymin=potential_at_center * 1e-6, ymax=(Detuning - potential_at_center + Actual_Laser) * 1e-6, 
#             color='black', linestyle='-', linewidth=2, label = f"{Detuning * 1e-6} MHz")
# # Black Arrow
# ax1.annotate('', xy=(0, (Detuning - potential_at_center + Actual_Laser) * 1e-6), 
#               xytext=(0, (Detuning - potential_at_center + Actual_Laser) * 1e-6 - 0.1), 
#               arrowprops=dict(facecolor='black', arrowstyle='->', mutation_scale=30, linewidth=2))
# # # Blue Arrow
# # ax1.annotate('', xy=(particles_xpos, (Actual_Laser - potential_at_pos1 +  Detuning + 2 * (potential_at_pos1 - potential_at_center)) * 1e-6), 
# #               xytext=(particles_xpos, (Actual_Laser - potential_at_pos1 +  Detuning + 2 * (potential_at_pos1 - potential_at_center)) * 1e-6 - 0.1), 
# #               arrowprops=dict(facecolor='blue', arrowstyle='->', mutation_scale=30, linewidth=2),color='blue')

# # Detuning Lines
# # ax1.vlines(particles_xpos-0.1, ymin=(-potential_at_pos1+Actual_Laser) * 1e-6, ymax=(Actual_Laser - potential_at_pos1 +  Detuning + 2 * (potential_at_pos1 - potential_at_center)) * 1e-6, 
# #             color='orange', linestyle='-', linewidth=2, label = f"{(Detuning - (-potential_at_center + potential_at_pos1))*1e-6:.0f} MHz")
# # ax1.text(-1.71, (Actual_Laser) * 1e-6 +17, f"{(Detuning + 2 * (potential_at_pos1 - potential_at_center) ) * 1e-6:.1f} MHz", verticalalignment='center', fontsize=15, color='black')

# ax1.vlines(0.1, ymin=(-potential_at_center +Actual_Laser) * 1e-6, ymax=(Detuning - potential_at_center + Actual_Laser) * 1e-6, 
#             color='orange', linestyle='-', linewidth=2, label = f"{(Detuning - (-potential_at_center + potential_at_pos1))*1e-6:.0f} MHz")
# ax1.text(0.2, (Actual_Laser) * 1e-6 +17, f"{(Detuning) * 1e-6:.1f} MHz", verticalalignment='center', fontsize=15, color='black')

# # ax2.vlines(particles_xpos-0.1, ymin=potential_at_center * 1e-6, ymax=potential_at_pos1 * 1e-6, 
# #             color='orange', linestyle='-', linewidth=2, label = f"{(Detuning - (-potential_at_center + potential_at_pos1))*1e-6:.0f} MHz")
# # ax2.text(-1.71, -10, f"{(potential_at_pos1 - potential_at_center ) * 1e-6:.1f} MHz", verticalalignment='center', fontsize=15, color='black')

# # Lower subplot
# ax2.plot(xvalues * 1e6, potential_values, color = 'black')
# ax2.set_xlabel('x Position (µm)', fontsize=20)
# ax2.set_ylabel('Potential (THz)', fontsize=20, loc="top")
# ax2.scatter([particles_xpos], [potential_at_pos1 * 1e-6], color='red', s=80, zorder=5)
# # ax2.vlines(particles_xpos, ymin=potential_at_pos1 * 1e-6, ymax=(Detuning + potential_at_pos1 + D1_Freq) * 1e-6, 
#             # color='blue', linestyle='-', linewidth=2, label = f"{(Detuning - (-potential_at_center + potential_at_pos1))*1e-6:.0f} MHz")
# ax2.vlines(0, ymin=potential_at_center * 1e-6, ymax=(Detuning + potential_at_center + D1_Freq) * 1e-6, 
#             color='black', linestyle='-', linewidth=2, label = f"{Detuning * 1e-6} MHz")

# # Set limits
# ax1.set_ylim([Actual_Laser*1e-6-1,24+Actual_Laser*1e-6])  # Set limits for the upper part
# ax2.set_ylim([-19, 6])  # Set limits for the lower part

# # ax1.tick_params(axis='both', labelsize=16)
# # ax2.tick_params(axis='both', labelsize=16)
# # ax2.legend(loc='upper right', fontsize=16)

# # Formatting the upper y-axis to show the offset
# # ax1.get_yaxis().get_major_formatter().set_useOffset(False)
# formatter = ScalarFormatter(useOffset=True)
# formatter.set_useOffset(Actual_Laser * 1e-6)  
# ax1.yaxis.set_major_formatter(formatter)

# ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
# ax2.yaxis.set_major_locator(MaxNLocator(integer=True))

# plt.show()


####### Time Step Distribution ##########
fig, ax = plt.subplots(figsize=(11, 7))

# Define the number of bins
num_bins = 30

min_rate = np.min(Timestep_Dist)  # Minimum positive value
max_rate = np.max(Timestep_Dist)
bin_edges = np.logspace(np.log10(min_rate), np.log10(max_rate), num_bins)

# Generate the histogram data with logarithmic bins
counts, bin_edges = np.histogram(Timestep_Dist, bins=bin_edges)
percentages = (counts / counts.sum()) * 100

# Create the bar plot with correct bar width
ax.bar(bin_edges[:-1], percentages, 
       width=np.diff(bin_edges), edgecolor='black')

# Set labels and title
ax.set_xlabel('Time Step size', fontsize=20)
ax.set_ylabel('Frequency (%)', fontsize=20)
ax.set_xscale('log')

# Show the plot
plt.title('Time Step Distribution', fontsize=24)
plt.show()



####### Scattering Rate Distribution ##########
num_bins = 30

# Generate the histogram data (using the same bin edges)
counts, bin_edges = np.histogram(scattering_rates, bins=num_bins)
counts_red, _ = np.histogram(scattering_rates_red, bins=bin_edges)

# Convert counts to percentages
percentages = (counts / counts.sum()) * 100
percentages_red = (counts_red / counts_red.sum()) * 100

# Create the figure and axis
fig, ax = plt.subplots(figsize=(11, 7))
ax.bar(bin_edges[:-1], percentages, width=np.diff(bin_edges), alpha=0.7, color='b', label='Blue Scattering Rates')
ax.bar(bin_edges[:-1], percentages_red, width=np.diff(bin_edges), alpha=0.7, color='r', label='Red Scattering Rates')
# Set labels and title
ax.set_xlabel('Scattering Rate', fontsize=20)
ax.set_ylabel('Frequency (%)', fontsize=20)
ax.set_title('Scattering Rate Distribution', fontsize=24)

# Add legend
ax.legend()

# Show the plot
plt.show()



### Collision Colour ######
from collections import Counter
flattened_collision_colours = [colour[0] for colour in collision_colour_list if colour]

# Count the occurrences of each color
collision_counts = Counter(flattened_collision_colours)

# Prepare data for plotting
colors = list(collision_counts.keys())
counts = list(collision_counts.values())

# Create the bar chart
plt.figure(figsize=(8, 6))
plt.bar(colors, counts, color=colors)

# Add labels and title
plt.xlabel('Collision Color', fontsize=14)
plt.ylabel('Count', fontsize=14)
plt.title('Collision Color Frequency', fontsize=16)

# Show the plot
plt.show()


### Scattering Times ######

scattering_times_flat = []

# Loop through each item in collision_times
for item in scattering_times_list:
    if isinstance(item, (list, np.ndarray)):  # Check if the item is a list or array
        scattering_times_flat.extend(item)     # Add all elements from the list/array to the flat list
    else:
        scattering_times_flat.append(item)     # If it's a single value, append it directly

# Convert to a numpy array after flattening
scattering_times_flat = np.array(scattering_times_flat)

# Create the figure for collision times histogram
fig2, ax2 = plt.subplots(figsize=(14, 8))

# Collision times histogram plot with weights
ax2.hist(1e3*scattering_times_flat, bins=40, edgecolor='black', 
          weights=np.ones_like(scattering_times_flat) * 100.0 / len(scattering_times_flat))
ax2.set_xlabel('Scattering Time (ms)', fontsize=18)
ax2.set_ylabel('Frequency (%)', fontsize=18)  # Adjusted to 'Frequency (%)' for clarity
ax2.set_title('Distribution of Scattering Times (Individual Particles)', fontsize=18)

plt.tight_layout()
plt.show()


### Collision Times ######

collision_times_flat = []

# Loop through each item in collision_times
for item in collision_times:
    if isinstance(item, (list, np.ndarray)):  # Check if the item is a list or array
        collision_times_flat.extend(item)     # Add all elements from the list/array to the flat list
    else:
        collision_times_flat.append(item)     # If it's a single value, append it directly

# Convert to a numpy array after flattening
collision_times_flat = np.array(collision_times_flat)

# Create the figure for collision times histogram
fig2, ax2 = plt.subplots(figsize=(14, 8))

# Collision times histogram plot with weights
ax2.hist(1e3*collision_times_flat, bins=40, edgecolor='black', 
          weights=np.ones_like(collision_times_flat) * 100.0 / len(collision_times_flat))
ax2.set_xlabel('Collision Time (ms)', fontsize=18)
ax2.set_ylabel('Frequency (%)', fontsize=18)  # Adjusted to 'Frequency (%)' for clarity
ax2.set_title('Distribution of Collision Times', fontsize=18)

plt.tight_layout()
plt.show()


### Collision Coords ######

collision_coords_flat = []

# Loop through each collision event in collision_coords_list
for event in collision_coords_list:
    for pair in event:
        for array in pair:
            # Flatten the individual arrays in the event
            collision_coords_flat.extend(array)

# Convert the flattened list into a numpy array
collision_coords_flat = np.array(collision_coords_flat)

# Ensure it's a 1D array for plotting
collision_coords_flat = collision_coords_flat.flatten()

# Now plot the histogram
fig2, ax2 = plt.subplots(figsize=(12, 8))

# Collision coordinates histogram plot with weights
ax2.hist(1e6 * collision_coords_flat, bins=20, edgecolor='black', 
         weights=np.ones_like(collision_coords_flat) * 100.0 / len(collision_coords_flat))
ax2.set_xlabel(f"Collision Coords ($\mu$m)", fontsize=18)
ax2.set_ylabel('Frequency (%)', fontsize=18)

plt.tight_layout()
plt.show()



### Collision Separation Distances ######
flat_collision_distance = []

# Loop through each list in collision_distance_list and flatten it
for sublist in collision_distance_list:
    for distance in sublist:
        flat_collision_distance.append(distance)

# Convert the flattened list into a numpy array
flat_collision_distance = np.array(flat_collision_distance)

# Now multiply by 1e6 (to convert to micrometers) and plot the histogram
fig2, ax2 = plt.subplots(figsize=(14, 8))

# Collision separation distance histogram plot with weights
ax2.hist(1e6 * flat_collision_distance, bins=20, edgecolor='black', 
         weights=np.ones_like(flat_collision_distance) * 100.0 / len(flat_collision_distance))
ax2.set_xlabel(f"Collision Separation Distance ($\mu$m)", fontsize=18)
ax2.set_ylabel('Frequency (%)', fontsize=18)  # Adjusted to 'Frequency (%)' for clarity

plt.tight_layout()
plt.show()


# ####### Doppler Distribution ##########

# # Generate the histogram data
# counts, bin_edges = np.histogram(1e-6*np.array(Detuning_record), bins=num_bins)

# # Convert counts to percentages
# percentages = (counts / counts.sum()) * 100

# # Create the figure and axis
# fig, ax = plt.subplots(figsize=(11, 7))
# ax.bar(bin_edges[:-1], percentages, width=np.diff(bin_edges), alpha=0.7)

# # Set labels and title
# ax.set_xlabel('Normalized Detuning', fontsize=20)
# ax.set_ylabel('Frequency (%)', fontsize=20)
# ax.set_title('Normalized Doppler Shifted Detuning Distribution', fontsize=24)

# # Show the plot
# plt.show()


# ####### Approach Velocities ##########

# flat_ApproachVelocities = []
# # Loop through each list in collision_distance_list and flatten it
# for sublist in ApproachVelocities_list:
#     for vel in sublist:
#         flat_ApproachVelocities.append(vel)

# # Convert the flattened list into a numpy array
# flat_ApproachVelocities = np.array(flat_ApproachVelocities)

# # Generate the histogram data
# counts, bin_edges = np.histogram(flat_ApproachVelocities, bins=num_bins)

# # Convert counts to percentages
# percentages = (counts / counts.sum()) * 100

# # Create the figure and axis
# fig, ax = plt.subplots(figsize=(11, 7))
# ax.bar(bin_edges[:-1], percentages, width=np.diff(bin_edges), alpha=0.7)

# # Set labels and title
# ax.set_xlabel('Approach Velocities (m/s)', fontsize=20)
# ax.set_ylabel('Frequency (%)', fontsize=20)
# ax.set_title('', fontsize=24)

# # Show the plot
# plt.show()



# # Optionally, plot the positions from a few of the runs
# for i in range(10):
    # plot_positions(i,position_histories[i], total_trapped_list[i], Time_Stamps[i], collision_times[i], collision_coords_list[i], w0_tweezer, collision_colour_list[i], SeparationHistory_list[i], scattering_times_list[i], scattering_coords_list[i])
    # plot_velocities(i,velocity_histories[i], total_trapped_list[i], Time_Stamps[i], collision_times[i], collision_coords_list[i], w0_tweezer, collision_colour_list[i])
    # plot_Separation(i, Time_Stamps[i], SeparationHistory_list[i])

import numpy as np
z = 8.0e-3
lamda = 784e-9
mfd = 5.0e-6
zr = np.pi/lamda * (mfd/2)**2
wz = mfd/2 * np.sqrt(1+ (z/zr)**2)
d = 2*wz
print('diameter =',d*1e3,'mm')

enlarged = d*7.5/1.8
print('enlarged beam =',enlarged*1e3,'mm')





