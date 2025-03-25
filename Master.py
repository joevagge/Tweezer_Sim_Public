#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 16:18:40 2024

@author: joe v
"""

"""
This script is used to run the simulation for a given laser detuning.

Key Parameters:
- Detuning: Laser detuning in Hz.
- num_runs: Number of simulation iterations.

Outputs:
- Plots: Histograms of trap occupation, collision times, etc.
- Data: Results stored in lists (e.g., `results`, `collision_times`).

There are many possible plots to output, depending on what data is required. Simply uncomment the desired sections. 
"""

import numpy as np 
import matplotlib.pyplot as plt
from Trap_Functions import *
from Trap_Simulation import *
import time
from datetime import datetime, timedelta
from matplotlib.colors import LinearSegmentedColormap


Detuning = 150e6 # [Hz].  # Set the cooling laser detuning 

num_runs = 20      # Set the number of iterations

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

    
    # This keeps track of how many iterations have been run and how long the simulation will take
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

# Filter results based on whether there was intially and even or odd number of particles in the trap
even_trapped_results = [result for result, total_trapped in zip(results, total_trapped_list) if total_trapped % 2 == 0]
odd_trapped_results = [result for result, total_trapped in zip(results, total_trapped_list) if total_trapped % 2 != 0]
average_particles_left_even = np.mean(even_trapped_results)
average_particles_left_odd = np.mean(odd_trapped_results)


##### Output Parameters  #####
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


################# ----- PLOTS ------ ###########################



###### Main Histogram #######

## This is the baseline histogram, it shows the distribution of number of atoms remaining, and prints some important parameters

fig1, axs = plt.subplots(2, 1, figsize=(14, 10))
axs[0].hist(results, edgecolor='black', weights=np.ones_like(results) * 100.0 / num_runs)
axs[0].set_ylabel('Frequency (%)', fontsize=18)
axs[0].set_title(f"Trap Occupation Rate After {num_runs} Iterations", fontsize=18)
axs[1].text(0.1, 0.5, param_text, fontsize=18, verticalalignment='center', fontfamily='monospace')
axs[1].axis('off')
plt.tight_layout()
plt.show()


######## Even vs Odd Histogram ######## 

## Histogram with Even and Odd Initial Paritcles no.

all_data = [results, even_trapped_results, odd_trapped_results]
labels = ['All Results', 'Even TotalTrapped', 'Odd TotalTrapped']

bins = np.linspace(min(min(results), min(even_trapped_results), min(odd_trapped_results)),
                    max(max(results), max(even_trapped_results), max(odd_trapped_results)), 21)

fig, ax = plt.subplots(figsize=(14, 8))
spacing = (bins[1] - bins[0]) * 0.3  # Adjust this value to increase or decrease spacing

for i, (data, label) in enumerate(zip(all_data, labels)):
    counts, _ = np.histogram(data, bins=bins)
    # Convert counts to percentages
    percentages = (counts / np.sum(counts)) * 100
    offset = i * (bins[1] - bins[0] + spacing)  # Calculate offset for spacing
    ax.bar(bins[:-1] + offset, percentages, width=bins[1] - bins[0], alpha=0.5, label=label, edgecolor='black')

ax.set_xlabel('Number of Particles Left', fontsize=18)
ax.set_ylabel('Frequency (%)', fontsize=18)
ax.set_title(f'Trap Occupation Rate After {num_runs} Iterations', fontsize=18)

ax.legend(loc='upper right', fontsize=14)
plt.tight_layout()
plt.show()


###### COM Distribution ##########

## Distribution of the velocity of the centre of mass between two particles 
fig, ax = plt.subplots(figsize=(10, 6))
bin_range = (-1, 1)
counts, bin_edges = np.histogram(COMvelocities, bins=40, range=bin_range)
percentages = (counts / counts.sum()) * 100
ax.bar(bin_edges[:-1], percentages, width=np.diff(bin_edges), edgecolor='black')
ax.set_xlabel('Center of Mass Velocity (m/s)', fontsize=14)
ax.set_ylabel('Frequency (%)', fontsize=14)
ax.set_title('COM Velocities', fontsize=16)
plt.show()


# ####### Detuning Distribution ##########

# The detuning experienced by the atom varies slightly, depending on where the atom is in the trap (how strong the AC stark shift is)
# This plot shows the distribution of apparent detunings.

fig, ax = plt.subplots(figsize=(11, 7))
counts, bin_edges = np.histogram(1e-6 * np.array(broadened_detunings), bins=30)
percentages = (counts / counts.sum()) * 100

negative_indices = bin_edges[:-1] < 0  # Split the bins into negative and non-negative detunings, i.e. red or blue
positive_indices = ~negative_indices

negative_percentage = percentages[negative_indices].sum()   # Calculate the total percentage for each group
positive_percentage = percentages[positive_indices].sum()

ax.bar(bin_edges[:-1][negative_indices], percentages[negative_indices],  # Plot negative detuning bars in red
        width=np.diff(bin_edges)[negative_indices], color='lightcoral', edgecolor='black', 
        label=f"Red Collisions: {negative_percentage:.0f}%")

ax.bar(bin_edges[:-1][positive_indices], percentages[positive_indices], # Plot non-negative detuning bars in the default color (blue)
        width=np.diff(bin_edges)[positive_indices], edgecolor='black',
        label=f"Blue Collisions: {positive_percentage:.0f}%")

ax.axvline(1e-6 * Detuning, linestyle='--', linewidth=4, color='orange',  # Add a vertical line at the default Detuning value
            label=rf"$\delta$ = {1e-6 * Detuning:.2f} MHz")

ax.set_xlabel('Detuning (MHz)', fontsize=20)
ax.set_ylabel('Frequency (%)', fontsize=20)
ax.set_title(f"New Model", fontsize=16)
ax.legend(fontsize=20, facecolor='none')
plt.show()

 

##### Pretty 3D Potential ##########

# This is just a graphical illustration of the Tweezer Potential in 3D

x = np.linspace(-5e-6, 5e-6, 1000)
z = np.linspace(-15e-6, 15e-6, 1000)
X, Z = np.meshgrid(x, z)
Y = np.zeros_like(X)  # Fix y at 0 for simplicity
Potential = np.zeros_like(X)

for i in range(len(x)):
    for j in range(len(z)):
        Potential[i, j] = potential(X[i, j], Y[i, j], Z[i, j], w0_tweezer, zR_tweezer, U0) / sc.h * 1e-6  # Convert potential to MHz

min_val, max_val = 0.0, 0.7  # Create a custom colormap
n = 256
orig_cmap = plt.cm.Blues_r
colors = orig_cmap(np.linspace(min_val, max_val, n))
cmap = LinearSegmentedColormap.from_list("mycmap", colors)

fig = plt.figure(figsize=(20,15), dpi=500)
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, Z, Potential, cmap=cmap, linewidth=0.1, edgecolor='white')  # Add edgecolor for grid mesh

ax.set_box_aspect([4, 12, 6])  # [X, Y, Z] -> Make X twice as wide as Y and Z
ax.grid(False)

ax.view_init(elev=20)  # Set the viewing angle

ax.set_xlabel(f"X Position (μm)", fontsize=20, labelpad=10)
ax.set_ylabel(f"Z Position (μm)", fontsize=20, labelpad=15)
ax.set_zlabel('Potential (MHz)', fontsize=20, labelpad=25)

ax.set_xticks(np.arange(X.min(), X.max() + 1e-6, 5e-6))
ax.set_xticklabels([f'{int(x * 1e6)}' for x in np.arange(X.min(), X.max() + 1e-6, 5e-6)], fontsize=15)

ax.set_yticks(np.arange(z.min(), z.max() + 5e-6, 5e-6))
ax.set_yticklabels([f'{int(z * 1e6)}' for z in np.arange(z.min(), z.max() + 5e-6, 5e-6)], fontsize=15)

cbar = fig.colorbar(surf, ax=ax, label='Potential (MHz)', fraction=0.026, pad=0.04, shrink=0.6)
cbar.ax.tick_params(labelsize=14)  # Increase the font size of the colorbar ticks
cbar.set_label('Potential (MHz)', fontsize=20) 
plt.show()


# ##### Initial Conditions Distributions #######

# Ensure all records have the same shape
consistent_length = 2  # Set this to the expected number of particles
# Pad or truncate records to have consistent shape
def adjust_shape(records, length):
    adjusted_records = []
    for record in records:
        arr = np.array(record)
        if arr.shape[0] < length:
            padded = np.pad(arr, ((0, length - arr.shape[0]), (0, 0)), mode='constant') # Pad with zeros if shorter
        else:
            padded = arr[:length] # Truncate if longer
        adjusted_records.append(padded)
    return np.array(adjusted_records)

positions_array = adjust_shape(Initial_Positions_Record, consistent_length)
velocities_array = adjust_shape(Initial_Velocities_Record, consistent_length)
trapped_array = Total_Trapped_Record

# Function to plot histograms as frequency percentages
def plot_histogram(data, bins, xlabel, ylabel, custom_ticks=None):
    plt.figure(figsize=(11, 7))
    counts, bin_edges = np.histogram(data, bins=bins)
    percentages = (counts / counts.sum()) * 100  # Convert to percentages
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    plt.bar(bin_centers, percentages, width=np.diff(bin_edges),edgecolor='black', align='center')
    plt.xlabel(xlabel, fontsize=20)
    plt.ylabel(ylabel, fontsize=20)    
    plt.yticks(fontsize=16)
    if custom_ticks:
        plt.xticks(custom_ticks, fontsize=16)  # Set custom ticks
    else:
        plt.xticks(fontsize=16)  # Use default ticks
    plt.show()

# Plot histogram of x positions
bins = np.array(np.linspace(-0.075, 0.075, 14))
plot_histogram(positions_array[:, 0].flatten()*1e6, bins=30,
                xlabel=rf'X Position ($\mu$m)', ylabel='Frequency (%)')
                # title='Histogram of Initial X Positions')

# Plot histogram of velocities
bins = np.array(np.linspace(-5, 5, 10))
plot_histogram(velocities_array.flatten(), bins=30,
                xlabel=rf'Velocity ($\mu$m/s)', ylabel='Frequency (%)')
                # title='Histogram of Initial Velocities')

# Histogram of the total trapped particles
# bins = np.array([3.5, 4.5, 5.5, 6.5])  # Bins for 4, 5, 6 with edges around integers
# custom_ticks = [4, 5, 6]
# plot_histogram(trapped_array, bins=bins,
#                 xlabel='Total Trapped Particles', ylabel='Frequency (%)',
#                 # title='Histogram of Total Trapped Particles',
#                 custom_ticks=custom_ticks) 


####### Time Step Distribution ##########

fig, ax = plt.subplots(figsize=(11, 7))
num_bins = 30
min_rate = np.min(Timestep_Dist)  
max_rate = np.max(Timestep_Dist)
bin_edges = np.logspace(np.log10(min_rate), np.log10(max_rate), num_bins)
counts, bin_edges = np.histogram(Timestep_Dist, bins=bin_edges)
percentages = (counts / counts.sum()) * 100
ax.bar(bin_edges[:-1], percentages, 
       width=np.diff(bin_edges), edgecolor='black')
ax.set_xlabel('Time Step size', fontsize=20)
ax.set_ylabel('Frequency (%)', fontsize=20)
ax.set_xscale('log')
plt.title('Time Step Distribution', fontsize=24)
plt.show()



####### Scattering Rate Distribution ##########

num_bins = 30
counts, bin_edges = np.histogram(scattering_rates, bins=num_bins)
counts_red, _ = np.histogram(scattering_rates_red, bins=bin_edges)
percentages = (counts / counts.sum()) * 100
percentages_red = (counts_red / counts_red.sum()) * 100
fig, ax = plt.subplots(figsize=(11, 7))
ax.bar(bin_edges[:-1], percentages, width=np.diff(bin_edges), alpha=0.7, color='b', label='Blue Scattering Rates')
ax.bar(bin_edges[:-1], percentages_red, width=np.diff(bin_edges), alpha=0.7, color='r', label='Red Scattering Rates')
# Set labels and title
ax.set_xlabel('Scattering Rate', fontsize=20)
ax.set_ylabel('Frequency (%)', fontsize=20)
ax.set_title('Scattering Rate Distribution', fontsize=24)
ax.legend()
plt.show()



####### Collision Colour #######

## Comparison of number of red-detuned vs blue-detuned collisions

from collections import Counter
flattened_collision_colours = [colour[0] for colour in collision_colour_list if colour]
collision_counts = Counter(flattened_collision_colours)
colors = list(collision_counts.keys())
counts = list(collision_counts.values())
plt.figure(figsize=(8, 6))
plt.bar(colors, counts, color=colors)
plt.xlabel('Collision Color', fontsize=14)
plt.ylabel('Count', fontsize=14)
plt.title('Collision Color Frequency', fontsize=16)
plt.show()


####### Scattering Time Distribution #######

scattering_times_flat = []
for item in scattering_times_list:
    if isinstance(item, (list, np.ndarray)):  
        scattering_times_flat.extend(item)     
    else:
        scattering_times_flat.append(item)     
scattering_times_flat = np.array(scattering_times_flat)
fig2, ax2 = plt.subplots(figsize=(14, 8))
ax2.hist(1e3*scattering_times_flat, bins=40, edgecolor='black', 
          weights=np.ones_like(scattering_times_flat) * 100.0 / len(scattering_times_flat))
ax2.set_xlabel('Scattering Time (ms)', fontsize=18)
ax2.set_ylabel('Frequency (%)', fontsize=18)  # Adjusted to 'Frequency (%)' for clarity
ax2.set_title('Distribution of Scattering Times (Individual Particles)', fontsize=18)

plt.tight_layout()
plt.show()


####### Collision Times #######

collision_times_flat = []
for item in collision_times:
    if isinstance(item, (list, np.ndarray)): 
        collision_times_flat.extend(item)     
    else:
        collision_times_flat.append(item)     
collision_times_flat = np.array(collision_times_flat)
fig2, ax2 = plt.subplots(figsize=(14, 8))
ax2.hist(1e3*collision_times_flat, bins=40, edgecolor='black', 
          weights=np.ones_like(collision_times_flat) * 100.0 / len(collision_times_flat))
ax2.set_xlabel('Collision Time (ms)', fontsize=18)
ax2.set_ylabel('Frequency (%)', fontsize=18)  # Adjusted to 'Frequency (%)' for clarity
ax2.set_title('Distribution of Collision Times', fontsize=18)
plt.tight_layout()
plt.show()


####### Collision Coords #######

collision_coords_flat = []
for event in collision_coords_list:
    for pair in event:
        for array in pair:
            # Flatten the individual arrays in the event
            collision_coords_flat.extend(array)
collision_coords_flat = np.array(collision_coords_flat)
collision_coords_flat = collision_coords_flat.flatten()
fig2, ax2 = plt.subplots(figsize=(12, 8))
ax2.hist(1e6 * collision_coords_flat, bins=20, edgecolor='black', 
         weights=np.ones_like(collision_coords_flat) * 100.0 / len(collision_coords_flat))
ax2.set_xlabel(f"Collision Coords ($\mu$m)", fontsize=18)
ax2.set_ylabel('Frequency (%)', fontsize=18)

plt.tight_layout()
plt.show()



####### Collision Separation Distances #######

flat_collision_distance = []
for sublist in collision_distance_list:
    for distance in sublist:
        flat_collision_distance.append(distance)
flat_collision_distance = np.array(flat_collision_distance)
fig2, ax2 = plt.subplots(figsize=(14, 8))
ax2.hist(1e6 * flat_collision_distance, bins=20, edgecolor='black', 
         weights=np.ones_like(flat_collision_distance) * 100.0 / len(flat_collision_distance))
ax2.set_xlabel(f"Collision Separation Distance ($\mu$m)", fontsize=18)
ax2.set_ylabel('Frequency (%)', fontsize=18)  # Adjusted to 'Frequency (%)' for clarity
plt.tight_layout()
plt.show()


# ####### Doppler Distribution ##########

# Distribution of the doppler shifted detuning experienced by a particle

counts, bin_edges = np.histogram(1e-6*np.array(Detuning_record), bins=num_bins)
percentages = (counts / counts.sum()) * 100
fig, ax = plt.subplots(figsize=(11, 7))
ax.bar(bin_edges[:-1], percentages, width=np.diff(bin_edges), alpha=0.7)
ax.set_xlabel('Normalized Detuning', fontsize=20)
ax.set_ylabel('Frequency (%)', fontsize=20)
ax.set_title('Normalized Doppler Shifted Detuning Distribution', fontsize=24)
plt.show()


# ####### Approach Velocities ##########

flat_ApproachVelocities = []
for sublist in ApproachVelocities_list:
    for vel in sublist:
        flat_ApproachVelocities.append(vel)

flat_ApproachVelocities = np.array(flat_ApproachVelocities)
counts, bin_edges = np.histogram(flat_ApproachVelocities, bins=num_bins)
percentages = (counts / counts.sum()) * 100
fig, ax = plt.subplots(figsize=(11, 7))
ax.bar(bin_edges[:-1], percentages, width=np.diff(bin_edges), alpha=0.7)
ax.set_xlabel('Approach Velocities (m/s)', fontsize=20)
ax.set_ylabel('Frequency (%)', fontsize=20)
ax.set_title('', fontsize=24)
plt.show()



####### Plot the position as a function of time from a few of the iterations #######

for i in range(10):
    plot_positions(i,position_histories[i], total_trapped_list[i], Time_Stamps[i], collision_times[i], collision_coords_list[i], w0_tweezer, collision_colour_list[i], SeparationHistory_list[i], scattering_times_list[i], scattering_coords_list[i])
    plot_velocities(i,velocity_histories[i], total_trapped_list[i], Time_Stamps[i], collision_times[i], collision_coords_list[i], w0_tweezer, collision_colour_list[i])
    # plot_Separation(i, Time_Stamps[i], SeparationHistory_list[i])




