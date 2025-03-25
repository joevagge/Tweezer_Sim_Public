#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 15:24:03 2024

@author: joe
"""

"""
This script calculates and plots interatomic potential energy curves (PECs) for atoms in a tweezer trap.


Key Functions:
- CouplingPotential(R): Calculates the coupling potential at distance R.
- PEC_Dynamics(R_initial, v_initial, tau): Simulates motion along the PEC.

"""

import numpy as np 
import scipy.constants as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import random
from Trap_Functions import * 


# Constants
D1_Freq = 377107.463380e9 # Transition frequency 
D2_Freq = 384230.4844685e9 # Transition frequency 
tdme = 2.992 * sc.eV* sc.physical_constants['Bohr radius'][0] # [C·m] Transition Dipole Matrix Element
Delta = 377.107463380 * 1e12  # [Hz], converted from THz to Hz
eps0 = sc.epsilon_0
linewidth = 2 * sc.pi * 5.7500 * 1e6  # [Hz]
m = 59*sc.atomic_mass  

# Define the CouplingPotential function in Hz
def CouplingPotential(R):
    # Calculate the potential U in Joules
    U = ((1) * tdme**2) / (4 * sc.pi * eps0 * R**3)
    
    # Convert the potential U from Joules to frequency in Hz
    f_Hz = U / sc.h
    return f_Hz

# Define the potential functions
def SSPotential(Delta, U):
    return Delta - np.sqrt(Delta**2 + U**2)

def PPPotential(Delta, U):
    return Delta + np.sqrt(Delta**2 + U**2)

def SPPotential(Delta, U):
    return Delta + U

def PSPotential(Delta, U):
    return Delta - U

def difference(r, Laser_Freq, Delta):
    u = CouplingPotential(r)
    sp_potential = Delta + u
    ss_potential = Delta - np.sqrt(Delta**2 + u**2)
    return np.abs((sp_potential - ss_potential) - Laser_Freq)

def Red_difference(r, Laser_Freq, Delta):
    u = CouplingPotential(r)
    ps_potential = Delta - u
    ss_potential = Delta - np.sqrt(Delta**2 + u**2)
    return np.abs((ps_potential - ss_potential) - Laser_Freq)

# This function is just for calculating the separation distances at which a given frequency of lght becomes resonant with the transition, plotted with the green vertical lines
def Collision_Separation_Calculator(R, Laser_Freq, Delta):
    # Evaluate the difference function
    differences = np.array([difference(r, Laser_Freq, Delta) for r in R])

    # Find the R value where the difference is minimized
    min_diff_index = np.argmin(differences)
    min_R = R[min_diff_index]

    # Calculate the y-values at min_R
    sp_potential_at_min_R = SPPotential(Delta, CouplingPotential(min_R))
    ps_potential_at_min_R = PSPotential(Delta, CouplingPotential(min_R))
    ss_potential_at_min_R = SSPotential(Delta, CouplingPotential(min_R))
    
    return min_R, sp_potential_at_min_R, ps_potential_at_min_R, ss_potential_at_min_R

# This is for plotting the lines calculated in the function above
def add_line_to_plot(ax, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, linewidth):
    with mpl.rc_context({'path.sketch': (4, 70, 1)}):
        ax.plot([min_R * 1e6, min_R * 1e6], 
                [ss_potential_at_min_R * 1e-12, sp_potential_at_min_R * 1e-12], 
                color=color, linestyle=linestyle, 
                label=f'$\delta$ = {detuning*1e-6:.0f}$\,$MHz $\Rightarrow$ {min_R*1e6:.2}$\,$µm',
                linewidth=linewidth)                

# Formatting y-tick labels
def format_func(value, tick_number):
    return f'{value:.3f}'  # Customize this formatting as needed


####### The following code produces some sample data in order to plot the interatomic PECs #######

R = np.logspace(-9,  -6, 500)
Detuning_1 = 50e6
Detuning_2 = 50e6
Detuning_3 = 50e6
Laser_Freq_1 = D1_Freq + Detuning_1
Laser_Freq_2 = D1_Freq + Detuning_2
Laser_Freq_3 = D1_Freq + Detuning_3
Laser_Frequencies = [Laser_Freq_1, Laser_Freq_2, Laser_Freq_3]
results = [Collision_Separation_Calculator(R, freq, Delta) for freq in Laser_Frequencies]
SP_potential = SPPotential(Delta, CouplingPotential(R))
PS_potential = PSPotential(Delta, CouplingPotential(R))
SS_potential = SSPotential(Delta, CouplingPotential(R))
colors = ['darkgreen', 'darkgreen', 'darkgreen']
linestyles = ['solid', 'dashed', 'dotted']


####### Plotting Interatomic PECs with broken y-axis ###########
chosen_linewidth = 3
fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(10, 11), dpi=500) 

# Upper subplot
ax1.plot(R * 1e6, SP_potential * 1e-12, color='blue',linewidth = chosen_linewidth)
ax1.plot(R * 1e6, PS_potential * 1e-12, color='red',linewidth = chosen_linewidth)

# Lower subplot
ax2.plot(R * 1e6, SP_potential * 1e-12, label=f"SP $= \Delta + U$", color='blue',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, PS_potential * 1e-12, label=f"PS $= \Delta - U$", color='red',linewidth = chosen_linewidth)
ax2.plot(R * 1e6, SS_potential * 1e-12, label=r"SS $= \Delta - \sqrt{\Delta^2 + U^2}$", color='orange',linewidth = chosen_linewidth)
ax2.set_xlabel('Distance (µm)', fontsize=20)
ax2.set_ylabel('Potential (THz)', fontsize=20, loc="top")
ax2.yaxis.set_label_coords(-0.16, 1.23)
# Set limits
ax1.set_ylim([(1 - 3e-8) * (SP_potential[-1] * 1e-12), (1 + 3e-8) * (SP_potential[-1] * 1e-12)])  # Set limits for the upper subplot
ax1.set_xlim([6e-2, 1])  # Set x-axis limits in micrometers
ax2.set_ylim([-0.01, 0.04])  # Set x-axis limits for the lower subplot

# Add lines for each value of detuning
for i, (min_R, sp_potential_at_min_R, ps_potential_at_min_R, ss_potential_at_min_R) in enumerate(results):
    detuning = Laser_Frequencies[i] - Delta
    color = colors[i % len(colors)]  # Cycle through colors
    linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
    
    add_line_to_plot(ax1, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, chosen_linewidth)
    add_line_to_plot(ax2, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, chosen_linewidth)
    
ax1.tick_params(axis='both', labelsize=16)
ax2.tick_params(axis='both', labelsize=16)
ax2.legend(loc='upper right', fontsize=16)
ax1.get_yaxis().get_major_formatter().set_useOffset(False)
ax1.set_xscale('log')
ax2.set_xscale('log')
plt.show()


###### Same plot as above, but without page break ######
fig, ax = plt.subplots(figsize=(10, 10))
ax.plot(R * 1e6, SS_potential * 1e-12, label=r"SS $= \Delta - \sqrt{\Delta^2 + U^2}$", color='orange')
ax.plot(R * 1e6, SP_potential * 1e-12, label=f"SP $= \Delta + U$", color='blue')
ax.plot(R * 1e6, PS_potential * 1e-12, label=f"PS $= \Delta - U$", color='red')

# Add wavy transition lines
for i, (min_R, sp_potential_at_min_R, ps_potential_at_min_R, ss_potential_at_min_R) in enumerate(results):
    detuning = Laser_Frequencies[i] - Delta
    color = colors[i % len(colors)]  # Cycle through colors
    linestyle = linestyles[i % len(linestyles)]  # Cycle through linestyles
    
    add_line_to_plot(ax, min_R, ss_potential_at_min_R, sp_potential_at_min_R, detuning, color, linestyle, chosen_linewidth)

ax1.set_xlim([6e-2, 1])  # Set x-axis limits in micrometers
plt.xlabel('Distance (µm)', fontsize=16)
plt.ylabel('Potential (THz)', fontsize=16)
plt.title('Potentials as a Function of Distance (R)', fontsize=18)
plt.legend(fontsize=14)
plt.xscale('log')
plt.show()



#### Now Modelling the dynamics of a particle moving along the PEC ######

## This function describes the motion along the repulsive (blue) PEC
def PEC_Dynamics(R_initial, v_initial, tau):

    delta_t = 0.1 * tau
    v_approach = v_initial
    R = R_initial
    t = 0  # Initial time
    
    while True:
        p = random.uniform(0, 1)
        if p < delta_t / tau:
            break 

        F = -(3 * tdme**2) / (4 * sc.pi * eps0 * R**4)  
        a = F / m  
        v_approach += a * delta_t
        R -= v_approach * delta_t
        t += delta_t

        if t > 500e-9:
            print('collision took too long')
            break        

    return R, v_approach

## This function describes the motion along the attractive (red) PEC
def PEC_Dynamics_Red(R_initial, v_initial, tau):

    delta_t = 0.1 * tau
    v_approach = v_initial
    R = R_initial
    t = 0  # Initial time
    
    while True:
        p = random.uniform(0, 1)
        if p < delta_t / tau:
            break 

        F = (3 * tdme**2) / (4 * sc.pi * eps0 * R**4)  
        a = F / m 
        v_approach += a * delta_t
        R -= v_approach * delta_t

        # Increment time
        t += delta_t
        
        if t > 500e-9:
            print('collision took too long')
            break

    return R, v_approach

from datetime import datetime, timedelta



##### The rest of this script is used to test the above functions ####

def PEC_Dynamics_Test(R_initial, v_initial, tau):
    tvalues = []
    Rvalues = []
    Uvalues = []
    Vvalues = []
    
    delta_t = 0.001 * tau
    v_approach = v_initial
    R = R_initial
    t = 0  # Initial time
    
    start_time = datetime.now()
    
    while True:
        p = random.uniform(0, 1)
        if p < delta_t / tau:
            end_time = (datetime.now() - start_time).total_seconds()

            break 

        U = SPPotential(Delta, CouplingPotential(R)) - SPPotential(Delta, CouplingPotential(1e-3))
        F = (3 * tdme**2) / (4 * sc.pi * eps0 * R**4)  
        a = F / m 
        
        # Update velocity and position
        v_approach += a * delta_t
        R -= v_approach * delta_t
        
        if R < 0:
            end_time = (datetime.now() - start_time).total_seconds()
            break
        
        # Append results for each time step
        tvalues.append(t)
        Rvalues.append(R)
        Uvalues.append(U * 1e-6)  # (in MHz)
        Vvalues.append(v_approach)
        t += delta_t

    return tvalues, Uvalues, Rvalues, Vvalues, v_approach, end_time



Natural_Decay_Time = 26e-9


# Run multiple trials for each initial velocity
v_inits = np.linspace(-0.3, 0.3, 50)  # Range of initial velocities
num_trials = 100  # Number of simulations per initial velocity
final_velocities = []
total_times = []
end_times = []


# Loop over each initial velocity and collect final velocities and lifetimes from multiple trials
for v_init in v_inits:
    for _ in range(num_trials):
        # Run PEC_Dynamics and track lifetime and calculation duration
        times, _, _, _, v_final, calc_duration = PEC_Dynamics_Test(50e-8, v_init, Natural_Decay_Time)
        
        # Append final velocity data
        final_velocities.append((v_init, v_final))
        
        # Track particle lifetime (last recorded time in `times`) in nanoseconds
        if len(times) > 0:
            total_times.append(times[-1] * 1e9)  # Convert to nanoseconds
        else:
            total_times.append(np.nan)  # Handle missing values
        
        # Append the calculation time in seconds
        end_times.append(calc_duration)

# Convert results to numpy array for easier plotting
final_velocities = np.array(final_velocities)
initial_vels = final_velocities[:, 0]
final_vels = final_velocities[:, 1]

# Plot histogram of lifetimes as percentage
plt.hist(total_times, bins=50, weights=np.ones(len(total_times)) / len(total_times) * 100)  # Normalize to percentage
plt.xlabel("Lifetime (ns)")
plt.ylabel("Frequency (%)")
plt.title("Lifetime Distribution as Percentage")
plt.show()

# Plot histogram of calculation times
plt.hist(end_times, bins=50)
plt.xlabel("Calculation Time (s)")
plt.ylabel("Frequency")
plt.title("Distribution of Calculation Durations")
plt.show()

# Filter out high velocities for the 2D histogram
valid_mask = np.abs(final_vels) <= 10

# Apply mask to filter the initial and final velocity arrays
filtered_initial_vels = initial_vels[valid_mask]
filtered_final_vels = final_vels[valid_mask]

# Plot 2D histogram of initial and final velocities (excluding outliers)
plt.figure(dpi=300)
plt.hist2d(filtered_initial_vels, filtered_final_vels, bins=[50, 50], cmap="viridis")
plt.colorbar(label="Count")
plt.xlabel("Initial Velocity (m/s)")
plt.ylabel("Final Velocity (m/s)")
plt.title("Distribution of Final Velocities (Excluding Outliers)")
plt.show()


##### The plots below are for individual samples, it is worth running them a few times to get a sense for the variability.

# Run simulations with initial velocities
dist = 50e-9
t_01, U_01, R_01, V_01, _, _ = PEC_Dynamics_Test(dist, 0.1, Natural_Decay_Time)
t__01, U__01, R__01, V__01, _, _ = PEC_Dynamics_Test(dist, -0.1, Natural_Decay_Time)
t_02, U_02, R_02, V_02, _, _ = PEC_Dynamics_Test(dist, 0.2, Natural_Decay_Time)
t__02, U__02, R__02, V__02, _, _ = PEC_Dynamics_Test(dist, -0.2, Natural_Decay_Time)
t_0, U_0, R_0, V_0, _, _ = PEC_Dynamics_Test(dist, 0, Natural_Decay_Time)

# Plot detuning results
plt.plot(t_01, U_01, label="Initial v = 0.1 m/s")
plt.plot(t__01, U__01, label="Initial v = -0.1 m/s")
plt.plot(t_02, U_02, label="Initial v = 0.2 m/s")
plt.plot(t__02, U__02, label="Initial v = -0.2 m/s")
plt.xlabel('Time (s)')
plt.ylabel('Detuning (MHz)')
plt.legend()
plt.show()

# Plot detuning results
plt.plot(t_01, V_01, label="Initial v = 0.1 m/s")
plt.plot(t__01, V__01, label="Initial v = -0.1 m/s")
plt.plot(t_02, V_02, label="Initial v = 0.2 m/s")
plt.plot(t__02, V__02, label="Initial v = -0.2 m/s")
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.legend()
plt.show()

# Optional: Plot separation distances over time
plt.plot(t_01, R_01, label="Initial v = 0.1 m/s")
plt.plot(t__01, R__01, label="Initial v = -0.1 m/s")
plt.plot(t_02, R_02, label="Initial v = 0.2 m/s")
plt.plot(t__02, R__02, label="Initial v = -0.2 m/s")
plt.xlabel('Time (s)')
plt.ylabel('Separation Distance R (m)')
plt.legend()
plt.title('Time Evolution of Particle Separation')
plt.show()

