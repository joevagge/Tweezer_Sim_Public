#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:38:09 2024

@author: joe
"""

"""
This script holds most of the functions used in the simulation
"""

import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt
from scipy.special import logsumexp
from Potential_Functions import * 


# Constants
c = sc.c
h = sc.h
kB = sc.k
eps0 = sc.epsilon_0 # [F⋅m−1]

D1_Freq = 377107.463380e9 # Transition frequency 
D2_Freq = 384230.4844685e9 # Transition frequency 

lambdatweezer = 852e-9 # Tweezing laser wavelength [m]
tweezer_freq = sc.c / lambdatweezer

m = 87 * sc.atomic_mass  # Mass of Rubidium 87
gamma_D1 = 36.129e6 /(2*sc.pi)# Hz
gamma_D2 = 38.117e6 /(2*sc.pi) # Hz

decay_time = 1/gamma_D2

Natural_Decay_Time = 26e-9

ActualTime = 8e-3 # Total period of the experiment in seconds

Temp = 20e-6
beta = m/(2*163e-6) # 4.43149996301227e-22
w0_tweezer = 1.24e-6 #[m]  Tweezer waist
w0_cooling = 3e-3 #[m]  Cooling laser waist
zR_tweezer = sc.pi * w0_tweezer**2 / lambdatweezer
P_tweezer = 8e-3 # Laser Power [W]
P_cooling = 5e-3 # Laser Power per beam [W]
I0_tweezer = 2 * P_tweezer / (sc.pi * w0_tweezer**2)
I0_cooling = 2 * P_cooling / (sc.pi * w0_cooling**2)

def U_dip(laser_freq, gamma, I_r, transition_freq):
    prefactor = (3 * sc.pi * sc.c**2) / (2 * (2 * sc.pi*transition_freq)**3)
    detuning_factor = ((gamma  / ((laser_freq - transition_freq)))+(gamma / ((laser_freq + transition_freq))))
    U_dip_value = - prefactor * detuning_factor * I_r
    return U_dip_value

U0_D1 = U_dip(tweezer_freq, gamma_D1 , I0_tweezer, D1_Freq)
U0_D2 = U_dip(tweezer_freq, gamma_D2, I0_tweezer, D2_Freq)
U0 = U0_D1 + U0_D2
Trap_Depth_mK = 1e3 * U0 / sc.k
Trap_Depth_MHz = 1e-6 * U0 / sc.h 
201/61.9715031759222
TrapFreq = np.sqrt(4 * U0 / (m * w0_tweezer**2))
PosStdDev = np.sqrt((sc.k * Temp )/(m * TrapFreq**2 )) 

maxDeltaT = w0_tweezer/(0.25*30)

def AccelVec(PosVec, VelVec, beta, U0, zR, w0, m):
    x, y, z = PosVec
    vxt, vyt, vzt = VelVec

    # Calculate the exponent for the potential term
    exponent = (2 * (x**2 + y**2)) / (w0**2 * (1+(z/zR)**2))
    
    # Clamping the exponent to prevent overflow
    max_exponent = 600
    clamped_exponent = np.clip(exponent, None, max_exponent)
    
    # Compute the exponential term
    exp_term = np.exp(clamped_exponent)

    # Compute the denominator, adding a small value to avoid division by zero
    denominator = exp_term * w0**2 * (1 + (z/zR)**2)**2
    
    zdenominator1 = exp_term * w0**2 * (1 + (z/zR)**2)**3 * zR**2
    zdenominator2 = exp_term * (1 + (z/zR)**2)**2 * zR**2

    # If the denominator is still zero, log and return zeros
    if denominator <= 1e-100:
        print("Warning: Denominator is zero. Adjusting to avoid division by zero.")
        return np.array([0, 0, 0])

    # Compute acceleration components
    accel_x = -(1/m) * (beta * vxt + (4 * U0 * x) / denominator)
    accel_y = -(1/m) * (beta * vyt + (4 * U0 * y) / denominator)
    # accel_z = -(1/m) * (beta * vzt + (4 * U0 * z) / denominator)

    term1 = (4 * U0 * (x**2 + y**2) * z) / zdenominator1  
    term2 = (2 * U0 * z) / zdenominator2                   
    accel_z = -(1/m) * (beta * vzt - term1 + term2)    
    
    return np.array([accel_x, accel_y, accel_z])


def DetuningKEVelocity(Detuning, m):
    phi = np.random.uniform(0, 2 * np.pi)
    theta = np.arccos(np.random.uniform(-1, 1))  
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)
    
    direction = np.array([x, y, z])
    magnitude = np.sqrt(np.abs(2 * h * Detuning / m))
    return magnitude * direction
 
def AbsorptionEventKick(k, m):

    phi = np.random.uniform(0, 2 * sc.pi) 
    cos_theta = np.random.uniform(-1, 1) 
    sin_theta = np.sqrt(1 - cos_theta**2)  
    
    x = sin_theta * np.cos(phi)
    y = sin_theta * np.sin(phi)
    z = cos_theta
    
    magnitude = sc.hbar * k / m  
    
    return magnitude * np.array([x, y, z])



def plot_positions(num, PositionHistory, TotalTrapped, TimeStamps, collision_time, collision_coords, w0, collision_colour, SeparationHistory_list, scattering_times, scattering_coords):
    fig, axs = plt.subplots(4, 1, figsize=(12, 24), dpi = 300)  # 4 subplots now

    titles = [
        rf'({num}) X Position Over Time', 
        rf'({num}) Y Position Over Time', 
        rf'({num}) Z Position Over Time', 
        rf'({num}) Separation Distance Over Time'
    ]
    
    y_labels = ['X Position (μm)', 'Y Position (μm)', 'Z Position (μm)', 'Separation Distance (μm)']

    # Determine x-limit range
    x_min = (collision_time[-1] - 0.00005) if collision_time else 0.0025 - 0.00004
    x_max = (collision_time[-1] + 0.00005) if collision_time else TimeStamps[-1]
    
    # Filter TimeStamps to only include those within the x-limit
    filtered_indices = [i for i, t in enumerate(TimeStamps) if x_min <= t <= x_max]
    filtered_timestamps = [TimeStamps[i] for i in filtered_indices]

    # Loop through X, Y, Z coordinates
    for idx in range(3):  # idx corresponds to x=0, y=1, z=2
        for i in range(TotalTrapped):
            positions_to_plot = [
                PositionHistory[i][t][idx] * 1e6 if np.linalg.norm(PositionHistory[i][t]) <= 200 * w0 else np.nan
                for t in filtered_indices
            ]

            # Plot the particle positions
            axs[idx].plot(filtered_timestamps, positions_to_plot, 
                          'g', label=f'Particle {i+1}')

        # Add collision markers for this axis
        for ct in collision_time:
            if x_min <= ct <= x_max:  # Only plot collisions within x-limit
                collision_idx = collision_time.index(ct)
                ct_index = np.argmin(np.abs(np.array(TimeStamps) - ct))  # Closest index to collision time
                if ct_index in filtered_indices:
                    pos_at_collision = PositionHistory[i][ct_index][idx] * 1e6
                    color = 'ro' if collision_colour[collision_idx] == 'red' else 'bo'
                    axs[idx].plot(ct, pos_at_collision, color)
                    
        # for ct in scattering_times:
        #     if x_min <= ct <= x_max:  # Only plot collisions within x-limit
        #         collision_idx = scattering_times.index(ct)
        #         ct_index = np.argmin(np.abs(np.array(TimeStamps) - ct))  # Closest index to collision time
        #         if ct_index in filtered_indices:
        #             pos_at_collision = PositionHistory[i][ct_index][idx] * 1e6
        #             color = 'yo' 
        #             axs[idx].plot(ct, pos_at_collision, color)

        # Add titles, labels, grid, and y-axis limits
        axs[idx].set_xlabel("Time (ms)")
        axs[idx].set_ylabel(y_labels[idx])
        axs[idx].set_title(titles[idx], fontsize=20, fontweight='bold')
        axs[idx].grid(True)
        axs[idx].set_ylim(-1 * w0 * 1e6, 1 * w0 * 1e6)  
        if len(collision_time) > 0:
            axs[idx].set_xlim((collision_time[-1]) - 0.00002, (collision_time[-1]) + 0.00002)   
        else:
            axs[idx].set_xlim(0, TimeStamps[-1])

    # Plot separation distance
    if SeparationHistory_list:
        separation_to_plot = [SeparationHistory_list[i] * 1e6 for i in filtered_indices]
        axs[3].plot(filtered_timestamps, separation_to_plot, 'k', label='Separation Distance')

        # Add collision markers for the separation distance
        for ct in collision_time:
            if x_min <= ct <= x_max:
                collision_idx = collision_time.index(ct)
                ct_index = np.argmin(np.abs(np.array(TimeStamps) - ct))
                if ct_index in filtered_indices:
                    sep_at_collision = SeparationHistory_list[ct_index] * 1e6
                    color = 'ro' if collision_colour[collision_idx] == 'red' else 'bo'
                    axs[3].plot(ct, sep_at_collision, color)

    axs[3].set_xlabel("Time (ms)")
    axs[3].set_ylabel(y_labels[3])
    axs[3].set_title(titles[3], fontsize=20, fontweight='bold')
    if len(collision_time) > 0:
        axs[3].set_xlim((collision_time[-1]) - 0.00002, (collision_time[-1]) + 0.00002)   
    else:
        axs[3].set_xlim(0, TimeStamps[-1])
    axs[3].grid(True)

    # Adjust layout
    plt.tight_layout()
    plt.show()
    
def plot_Separation(num, TimeStamps, SeparationHistory_list):
    # Create a single plot
    fig, ax = plt.subplots(1, 1, figsize=(12, 8), dpi=300)
    
    # Ensure SeparationHistory_list is not empty
    if not SeparationHistory_list:
        raise ValueError("SeparationHistory_list is empty. Cannot plot separation.")
    if len(SeparationHistory_list) > 1:
        raise ValueError("Can't plot separations for multiple particles")

    # Multiply by 1e6 to convert to micrometers
    separation_to_plot = [sep * 1e6 for sep in SeparationHistory_list]
    
    # Plot the separation data
    ax.plot(TimeStamps, separation_to_plot, 'k', label='Separation Distance')

    # Set labels, title, and grid
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Separation Distance (μm)")
    ax.set_title(f"({num}) Separation Distance Over Time", fontsize=20, fontweight='bold')
    ax.set_xlim(0, TimeStamps[-1])
    ax.grid(True)

    # Adjust layout for better appearance
    plt.tight_layout()
    plt.show()


def plot_velocities(num, VelocityHistory, TotalTrapped, TimeStamps, collision_time, collision_coords, w0, collision_colour):
    fig, axs = plt.subplots(3, 1, figsize=(12, 18))  # 3 subplots for X, Y, Z velocities

    titles = [
        rf'({num}) X Velocity Over Time', 
        rf'({num}) Y Velocity Over Time', 
        rf'({num}) Z Velocity Over Time'
    ]
    y_labels = ['X Velocity (μm/s)', 'Y Velocity (μm/s)', 'Z Velocity (μm/s)']

    # Determine x-limit range
    x_min = (collision_time[-1] - 0.00002) if len(collision_time) > 0 else 0
    x_max = (collision_time[-1] + 0.00002) if len(collision_time) > 0 else TimeStamps[-1]

    # Filter TimeStamps and create indices for filtering
    filtered_indices = [i for i, t in enumerate(TimeStamps) if x_min <= t <= x_max]
    filtered_timestamps = [TimeStamps[i] for i in filtered_indices]

    for idx in range(3):  # Loop for x, y, z velocities
        for i in range(TotalTrapped):
            velocities_to_plot = [
                VelocityHistory[i][t][idx] * 1e6 if np.linalg.norm(VelocityHistory[i][t]) <= 60 * w0 else np.nan
                for t in filtered_indices
            ]

            axs[idx].plot(filtered_timestamps, velocities_to_plot, 'b', label=f'Particle {i+1}')

        # Add collision markers
        for ct in collision_time:
            if x_min <= ct <= x_max:  # Only plot collisions within x-limit
                collision_idx = collision_time.index(ct)
                ct_index = np.argmin(np.abs(np.array(TimeStamps) - ct))  # Closest timestamp to collision
                if ct_index in filtered_indices:  # Ensure the collision index is valid
                    velocity_at_collision = VelocityHistory[i][ct_index][idx] * 1e6
                    color = 'ro' if collision_colour[collision_idx] == 'red' else 'bo'
                    axs[idx].plot(ct, velocity_at_collision, color)

        axs[idx].set_xlabel("Time (ms)")
        axs[idx].set_ylabel(y_labels[idx])
        axs[idx].set_title(titles[idx], fontsize=20, fontweight='bold')
        axs[idx].grid(True)
        # axs[idx].set_ylim(-60 * w0 * 1e6, 60 * w0 * 1e6)

    plt.tight_layout()
    plt.show()



def maxwell_boltzmann_random_velocity(T, m):
    v = np.sqrt(2 * kB * T / m) * np.random.randn(3)
    return v

def potential(x, y, z, w0, zR, U0):
    w_z = w0 * np.sqrt(1 + (z / zR)**2)
    intensity = (w0 / w_z)**2 * np.exp(-2 * (x**2 + y**2) / w_z**2)
    return -U0 * intensity

def OtherScatteringRate_D1(laser_freq, gamma, P, w0, R):
    prefactor = (3 * sc.pi * sc.c**2) / (2 * 2 * sc.pi*D1_Freq**3 * h / (2 * sc.pi))
    lambda_ratio = (laser_freq / D1_Freq)**3
    detuning_factor = (gamma_D1 / (SPPotential(transition_freq, CouplingPotential(R)) - laser_freq) +
                       gamma_D1 / (SPPotential(transition_freq, CouplingPotential(R)) + laser_freq))**2
    return prefactor * lambda_ratio * (2 * P) / (sc.pi * w0**2) * detuning_factor


#### These functions describe the scattering during light assisted collisions ###

def absorption_scattering_D1(x, y, z, zR, laser_freq, P, w0, R):
    wz = w0 * np.sqrt(1 + (z / zR)**2)
    I_sat = (sc.pi * sc.h * sc.c * gamma_D1 * D1_Freq **3) / (3 * sc.c**3) 
    # print(rf'I_sat ={I_sat:.1e}')
    I0 = 2 * P / (sc.pi * w0**2)  
    # print(rf'I0 = {I0:.1e}')
    I = I0 * (w0 / wz)**2 * np.exp((-2 * (x**2 + y**2)) / wz**2)    
    s = I / I_sat  
    # print(rf's 2 = {s:.1e}')
    Delta = laser_freq - SPPotential(D1_Freq, CouplingPotential(R)) 
    # print(rf'Delta 1 Blue= {Delta:.1e}')
    R_abs = (gamma_D1/2) * ( s / (1 + 4*(Delta/gamma_D1)**2))
    R_abs = R_abs * 3 * 2 # (for retro reflected beams)
    return R_abs

def absorption_scattering_D2(x, y, z, zR, laser_freq, P, w0, R):
    wz = w0 * np.sqrt(1 + (z / zR)**2)    
    I_sat = (sc.pi * sc.h * sc.c * gamma_D2 * D2_Freq **3) / (3 * sc.c**3) 
    # print(rf'I_sat ={I_sat:.1e}')
    I0 = 2 * P / (sc.pi * w0**2)    
    # print(rf'I0 = {I0:.1e}')
    I = I0 * (w0 / wz)**2 * np.exp((-2 * (x**2 + y**2)) / wz**2)    
    s = I / I_sat   
    # print(rf's 2 = {s:.1e}')
    Delta = laser_freq - SPPotential(D2_Freq, CouplingPotential(R))  
    R_abs = (gamma_D2/2) * ( s / (1 + 4*(Delta/gamma_D2)**2))
    R_abs = R_abs * 3 * 2 # (for retro reflected beams)
    return R_abs

def absorption_scattering_D1_Red(x, y, z, zR, laser_freq, P, w0, R):
    wz = w0 * np.sqrt(1 + (z / zR)**2)    
    I_sat = (sc.pi * sc.h * sc.c * gamma_D1 * D1_Freq **3) / (3 * sc.c**3) 
    # print(rf'I_sat ={I_sat:.1e}')
    I0 = 2 * P / (sc.pi * w0**2)    
    # print(rf'I0 = {I0:.1e}')
    I = I0 * (w0 / wz)**2 * np.exp((-2 * (x**2 + y**2)) / wz**2)    
    s = I / I_sat   
    # print(rf's 2 = {s:.1e}')
    Delta = laser_freq - PSPotential(D1_Freq, CouplingPotential(R))  
    # print(rf'Delta 1 Red= {Delta:.1e}')
    R_abs = (gamma_D1/2) * ( s / (1 + 4*(Delta/gamma_D1)**2))
    R_abs = R_abs * 3 * 2 # (for retro reflected beams)
    return R_abs 

def absorption_scattering_D2_Red(x, y, z, zR, laser_freq, P, w0, R):
    wz = w0 * np.sqrt(1 + (z / zR)**2)    
    I_sat = (sc.pi * sc.h * sc.c * gamma_D2 * D2_Freq **3) / (3 * sc.c**3) 
    # print(rf'I_sat ={I_sat:.1e}')
    I0 = 2 * P / (sc.pi * w0**2)    
    # print(rf'I0 = {I0:.1e}')
    I = I0 * (w0 / wz)**2 * np.exp((-2 * (x**2 + y**2)) / wz**2)    
    s = I / I_sat   
    # print(rf's 2 = {s:.1e}')
    Delta = laser_freq - PSPotential(D2_Freq, CouplingPotential(R))  
    # print(rf'Delta 2 = {Delta:.1e}')
    R_abs = (gamma_D2/2) * ( s / (1 + 4*(Delta/gamma_D2)**2))
    R_abs = R_abs * 3 * 2 # (for retro reflected beams)
    return R_abs

#### These 'single atom' ones are for regular scattering events ###

def absorption_scattering_D2_single_atom(x, y, z, zR, laser_freq, P, w0):
    wz = w0 * np.sqrt(1 + (z / zR)**2)    
    I_sat = (sc.pi * sc.h * sc.c * gamma_D2 * D2_Freq **3) / (3 * sc.c**3) 
    # print(rf'I_sat ={I_sat:.1e}')
    I0 = 2 * P / (sc.pi * w0**2)    
    # print(rf'I0 = {I0:.1e}')
    I = I0 * (w0 / wz)**2 * np.exp((-2 * (x**2 + y**2)) / wz**2)    
    s = I / I_sat   
    # print(rf's 2 = {s:.1e}')
    Delta = laser_freq - D2_Freq
    # print(rf'Delta 2 = {Delta:.1e}')
    R_abs = (gamma_D2/2) * ( s / (1 + 4*(Delta/gamma_D2)**2))
    R_abs = R_abs * 3 * 2 # (for retro reflected beams)
    return R_abs / 1000

def absorption_scattering_D1_single_atom(x, y, z, zR, laser_freq, P, w0):
    wz = w0 * np.sqrt(1 + (z / zR)**2)    
    I_sat = (sc.pi * sc.h * sc.c * gamma_D1 * D1_Freq **3) / (3 * sc.c**3) 
    # print(rf'I_sat ={I_sat:.1e}')
    I0 = 2 * P / (sc.pi * w0**2)    
    # print(rf'I0 = {I0:.1e}')
    I = I0 * (w0 / wz)**2 * np.exp((-2 * (x**2 + y**2)) / wz**2)    
    s = I / I_sat   
    # print(rf's 2 = {s:.1e}')
    Delta = laser_freq - D1_Freq
    # print(rf'Delta 2 = {Delta:.1e}')
    R_abs = (gamma_D1/2) * ( s / (1 + 4*(Delta/gamma_D1)**2))
    R_abs = R_abs * 3 * 2 # (for retro reflected beams)
    return R_abs / 1000

def doppler_shift(freq, velocity):
    return freq * np.sqrt((1 + velocity/sc.c) / (1 - velocity/sc.c))