#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 16:18:40 2024

@author: joe v
"""

import numpy as np 
import scipy.constants as sc
import matplotlib.pyplot as plt
from Trap_Functions import *
from Potential_Functions import * 
import random
import time

COMvelocities = []
broadened_detunings = []
scattering_rates = []
scattering_rates_red = []
successful_scattering_rates = []
Initial_Positions_Record = []
Initial_Velocities_Record = []
Total_Trapped_Record = []
Detuning_record = []
Doppler_Shifted_Freq_record = []
Scattering_Dist = []
Timestep_Dist = []

def run_simulation(Detuning):
        
    if Detuning == 0:
        Detuning = 1e4 # Actual zero causes sim to fail, 1e4 is only 0.01 MHz 
    
    start_time = datetime.now()
    
    success = 'success'
    
    CurrentTime = 0
    FWHM = 5.75e6
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    Laser_Freq = D1_Freq + Detuning
    Laser_lamda = sc.c /Laser_Freq
    zR_cooling = sc.pi * w0_cooling**2 * Laser_Freq / sc.c
    k = 2 * sc.pi / Laser_Freq
        
    TimeStamps = []
    
    collision_count = 0
    collision_time = []
    scattering_times = []
    collision_coords = []
    scattering_coords = []
    collision_distance = []
    collision_colour = []
    ApproachVelocities = []
    involved_in_collisions = set()
    StepCounter = 1

    # Initialize number of trapped particles
    TotalTrapped = np.random.choice([1, 2], size=1, p=[0.5, 0.5])[0]

    # Initialize positions and velocities
    TweezerPosDist = np.random.normal(0, PosStdDev, (TotalTrapped, 3)) 
    TweezerVelDist = np.array([maxwell_boltzmann_random_velocity(Temp, m) for _ in range(TotalTrapped)])
    
    Initial_Positions_Record.append(TweezerPosDist.copy())
    Initial_Velocities_Record.append(TweezerVelDist.copy())
    Total_Trapped_Record.append(TotalTrapped)

    # Initialize position history arrays for plotting
    PositionHistory = [[] for _ in range(TotalTrapped)]
    VelocityHistory = [[] for _ in range(TotalTrapped)]
    SeparationHistory = []
        
    while CurrentTime < ActualTime:
        
        DeltaT = maxDeltaT
        
        StepCounter = StepCounter +1 
        if StepCounter > 1e7:  
            print(rf"StepCounter capped at {CurrentTime*1e3:.3f}ms . Exiting simulation.")
            end_time = datetime.now() - start_time
            print("Stuck for ",end_time)
            success = 'failure'
            break
        
        if TotalTrapped < 2:
            SeparationHistory.append(np.nan)

        else:
            distances = np.linalg.norm(TweezerPosDist[:, np.newaxis] - TweezerPosDist, axis=2)  # These are the separation distances between pairs of atoms
            Coupling_Potentials = CouplingPotential(distances)
            close_pairs = np.argwhere((distances < 7.5e-4) & (distances > 0))
            
            for pair in close_pairs:
                
                i, j = pair
                
                if i < j:
                    SeparationHistory.append(distances[i, j])
                
                COMvelocity = (TweezerVelDist[i] + TweezerVelDist[j]) / 2
                COMvelocities.append(COMvelocity)
                
                if i < j and Coupling_Potentials[i, j] > gamma_D2:
                    
                    COMvelocity_normal = np.linalg.norm(COMvelocity) # Centre of masss velocity
                    Doppler_Shifted_Freq = doppler_shift(Laser_Freq, COMvelocity_normal)
                    Doppler_Shifted_Freq_record.append(Doppler_Shifted_Freq)
                    Detuning = Doppler_Shifted_Freq - D1_Freq # Detuning from transition frequency
                    
                    Detuning_record.append(Detuning)
                    potential_at_center = potential(0, 0, 0, w0_tweezer, zR_tweezer, U0)/ sc.h  # The tweezer potential at the. entre of the trap
                    detuning_i = Detuning - (potential_at_center - potential(TweezerPosDist[i][0], TweezerPosDist[i][1], TweezerPosDist[i][2], w0_tweezer, zR_tweezer, U0)/ sc.h) # i and j here are the two atoms in a given pair
                    detuning_j = Detuning - (potential_at_center - potential(TweezerPosDist[j][0], TweezerPosDist[j][1], TweezerPosDist[j][2], w0_tweezer, zR_tweezer, U0)/ sc.h)
                    average_detuning = (detuning_i + detuning_j) / 2  # The detuning experienced by an atom depends on where in the trap it lies, so may be different for each
                    
                    avg_position = (TweezerPosDist[i] + TweezerPosDist[j]) / 2
                    D1_scattering = absorption_scattering_D1(avg_position[0], avg_position[1], avg_position[2], zR_cooling, Laser_Freq, P_cooling, w0_cooling, distances[i, j])
                    D2_scattering = absorption_scattering_D2(avg_position[0], avg_position[1], avg_position[2], zR_cooling, Laser_Freq, P_cooling, w0_cooling, distances[i, j])
                    Absorption_Scattering = D1_scattering + D2_scattering # This is for the scattering onto the repulsive (blue) potential
                    
                    D1_scattering_red = absorption_scattering_D1_Red(avg_position[0], avg_position[1], avg_position[2], zR_cooling, Laser_Freq, P_cooling, w0_cooling, distances[i, j])
                    D2_scattering_red = absorption_scattering_D2_Red(avg_position[0], avg_position[1], avg_position[2], zR_cooling, Laser_Freq, P_cooling, w0_cooling, distances[i, j])
                    Absorption_Scattering_red = D1_scattering_red + D2_scattering_red    # This is for the scattering onto the attractive (red) potential        
                    
                    DeltaT = min(0.1 / Absorption_Scattering, 0.1 / Absorption_Scattering_red, maxDeltaT)
                    Timestep_Dist.append(DeltaT)
                    
                    scattering_rates.append(Absorption_Scattering)
                    scattering_rates_red.append(Absorption_Scattering_red)
                                                        
                    if  (len(collision_time) == 0 or CurrentTime - collision_time[-1] > decay_time):  # The time between collisions cannot be shorter than the decay time
                        p = random.uniform(0,1)
                        
                        biggest_scattering = max(Absorption_Scattering, Absorption_Scattering_red)
                        if p < biggest_scattering * DeltaT:  # This probabilistically determines if a scattering event occurs 
                            start_collision = time.time()
                            collision_count += 1
                            collision_time.append(CurrentTime)
                            collision_coords.append([TweezerPosDist[i].copy(), TweezerPosDist[j].copy()])
                            collision_distance.append(distances[i, j])
                            involved_in_collisions.update([i, j])
                            
                            separation_vector = TweezerPosDist[i] - TweezerPosDist[j]
                            distance = np.linalg.norm(separation_vector)
                            separation_unit_vector = separation_vector / distance
                            relative_velocity = TweezerVelDist[i] - TweezerVelDist[j]
                            approach_velocity = np.dot(relative_velocity, separation_unit_vector)
                            ApproachVelocities.append(approach_velocity)
                            
                            old_approach_vector = approach_velocity * separation_unit_vector
                            TweezerVelDist[i] -= old_approach_vector / 2
                            TweezerVelDist[j] += old_approach_vector / 2
                            
                            p = random.uniform(0,1)
                            if p < Absorption_Scattering_red/(Absorption_Scattering+Absorption_Scattering_red): # Determines if the scattering event is 
                                R, v_approach = PEC_Dynamics_Red(distance, approach_velocity, Natural_Decay_Time)
                                collision_colour.append('red')
                            else:
                                R, v_approach = PEC_Dynamics(distance, approach_velocity, Natural_Decay_Time)
                                collision_colour.append('blue')                   
                            
                            a = 0.5 # This determines how the kinetic energy kick is distributed between the atoms.
                            b = 0.5 # Currently set to 50/50, can also set a = random.uniform(0,1), and b = 1 - a
                            
                            new_approach_vector = v_approach * separation_unit_vector
                            TweezerVelDist[i] += new_approach_vector * a
                            TweezerVelDist[j] -= new_approach_vector * b

                            kick = AbsorptionEventKick(k, m) + AbsorptionEventKick(k, m)

                            TweezerPosDist[i] += TweezerVelDist[i] * DeltaT + kick
                            TweezerPosDist[j] += TweezerVelDist[j] * DeltaT + kick
                                             
        for i in range(TotalTrapped):  # Single atom scattering events
               
            D1_single_scattering = absorption_scattering_D1_single_atom(TweezerPosDist[i][0], TweezerPosDist[i][1], TweezerPosDist[i][2], zR_cooling, Laser_Freq, P_cooling, w0_cooling)
            D2_single_scattering = absorption_scattering_D2_single_atom(TweezerPosDist[i][0], TweezerPosDist[i][1], TweezerPosDist[i][2], zR_cooling, Laser_Freq, P_cooling, w0_cooling)
            single_scattering = D1_single_scattering + D2_single_scattering
        
            DeltaT = min(0.1 / single_scattering, maxDeltaT)
                
            if  (len(scattering_times) == 0 or CurrentTime - scattering_times[-1] > decay_time):
                p = random.uniform(0,1)
                if p < single_scattering * DeltaT:
                    scattering_times.append(CurrentTime)
                    scattering_coords.append(TweezerPosDist[i].copy())
                    kick = AbsorptionEventKick(k, m) + AbsorptionEventKick(k, m)
                    TweezerVelDist[i] += np.sqrt(1000) * kick
            
            accel_vec = AccelVec(TweezerPosDist[i], TweezerVelDist[i], beta, U0, zR_tweezer, w0_tweezer, m)
            
            TweezerVelDist[i] += accel_vec * DeltaT
            TweezerPosDist[i] += TweezerVelDist[i] * DeltaT
                
            PositionHistory[i].append(TweezerPosDist[i].copy())
            VelocityHistory[i].append(TweezerVelDist[i].copy())
            
        
        
        CurrentTime += DeltaT
        TimeStamps.append(CurrentTime)
        
        if TotalTrapped >= 2 and len(SeparationHistory) < len(TimeStamps):
            SeparationHistory.append(np.nan)
            
            
        

    # Store final positions
    final_positions = TweezerPosDist.copy()

    within_beam_waist_count = np.sum(np.linalg.norm(final_positions, axis=1) <= 1.5 * w0_tweezer)
    end_time = datetime.now() - start_time
    print(success, end_time)
    
    return within_beam_waist_count, PositionHistory, TotalTrapped, ActualTime, collision_time, collision_coords, involved_in_collisions, collision_count, beta, TimeStamps, VelocityHistory, collision_distance, SeparationHistory, collision_colour, ApproachVelocities, scattering_coords, scattering_times  
