#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 11:09:10 2024

@author: joe
"""
#### Old kick ######

                       if p < biggest_scattering * DeltaT:
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
                           
                           Detuning_Guess = SPPotential(Delta, CouplingPotential(distance)) - SPPotential(Delta, CouplingPotential(1e-3))
                           # print(rf'Detuning_Guess = {average_detuning*1e-6:.2f}MHz')
                           
                           Slowing_Time, New_Delta_Hz, New_Distance, New_Detuning_Hz  = Hill_Kick_Time(approach_velocity, distance)
                           
                           Time_Left = Natural_Decay_Time - Slowing_Time
                           
                           Final_Separation = Hill_Roll_Distance(New_Distance, Time_Left)
                           
                           Final_Detuning = SPPotential(Delta, CouplingPotential(New_Distance)) - SPPotential(Delta, CouplingPotential(Final_Separation))
                           
                           # print(rf'Final_Separation = {Final_Separation*1:.2e}m')
                           # print(rf'average_detuning = {average_detuning*1e-6:.2f}MHz')
                           # print(rf'Top Detuning = {(SPPotential(Delta, CouplingPotential(New_Distance))-Delta)*1e-6:.2f}MHz')
                           # print(rf'Bottom Detuning = {(SPPotential(Delta, CouplingPotential(Final_Separation))-Delta)*1e-6:.2f}MHz')
                           # # print(rf'difference = {difference*1e-6:.2f}MHz')
                           # print(rf'Final_Detuning = {Final_Detuning*1e-6:.2f}MHz')
                           
                           successful_scattering_rates.append(Absorption_Scattering)
                           broadened_detunings.append(average_detuning)
                           
                           
                           a = random.uniform(0,1) #0.5
                           b = 1 - a #0.5
                           
                           kick = OneDAbsorptionKick(k, m) + RandomEmissionKick(k, m) + DetuningKEVelocity(average_detuning, m)
                           
                           p = random.uniform(0,1)
                           if p < Absorption_Scattering_red/(Absorption_Scattering+Absorption_Scattering_red):
                               a = 1000
                               b = 1000
                               collision_colour.append('red')
                               
                           else:
                               collision_colour.append('blue')

                           
                           TweezerVelDist[i] += (a * kick)
                           TweezerVelDist[j] += (b * kick)
                           
                           TweezerPosDist[i] += TweezerVelDist[i] * DeltaT
                           TweezerPosDist[j] += TweezerVelDist[j] * DeltaT
                           
                           
                           
                           
                           
                           
                           
                         
#### Rolling up the hill Calculation ####


def Hill_Kick_Time(v_app, sep_distance):
    
    R_1 = sep_distance
    
    if v_app >= 0:
        
        U_R1 = SPPotential(Delta, CouplingPotential(R_1)) * sc.h
        U_R2 = U_R1 + 0.5 * m * v_app**2
        
        R_2 = R_from_SPPotential(U_R2, Delta)
        # print(rf'R_2 = {R_2*1:.2e}m')
        r = abs(R_1 - R_2)
        
    else:
        
        U_R1 = SPPotential(Delta, CouplingPotential(R_1)) * sc.h
        U_R2 = U_R1 - 0.5 * m * v_app**2
        
        R_2 = R_from_SPPotential(U_R2, Delta)
        r = abs(R_1 - R_2)
        
    UR2_Hz = U_R2 / sc.h
    New_Detuning_Hz = (UR2_Hz - Delta) # This is the detuning after rolling up the hill
            
    R_2 = R_from_SPPotential(U_R2, Delta)
    
    # print('r1',R_1)
    # print(rf'v_app = {v_app*1:.2f}m/s')
    # print(rf'New_Distance = {R_2*1:.2e}m')
    # print(rf'New_Detuning_Hz = {New_Detuning_Hz*1:.2e}Hz')
    
    t = abs((2*r)/0.4)
            
    return t, New_Detuning_Hz, R_2, New_Detuning_Hz

def Hill_Roll_Distance(R_start, time_left):

    final_separation = np.linspace(2e-9,2e-4,5000)
    Timescales = []

    # Calculate T
    for R_finish in final_separation:
        T = calculate_T(R_start, R_finish)
        Timescales.append(T)
        
    # Find the closest T value to target_T
    closest_index = np.argmin(np.abs(np.array(Timescales) - time_left))
    R_finish_closest = final_separation[closest_index]
    T_closest = Timescales[closest_index]
    
    # print(rf'R_start = {R_start*1e9:.2f}nm')
    # print(rf'time_left = {time_left*1e9:.2f}ns')
    # print(rf'T_closest = {T_closest*1e9:.2f}ns')
    # print(rf'R_finish_closest = {R_finish_closest*1e6:.2f}um')
    
    return R_finish_closest


Slowing_Time, New_Delta_Hz, New_Distance, New_Detuning_Hz  = Hill_Kick_Time(-0.04, 50e-9)
Time_Left = Natural_Decay_Time - Slowing_Time

Hill_Roll_Distance(New_Distance, Time_Left)

(SPPotential(Delta, CouplingPotential(3.44e-08)) - SPPotential(Delta, CouplingPotential(1e-3)))*1e-6

(SPPotential(Delta, CouplingPotential(Hill_Roll_Distance(New_Distance, Time_Left))) - Delta)*1e-6




distance = 50e-9
approach_velocities = np.linspace(-0.4,0.4,300)

final_detunings = []

for approach_velocity in approach_velocities:

    Detuning_Guess = SPPotential(Delta, CouplingPotential(distance)) - SPPotential(Delta, CouplingPotential(1e-3))
    
    Slowing_Time, New_Delta_Hz, New_Distance, New_Detuning_Hz  = Hill_Kick_Time(approach_velocity, distance)
    
    Time_Left = Natural_Decay_Time - Slowing_Time
    
    Final_Separation = Hill_Roll_Distance(New_Distance, Time_Left)
    
    Final_Detuning = SPPotential(Delta, CouplingPotential(New_Distance)) - SPPotential(Delta, CouplingPotential(Final_Separation))
    
    final_detunings.append(Final_Detuning)
    

initial_detunings = np.zeros_like(final_detunings) + SPPotential(Delta, CouplingPotential(50e-9)) - SPPotential(Delta, CouplingPotential(1e-2))

plt.figure(figsize=(10, 6))
plt.plot(approach_velocities, final_detunings, color='blue', linestyle='-', linewidth=2, label='Final Detuning')
plt.plot(approach_velocities, initial_detunings, color='orange', linestyle='--', linewidth=2, label='Initial Detuning')

# Add labels and title
plt.xlabel("Approach Velocity (m/s)", fontsize=14)
plt.ylabel("Detuning (Hz)", fontsize=14)
plt.title("Detuning vs. Approach Velocity", fontsize=16)

# Add grid, legend, and make layout tighter
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=12)
plt.tight_layout()                         


