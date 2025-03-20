#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 12:49:06 2024
"""

import numpy as np
import matplotlib.pyplot as plt
from Trap_Simulation import *  
import time
from multiprocessing import Pool, Manager
import os
import sys

"""
This script runs the simulation over a range of laser detunings.

Key Parameters:
- Detunings: Array of detuning values to simulate.
- num_runs: Number of iterations per detuning.

Outputs:
- Data: Saved to a .npz file (e.g., `Scan_Data/neg200to400.npz`).
"""


file_name = 'Scan_Data/neg200to400.npz' 

if os.path.exists(file_name):
    print(f"Error: '{file_name}' already exists. Choose another file name.")
    sys.exit(1)  

Detunings = np.arange(-200, 400e6, 5e6)
num_runs = 250

def run_simulation_for_params(params, progress_counter):
    Detuning, num_runs = params
    results = []
    collision_counts = []

    for i in range(num_runs):
        within_beam_waist_count, _, _, _, _, _, _, collision_count, _, _, _, _, _, _, _, _, _ = run_simulation(Detuning)
        results.append(within_beam_waist_count)
        collision_counts.append(collision_count)  # Store collision counts
        progress_counter.value += 1  # Directly increment without locking
    
    return results, collision_counts

def wrapper_run_simulation(params_and_counter):
    params, progress_counter = params_and_counter
    return run_simulation_for_params(params, progress_counter)

def main():
    manager = Manager()
    progress_counter = manager.Value('i', 0)
    
    params = [(Detuning, num_runs) for Detuning in Detunings] 
    total_simulations = len(params) * num_runs

    start_time = time.time()
    last_printed = 0  

    results_list = []
    for param in params:
        results, collision_counts = run_simulation_for_params(param, progress_counter)
        results_list.append((results, collision_counts))

        completed_simulations = progress_counter.value
        if completed_simulations - last_printed >= 1: 
            elapsed_time = time.time() - start_time
            remaining_simulations = total_simulations - completed_simulations
            average_time_per_simulation = elapsed_time / completed_simulations if completed_simulations > 0 else 0
            estimated_time_remaining = average_time_per_simulation * remaining_simulations
            
            finish_timestamp = time.time() + estimated_time_remaining
            finish_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(finish_timestamp))
            
            hours, remainder = divmod(estimated_time_remaining, 3600)
            minutes, seconds = divmod(remainder, 60)
            
            print(f"Completed {completed_simulations}/{total_simulations} simulations")
            print(f"Estimated Finish Time: {finish_time}")
            print(f"Time Remaining: {int(hours)}:{int(minutes)}:{int(seconds)}")
            
            last_printed = completed_simulations  

    no_atom_likelihoods = []
    one_atom_likelihoods = []
    two_atom_likelihoods = []
    mean_occupations = []
    average_number_of_collisions = []

    for Detuning, (results, collision_counts) in zip(Detunings, results_list):
        no_atom_likelihood = (results.count(0) / num_runs) * 100
        one_atom_likelihood = (results.count(1) / num_runs) * 100
        two_atom_likelihood = (results.count(2) / num_runs) * 100
        mean_occupation = np.mean(results)
        no_atom_likelihoods.append(no_atom_likelihood)
        one_atom_likelihoods.append(one_atom_likelihood)
        two_atom_likelihoods.append(two_atom_likelihood)
        mean_occupations.append(mean_occupation)

        if isinstance(collision_counts[0], (list, np.ndarray)):
            collision_counts_flat = [item for sublist in collision_counts for item in sublist]
        else:
            collision_counts_flat = collision_counts

        mean_collisions = np.mean(collision_counts_flat)
        average_number_of_collisions.append(mean_collisions)

    Detunings_arr = np.array(Detunings)
    no_atom_likelihoods_arr = np.array(no_atom_likelihoods)
    one_atom_likelihoods_arr = np.array(one_atom_likelihoods)
    two_atom_likelihoods_arr = np.array(two_atom_likelihoods)
    mean_occupations_arr = np.array(mean_occupations)
    average_number_of_collisions_arr = np.array(average_number_of_collisions)
    
    param_names = [
        "Iterations",
        "Initial Particle Temp",
        "Trap Depth",
        "Trap Width",
        "Time Span",
        "Initial Position StdDev",
        "Damping Coefficient"
    ]

    param_values = [
        f"{num_runs}",
        f"{Temp * 1e6:.2f} $\mu$K",
        f"{np.round(Trap_Depth_mK, 3)} mK / {np.round(Trap_Depth_MHz, 3)} MHz",
        f"{w0_tweezer * 1e6:.2f} $\mu$m",
        f"{ActualTime*1e3} ms",
        f"{PosStdDev * 1e6:.2f} $\mu$m",
        f"{beta:.0e} Kg/s ({(m/(beta+1e-30))*1e6:.0f} $\mu$s)"
    ]

    totaltime = time.time() - start_time
    hours, remainder = divmod(totaltime, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"Time Taken: {int(hours)}:{int(minutes)}:{int(seconds)}")
    
    np.savez(file_name,
         Detunings=Detunings_arr,
         no_atom_likelihoods=no_atom_likelihoods_arr,
         one_atom_likelihoods=one_atom_likelihoods_arr,
         two_atom_likelihoods=two_atom_likelihoods_arr,
         mean_occupations=mean_occupations_arr,
         average_number_of_collisions=average_number_of_collisions_arr,
         param_names=param_names,
         param_values=param_values,
         U0=U0,
         beta = beta)

if __name__ == '__main__':
    main()
