
"""
# Tweezer_Sim_Public
Monte Carlo simulation of tweezer trap dynamics to investigate deterministic single-atom loading

# Monte Carlo Simulation of Atoms in a Tweezer Trap

## Overview

This project implements a Monte Carlo simulation to model the dynamics of atoms in an optical tweezer trap. The primary experimental variable is the cooling laser detuning, and we wish to investigate how it effects the loading efficiency.

## Structure

Master.py: Runs the simulation for a given cooling laser detuning, generating data and plots which help to characterise the simulation.

Scan.py: Runs the simulation over a range of detunings.

Scan_Plot.py: Processes the results from `Scan.py`, plotting how loading efficiency varies with detuning.

Potential_Functions: Contains functions for calculating potential energy curves.

Trap_Functions.py: Defines  functions related to the tweezer trap.

Trap_Simulation: Handles core simulation logic.

"""

