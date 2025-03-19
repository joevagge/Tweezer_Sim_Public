#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 11:32:12 2024

@author: joe
"""

import numpy as np
import scipy.constants as sc
import matplotlib.pyplot as plt

##### Load Red Data #########

red_data = np.load('Scan_Data/1608_Red_long_163e6_5050.npz')
red_Detunings_arr = red_data['Detunings']
red_no_atom_likelihoods_arr = red_data['no_atom_likelihoods']
red_one_atom_likelihoods_arr = red_data['one_atom_likelihoods']
red_two_atom_likelihoods_arr = red_data['two_atom_likelihoods']
red_mean_occupations_arr = red_data['mean_occupations']
red_average_number_of_collisions_arr = red_data['average_number_of_collisions']
red_param_names = red_data['param_names']
red_param_values = red_data['param_values']
red_U0 = red_data['U0']

##### Load Blue Data 1 ######### 0612_5050

data1 = np.load('Scan_Data/0302_neg150to30.npz')
Detunings_arr1 = data1['Detunings']
no_atom_likelihoods_arr1 = data1['no_atom_likelihoods']
one_atom_likelihoods_arr1 = data1['one_atom_likelihoods']
two_atom_likelihoods_arr1 = data1['two_atom_likelihoods']
mean_occupations_arr1 = data1['mean_occupations']
average_number_of_collisions_arr1 = data1['average_number_of_collisions']
param_names1 = data1['param_names']
param_values1 = data1['param_values']
U01 = data1['U0']

##### Load Blue Data 2 #########

data2 = np.load('Scan_Data/1112_5050.npz')
Detunings_arr2 = data2['Detunings']
no_atom_likelihoods_arr2 = data2['no_atom_likelihoods']
one_atom_likelihoods_arr2 = data2['one_atom_likelihoods']
two_atom_likelihoods_arr2 = data2['two_atom_likelihoods']
mean_occupations_arr2 = data2['mean_occupations']
average_number_of_collisions_arr2 = data2['average_number_of_collisions']
param_names2 = data2['param_names']
param_values2 = data2['param_values']
U02 = data2['U0']

##### Load Blue Data 3 #########

data3 = np.load('Scan_Data/1012_5050.npz')
Detunings_arr3 = data3['Detunings']
no_atom_likelihoods_arr3 = data3['no_atom_likelihoods']
one_atom_likelihoods_arr3 = data3['one_atom_likelihoods']
two_atom_likelihoods_arr3 = data3['two_atom_likelihoods']
mean_occupations_arr3 = data3['mean_occupations']
average_number_of_collisions_arr3 = data3['average_number_of_collisions']
param_names3 = data3['param_names']
param_values3 = data3['param_values']
U03 = data3['U0']


##### Load Blue Data 4 #########

data4 = np.load('Scan_Data/2108_line_broad_163e7.npz')
Detunings_arr4 = data4['Detunings']
no_atom_likelihoods_arr4 = data4['no_atom_likelihoods']
one_atom_likelihoods_arr4 = data4['one_atom_likelihoods']
two_atom_likelihoods_arr4 = data4['two_atom_likelihoods']
mean_occupations_arr4 = data4['mean_occupations']
average_number_of_collisions_arr4 = data4['average_number_of_collisions']
param_names4 = data4['param_names']
param_values4 = data4['param_values']
U04 = data4['U0']

##### Load Joint Data #########

Detunings_arrj = np.concatenate((Detunings_arr1, Detunings_arr2, Detunings_arr3))
no_atom_likelihoods_arrj = np.concatenate((no_atom_likelihoods_arr1, no_atom_likelihoods_arr2, no_atom_likelihoods_arr3))
one_atom_likelihoods_arrj = np.concatenate((one_atom_likelihoods_arr1, one_atom_likelihoods_arr2, one_atom_likelihoods_arr3))
two_atom_likelihoods_arrj = np.concatenate((two_atom_likelihoods_arr1, two_atom_likelihoods_arr2, two_atom_likelihoods_arr3))
mean_occupations_arrj = np.concatenate((mean_occupations_arr1, mean_occupations_arr2, mean_occupations_arr3))
average_number_of_collisions_arrj = np.concatenate((average_number_of_collisions_arr1, average_number_of_collisions_arr2, average_number_of_collisions_arr3))
param_namesj = param_names1 
param_valuesj = param_values1  
U0j = U01

# Detunings_arrj = np.concatenate((Detunings_arr2, Detunings_arr3))
# no_atom_likelihoods_arrj = np.concatenate((no_atom_likelihoods_arr2, no_atom_likelihoods_arr3))
# one_atom_likelihoods_arrj = np.concatenate((one_atom_likelihoods_arr2, one_atom_likelihoods_arr3))
# two_atom_likelihoods_arrj = np.concatenate((two_atom_likelihoods_arr2, two_atom_likelihoods_arr3))
# mean_occupations_arrj = np.concatenate((mean_occupations_arr2, mean_occupations_arr3))
# average_number_of_collisions_arrj = np.concatenate((average_number_of_collisions_arr2, average_number_of_collisions_arr3))
# param_namesj = param_names2
# param_valuesj = param_values2  
# U0j = U02

# Determine the maximum length of the parameter names
max_name_length = max(len(name) for name in param_names2)

# Create formatted parameter strings with centered colons
parameters = [
    f"{name.ljust(max_name_length)} : {value}"
    for name, value in zip(param_names2, param_values2)
]

# Create the parameter text with monospaced font
param_text = "\n".join(parameters)

# # Create the first figure with the original subplots
# fig1, axs1 = plt.subplots(2, 1, figsize=(14, 14), dpi=100)

# # Plot 1: Single Atom Occupation Rate
# # axs1[0].plot(Detunings_arr1 / (U01 / sc.h), one_atom_likelihoods_arr1, 'o-',label = 'Blue Detuned 1')
# # axs1[0].plot(Detunings_arr2 / (U02 / sc.h), one_atom_likelihoods_arr2, 'o-',label = 'Blue Detuned 2')
# axs1[0].plot(Detunings_arr3 / (U03 / sc.h), one_atom_likelihoods_arr3, 'o-',label = 'Blue Detuned 3')
# axs1[0].plot(-red_Detunings_arr / (red_U0 / sc.h), red_one_atom_likelihoods_arr, 'o-', color = 'lightcoral',label = 'Red Detuned')
# axs1[0].set_xlabel('Laser Detuning (Units of Trap Depth)', fontsize=20)
# axs1[0].set_ylabel('Likelihood of Getting One Atom in Trap (%)', fontsize=20)
# axs1[0].set_title('Single Atom Occupation Rate', fontsize=20, fontweight='bold')
# axs1[0].legend(loc='lower right', fontsize=20)
# axs1[0].grid(True)
# axs1[0].tick_params(axis='both', labelsize=16)
# axs1[0].set_xlim((Detunings_arr1[0] / (U01 / sc.h)),Detunings_arr1[-1] / (U01 / sc.h))
# axs1[0].set_ylim(0,100)

# # Plot 2: Mean Trap Occupation
# # axs1[1].plot(Detunings_arr1 / (U01 / sc.h), mean_occupations_arr1, 'o-')
# # axs1[1].plot(Detunings_arr2 / (U02 / sc.h), mean_occupations_arr2, 'o-')
# axs1[1].plot(Detunings_arr3 / (U03 / sc.h), mean_occupations_arr3, 'o-')
# axs1[1].plot(-red_Detunings_arr / (red_U0 / sc.h), red_mean_occupations_arr, 'o-', color = 'lightcoral')
# axs1[1].set_xlabel('Laser Detuning (Units of Trap Depth)', fontsize=20)
# axs1[1].set_ylabel('Mean Trap Occupation', fontsize=20)
# axs1[1].set_title('Mean Trap Occupation', fontsize=20, fontweight='bold')
# axs1[1].grid(True)
# axs1[1].set_xlim((Detunings_arr1[0] / (U01 / sc.h)),Detunings_arr1[-1] / (U01 / sc.h))
# axs1[1].set_ylim(0,1)

# # Add text with key parameters to the second subplot
# axs1[1].text(
#     0.31, 0.5, param_text,
#     fontsize=18,
#     verticalalignment='top',
#     horizontalalignment='left',
#     bbox=dict(
#         facecolor='white',
#         alpha=0.9,
#         edgecolor='black',
#         boxstyle='round,pad=1'
#     ),
#     transform=axs1[1].transAxes,
#     fontfamily='monospace'  # Ensures text alignment is maintained
# )

# # Adjust layout to make space for text box
# plt.tight_layout()
# plt.subplots_adjust(right=0.85)  # Adjust the right margin to make space for text
# plt.xticks(fontsize=16)
# plt.yticks(fontsize=16)

# # Show the first plot
# plt.show()


####### Individual Plot #########
fig1, axs1 = plt.subplots(1, 1, figsize=(11, 6), dpi=400)

# axs1.plot(Detunings_arr / (U0 / sc.h), one_atom_likelihoods_arr, 'o-',label = 'Blue Detuned')
# # axs1.plot(-red_Detunings_arr / (red_U0 / sc.h), red_one_atom_likelihoods_arr, 'o-', color = 'lightcoral',label = 'Red Detuned')
# axs1.set_xlabel('Laser Detuning (Units of Trap Depth)', fontsize=20)
# axs1.set_ylabel('Likelihood of Getting One Atom in Trap (%)', fontsize=20)
# axs1.set_title('Single Atom Occupation Rate', fontsize=20, fontweight='bold')
# axs1.legend(loc='lower right', fontsize=20)
# axs1.grid(True)
# axs1.tick_params(axis='both', labelsize=16)

axs1.plot(Detunings_arr1 / (U01 / sc.h), one_atom_likelihoods_arr1, 'o-',label = 'Blue Detuned corrected')
# axs1.plot(Detunings_arrj / (U0j / sc.h), one_atom_likelihoods_arrj, 'o-',label = 'Blue Detuned')
# axs1.plot(Detunings_arr1 / (U01 / sc.h), one_atom_likelihoods_arr1, 'o-', label=r"$\tau=663\,\mu\mathrm{s}$")
# axs1.plot(Detunings_arr2 / (U02 / sc.h), one_atom_likelihoods_arr2, 'o-',label = r"$\tau=163\,\mu\mathrm{s}$")
# axs1.plot(Detunings_arr3 / (U03 / sc.h), one_atom_likelihoods_arr3, 'o-',label = 'Blue Detuned')
# axs1.plot(Detunings_arr4 / (U04 / sc.h), one_atom_likelihoods_arr4, 'o-',label = r"$\tau=16.3\,\mu\mathrm{s}$")
# axs1.plot(-red_Detunings_arr / (red_U0 / sc.h), red_two_atom_likelihoods_arr, 'o-', color = 'lightcoral',label = 'Red Detuned')
axs1.set_xlabel('Laser Detuning (Units of Trap Depth)', fontsize=20)
axs1.set_ylabel('Likelihood of Getting Two Atoms in Trap (%)', fontsize=17)
# axs1.set_ylabel('Likelihood of Getting One Atom in Trap (%)', fontsize=20)
# axs1.set_title('Average No. of Collisions', fontsize=20, fontweight='bold')
# axs1.set_ylim(80,100)
# axs1.legend(fontsize=20)
# axs1.set_xlim((Detunings_arrj[0] / (U0j / sc.h)),Detunings_arrj[-1] / (U0j / sc.h))
# axs1.set_xlim(0.5,1.3)
axs1.grid(True)


# Adjust layout to make space for text box
plt.tight_layout()
plt.subplots_adjust(right=0.85)  # Adjust the right margin to make space for text
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.gca().set_facecolor('none')

# Show the first plot
plt.show()




# Create the second figure with the new subplots
fig2, axs2 = plt.subplots(4, 1, figsize=(11, 16))

# Plot 1: Zero Atom Occupation Rate
# axs2[0].plot(Detunings_arrj / (U0j / sc.h), no_atom_likelihoods_arrj, 'o-')
axs2[0].plot(Detunings_arr1 / (U01 / sc.h), no_atom_likelihoods_arr1, 'o-')
# axs2[0].plot(Detunings_arr2 / (U02 / sc.h), no_atom_likelihoods_arr2, 'o-')
# axs2[0].plot(Detunings_arr3 / (U03 / sc.h), no_atom_likelihoods_arr3, 'o-')
# axs2[0].plot(-red_Detunings_arr / (red_U0 / sc.h), red_no_atom_likelihoods_arr, 'o-', color = 'lightcoral')
# axs2[0].set_xlabel('Laser Detuning (Units of Trap Depth)', fontsize=20)
# axs2[0].set_ylabel('Likelihood of Getting No Atoms in Trap (%)', fontsize=16)
axs2[0].set_title('Zero Atom Occupation Rate', fontsize=20, fontweight='bold')
axs2[0].grid(True)

# Plot 2: One Atom Occupation Rate
# axs2[1].plot(Detunings_arrj / (U0j / sc.h), one_atom_likelihoods_arrj, 'o-')
axs2[1].plot(Detunings_arr1 / (U01 / sc.h), one_atom_likelihoods_arr1, 'o-')
# axs2[1].plot(Detunings_arr2 / (U02 / sc.h), one_atom_likelihoods_arr2, 'o-')
# axs2[1].plot(Detunings_arr3 / (U03 / sc.h), one_atom_likelihoods_arr3, 'o-')
# axs2[1].plot(-red_Detunings_arr / (red_U0 / sc.h), red_one_atom_likelihoods_arr, 'o-', color = 'lightcoral')
# axs2[1].set_xlabel('Laser Detuning (Units of Trap Depth)', fontsize=20)
axs2[1].set_ylabel('Likelihood of Getting One Atom in Trap (%)', fontsize=16)
axs2[1].set_title('Single Atom Occupation Rate', fontsize=20, fontweight='bold')
axs2[1].grid(True)

# Plot 3: Two Atom Occupation Rate
# axs2[2].plot(Detunings_arrj / (U0j / sc.h), two_atom_likelihoods_arrj, 'o-')
axs2[2].plot(Detunings_arr1 / (U01 / sc.h), two_atom_likelihoods_arr1, 'o-')
# axs2[2].plot(Detunings_arr2 / (U02 / sc.h), two_atom_likelihoods_arr2, 'o-')
# axs2[2].plot(Detunings_arr3 / (U03 / sc.h), two_atom_likelihoods_arr3, 'o-')
# axs2[2].plot(-red_Detunings_arr / (red_U0 / sc.h), red_two_atom_likelihoods_arr, 'o-', color = 'lightcoral')
# axs2[2].set_xlabel('Laser Detuning (Units of Trap Depth)', fontsize=20)
# axs2[2].set_ylabel('Likelihood of Getting Two Atoms in Trap (%)', fontsize=16)
axs2[2].set_title('Two Atom Occupation Rate', fontsize=20, fontweight='bold')
axs2[2].grid(True)

# # Plot 4: Average Number of Collisions
# axs2[3].plot(Detunings_arrj / (U0j / sc.h), average_number_of_collisions_arrj, 'o-')
axs2[3].plot(Detunings_arr1 / (U01 / sc.h), average_number_of_collisions_arr1, 'o-')
# axs2[3].plot(Detunings_arr2 / (U02 / sc.h), average_number_of_collisions_arr2, 'o-')
# # axs2[3].plot(Detunings_arr3 / (U03 / sc.h), average_number_of_collisions_arr3, 'o-')
# axs2[3].plot(-red_Detunings_arr / (red_U0 / sc.h), red_average_number_of_collisions_arr, 'o-', color = 'lightcoral')
axs2[3].set_xlabel('Laser Detuning (Units of Trap Depth)', fontsize=16)
# axs2[3].set_ylabel('Average Number of Collisions', fontsize=16)
axs2[3].set_title('Average Number of Collisions', fontsize=16, fontweight='bold')
# axs2[3].grid(True)

# Adjust layout for the second figure
plt.tight_layout()

# Show the second plot
plt.show()