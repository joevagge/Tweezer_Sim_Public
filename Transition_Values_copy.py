#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:19:11 2024

@author: joe
"""

# ### D2 Line ######

# D2_P_Base_to_S_Base = 384.230484468562e12
# D2_P_Base_to_F1 = 229.8518e6
# D2_P_Base_to_F3 = 193.7407e6
# D2_S_Base_to_F2 = 2.563005979089109e9
# D2_S_Base_to_F1 = 4.271676631815181e9

# AOM = (107.2 * 2) * 1e6
# Detuning = 0

# Blue_Repump_D2 = D2_P_Base_to_S_Base - D2_P_Base_to_F1 - D2_S_Base_to_F2 - AOM + Detuning
# BLue_Cooling_D2 = D2_P_Base_to_S_Base - D2_P_Base_to_F1 + D2_S_Base_to_F1 - AOM + Detuning
# print(rf'D2 Blue Cooling= {BLue_Cooling_D2*1e-12:.6f} THz')
# print(rf'D2 Blue Repump = {Blue_Repump_D2*1e-12:.6f} THz')

# print(rf'2fdelta = {(BLue_Cooling_D2 - Blue_Repump_D2)*1e-9:.12f} GHz')


# Red_Cooling_D2 = D2_P_Base_to_S_Base + D2_P_Base_to_F3 - D2_S_Base_to_F2 - AOM 
# print(rf'D2 Red Cooling = {Red_Cooling_D2*1e-12:.6f} THz')

### D1 Line ######

AOM = (107.2 * 2) * 1e6
Detuning = 0

D1_P_Base_to_S_Base = 377.107463380e12
D1_P_Base_to_F1 = 509.06e6
D1_P_Base_to_F2 = 305.44e6
D1_S_Base_to_F2 = 2.563005979089109e9
D1_S_Base_to_F1 = 4.271676631815181e9

# Cooling = F=1 -> F'=1
# Repump = F=2 -> F'=1

# Blue_Repump_D1 = D1_P_Base_to_S_Base - D1_S_Base_to_F2 - D1_P_Base_to_F1 - AOM + Detuning
# BLue_Cooling_D1 = D1_P_Base_to_S_Base + D1_S_Base_to_F1 - D1_P_Base_to_F1 - AOM + Detuning
# print(rf'D1 Blue Cooling= {BLue_Cooling_D1*1e-12:.6f} THz')
# print(rf'D1 Blue Repump = {Blue_Repump_D1*1e-12:.6f} THz')

# AP = (BLue_Cooling_D1 - Blue_Repump_D1)
# print(rf'AP = {AP*1e-12:.6f} THz')


# Cooling = F=2 -> F'=2
# Repump = F=1 -> F'=2

BLue_Cooling_D1 = D1_P_Base_to_S_Base - D1_S_Base_to_F2 + D1_P_Base_to_F2 - AOM + Detuning
Blue_Repump_D1 = D1_P_Base_to_S_Base + D1_S_Base_to_F1 + D1_P_Base_to_F2 - AOM + Detuning
print(rf'D1 Blue Cooling= {BLue_Cooling_D1*1e-12:.6f} THz')
print(rf'D1 Blue Repump = {Blue_Repump_D1*1e-12:.6f} THz')