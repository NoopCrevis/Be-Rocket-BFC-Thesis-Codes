# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 18:56:39 2025

@author: niels
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# Parameters
fluid = "N2O"
T_flow = 283.15  # K, Temperature of flow
rho_liq = CP.PropsSI('D', 'T', T_flow, 'Q', 0, fluid)  # Density of liquid N2O [kg/m3]
rho_gas = CP.PropsSI('D', 'T', T_flow, 'Q', 1, fluid)  # Density of gas N2O [kg/m3]

# Input parameters
mass_liq = 3.190902224912923 - 0.1825  # kg
mass_gas = 0.1885 + 0.1825  # kg
t_burn = 6.535386021327031  # s
V_liq = mass_liq / rho_liq  # Volume of liquid [m3]
V_gas = mass_gas / rho_gas  # Volume of gas [m3]
Q_liq = V_liq / (t_burn / 3600)  # m3/h
Q_gas = V_gas / (t_burn / 3600)  # m3/h

# Pressure range
P_in_values = np.linspace(25, 80, 100)  # Inlet pressure from 25 bar to 80 bar
P_out = 25  # Outlet pressure [bar]

Kv_liq_values = []
Kv_gas_values = []

for P_in in P_in_values:
    Delta_P_liq = P_in - P_out  # Pressure drop for liquid phase
    Delta_P_gas = P_in - P_out  # Pressure drop for gas phase
    
    # Liquid phase Kv calculation
    if Delta_P_liq > 0:
        Kv_liq = Q_liq * np.sqrt((1 / Delta_P_liq) * (rho_liq / 1000))
    else:
        Kv_liq = np.nan  # Avoid division by zero or negative values
    Kv_liq_values.append(Kv_liq)
    
    # Gas phase Kv calculation
    if P_out > P_in / 2:  # Subsonic flow
        Kv_gas = Q_gas / 514 * np.sqrt((rho_gas * T_flow) / (Delta_P_gas * P_out))
    else:  # Sonic flow
        Kv_gas = Q_gas / (257 * P_in) * np.sqrt(rho_gas * T_flow)
    Kv_gas_values.append(Kv_gas)

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(P_in_values, Kv_liq_values, label='Kv Liquid', linestyle='-', marker='o', markersize=3)
plt.plot(P_in_values, Kv_gas_values, label='Kv Gas', linestyle='-', marker='s', markersize=3)
plt.xlabel('Inlet Pressure (bar)')
plt.ylabel('Kv Value')
plt.title('Kv Value vs Inlet Pressure')
plt.legend()
plt.grid()
plt.show()
