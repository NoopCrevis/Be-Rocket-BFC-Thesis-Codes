# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 20:57:23 2025

@author: niels
"""

import CoolProp.CoolProp as CP

fluid = 'NitrousOxide'
P_in = 80e5     # 80 bar
T_in = 283.15   # 10°C in Kelvin, will start to have troubles close to saturation line
P_out = 40e5    # 40 bar

# Get inlet enthalpy
h_in = CP.PropsSI('H', 'P', P_in, 'T', T_in, fluid)

# Find outlet temperature assuming isenthalpic expansion
T_out = CP.PropsSI('T', 'P', P_out, 'H', h_in, fluid)

# Print results
print(f"Inlet Temperature: {T_in - 273.15:.2f} °C")
print(f"Outlet Temperature: {T_out - 273.15:.2f} °C")

mass_flow = 0.5  # kg/s

# Compute energy per second to reheat from outlet back to inlet temp
h_out = CP.PropsSI('H', 'P', P_out, 'T', T_out, fluid)
Q_dot = mass_flow * (h_in - h_out)

print(f"Heating Power Needed: {Q_dot:.10f} W")
