# -*- coding: utf-8 -*-
"""
@author: niels 
"""

import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP  # Importing CoolProp

# --- Tank Parameters ---
T_atm = 283.15 # [K]
P_tank = 80 # [Bar]

# --- Fluid and pipe parameters ---
rho = CP.PropsSI('D', 'T', T_atm, 'Q', 0, 'N2O')
c =  CP.PropsSI('A', 'T', T_atm, 'P', P_tank, 'N2O') # [m/s], speed of sound in liquid N2O
diameter = 0.008    # [m]
area = np.pi * (diameter / 2) ** 2  # [m²]

# --- Flow parameters ---
mass_flow_N2O = 0.465 # [kg/s]
vol_flow = mass_flow_N2O / rho  # [m³/s]
v0 = vol_flow / area  # [m/s], initial flow velocity 

# --- Simulation settings ---
t_max = 3  # [s] total sim time
time_steps = 1000
time = np.linspace(0, t_max, time_steps)

# --- Closure time sweep ---
closure_times = np.linspace(0.25, 1.5, 1000)  # [s]
profile_type = "s_curve"  # Options: "linear", "sine", "s_curve"
pressure_limit_bar = 80  # surge pressure limit in bar
temperature_limit_C = 36  # thermal decomposition threshold for N2O

# --- Velocity profile function ---
def velocity_profile(t, t_close, profile_type):
    x = t / t_close
    if profile_type == "linear":
        return np.clip(1 - x, 0, 1)
    elif profile_type == "sine":
        return np.where(t < t_close, np.cos(x * np.pi / 2), 0)
    elif profile_type == "s_curve":
        return np.where(t < t_close, 1 - 3 * x**2 + 2 * x**3, 0)
    else:
        raise ValueError("Invalid profile type")

# --- Simulation loop ---
max_pressures = []
max_temperatures = []

# Constants
T1_C = 20.0
T1_K = T1_C + 273.15
P1 = 80e5  # Pa
gamma = 1.15  # adiabatic index for liquid N2O

for t_close in closure_times:
    velocity = velocity_profile(time, t_close, profile_type) * v0
    dV_dt = np.gradient(velocity, time)
    delta_P = np.where(time <= t_close, rho * c * np.abs(dV_dt), 0)

    # Surge pressure
    max_delta_P_bar = np.max(delta_P) / 1e5
    max_pressures.append(max_delta_P_bar)

    # Temperature rise
    P2 = P1 + delta_P
    T2_K = T1_K * (P2 / P1)**((gamma - 1)/gamma)
    T2_C = T2_K - 273.15
    max_temp = T2_C.max()
    max_temperatures.append(max_temp)

# --- Combined Safety Condition: Pressure < 80 bar AND Temp < 36 °C ---
safe_mask = (np.array(max_pressures) < pressure_limit_bar) & (np.array(max_temperatures) < temperature_limit_C)
safe_idx = np.argmax(safe_mask)
safe_time_ms = closure_times[safe_idx] * 1000 if safe_mask.any() else None

# --- Plotting ---
fig, ax1 = plt.subplots(figsize=(10, 6))

color1 = 'tab:blue'
ax1.set_xlabel('Valve Closure Time [ms]')
ax1.set_ylabel('Max Pressure Surge Delta P [bar]', color=color1)
ax1.plot(closure_times * 1000, max_pressures, label='Max Pressure Surge', color=color1)
ax1.axhline(y=pressure_limit_bar, color='red', linestyle='--', label=f'Limit: {pressure_limit_bar} bar')
ax1.tick_params(axis='y', labelcolor=color1)

if safe_time_ms:
    ax1.axvline(x=safe_time_ms, color='green', linestyle='--', label=f'Safe closure ≥ {safe_time_ms:.1f} ms')
    print(f"Safe valve closure time (for < {pressure_limit_bar} bar and < {temperature_limit_C}°C): {safe_time_ms:.1f} ms")
else:
    print("No safe closure time found below pressure and temperature limits.")

# Second Y-axis for temperature
ax2 = ax1.twinx()
color2 = 'tab:orange'
ax2.set_ylabel('Max Surge Temperature [°C]', color=color2)
ax2.plot(closure_times * 1000, max_temperatures, label='Max Surge Temperature', color=color2)
ax2.axhline(y=temperature_limit_C, color='darkorange', linestyle='--', label=f'Limit: {temperature_limit_C} °C')
ax2.tick_params(axis='y', labelcolor=color2)

# Combined legend
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.title(f'Surge Pressure and Temperature vs. Valve Closure Time ({profile_type.capitalize()} Profile)')
plt.grid(True)
plt.tight_layout()
plt.show()
