# -*- coding: utf-8 -*-
"""
Created on Mon May 12 09:48:06 2025
@author: niels
"""

import pandas as pd
import matplotlib.pyplot as plt

# Load the data
df = pd.read_csv("Valve_Air_Test_3_40_bar_05062025.txt")
#df = pd.read_csv("Valve_Air_Test_12_15052025.txt")

# Create a time axis based on number of samples (assuming uniform sampling)
df["Time"] = range(len(df))

# --- Calibration of Sensors ---

# Torque sensor
Torque_Cal = 1e-6 * df["Torque_raw"] - 0.007
#Torque_Cal = df["Torque_raw"] 


# Load Cell sensor
Load_Cal = 2e-5 * df["Load_raw"] - 0.0023

# Temperature sensor offsets
Temp_1_Cal = df["Temp1_C"] - 3.72
Temp_2_Cal = df["Temp2_C"] - 2.94
Temp_3_Cal = df["Temp3_C"] - 4.9
Temp_4_Cal = df["Temp4_C"] - 1.22

# Pressure sensor offsets
Press_1_Cal = df["Pressure1_bar"] - 2.39
Press_2_Cal = df["Pressure2_bar"] - 2.39
Press_3_Cal = df["Pressure3_bar"] - 2.39
Press_4_Cal = df["Pressure4_bar"] - 3.06

# --- Valve Angle Conversion ---

valve_closed_V = 1.290
valve_open_V = 2.280
valve_range_deg = 90.0

df["ValveAngle_deg"] = (df["ValvePosition_V"] - valve_closed_V) / (valve_open_V - valve_closed_V) * valve_range_deg
df["ValveAngle_deg"] = df["ValveAngle_deg"].clip(lower=0, upper=90)

# Subset
Under_Limit_Range = 1650
Upper_Limit_Range = 1700
df_subset = df.iloc[Under_Limit_Range:Upper_Limit_Range]

# --- Quantitive values ---

# --- Quantitive values ---

# 1. Maximum and Minimum Torque
#max_torque = Torque_Cal.max()
min_torque = Torque_Cal.min()
#print(f"Max Torque: {max_torque:.4f} Nm")
print(f"Min Torque: {min_torque:.4f} Nm")

# 2. Mass loss estimation
mass_full = Load_Cal.iloc[0]
print(f"Esstimated mass gain: {mass_full:.4f} kg")
#mass_empty = Load_Cal.iloc[-1]
#mass_loss = mass_full - mass_empty
#print(f"Estimated Mass Loss: {mass_loss:.4f} kg (start: {mass_full:.4f} kg, end: {mass_empty:.4f} kg)")


# 4. Temperature drop from highest to last in Temp4
temp4_max = Temp_4_Cal.max()
temp4_last = Temp_4_Cal.iloc[-1]
temp_drop = temp4_max - temp4_last
print(f"Temp4 Drop: {temp_drop:.2f} °C (max: {temp4_max:.2f} °C, last: {temp4_last:.2f} °C)")


# --- Full Time Plots ---

# Plot Torque sensor on the left Y-axis and valve angle on the right Y-axis
plt.figure(figsize=(12, 6))

# Left Y-axis: Torque sensor
ax1 = plt.gca()
ax1.plot(df["Time"], Torque_Cal, label="Torque (Nm)", color='red')
ax1.set_xlabel("Time (sample number)")
ax1.set_ylabel("Torque (Nm)", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True)
ax1.legend(loc="upper left")

# Right Y-axis: Valve Angle
ax2 = ax1.twinx()
ax2.plot(df["Time"], df["ValveAngle_deg"], label="Valve Angle (°)", color='purple', linewidth=2, linestyle='--')
ax2.set_ylabel("Valve Angle (°)", color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# Title and layout
plt.title("Torque and Valve Angle vs Time")
plt.tight_layout()
plt.show()


# Plot Load cell on the left Y-axis and valve angle on the right Y-axis
plt.figure(figsize=(12, 6))

# Left Y-axis: Load Cell
ax1 = plt.gca()
ax1.plot(df["Time"], Load_Cal, label="Mass (kg)", color='red')
ax1.set_xlabel("Time (sample number)")
ax1.set_ylabel("Mass (kg)", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True)
ax1.legend(loc="upper left")

# Right Y-axis: Valve Angle
ax2 = ax1.twinx()
ax2.plot(df["Time"], df["ValveAngle_deg"], label="Valve Angle (°)", color='purple', linewidth=2, linestyle='--')
ax2.set_ylabel("Valve Angle (°)", color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# Title and layout
plt.title("Mass and Valve Angle vs Time")
plt.tight_layout()
plt.show()


# Plot all 4 Pressures on the left Y-axis and valve angle on the right Y-axis
plt.figure(figsize=(12, 6))

# Left Y-axis: Pressures
ax1 = plt.gca()
ax1.plot(df["Time"], Press_1_Cal, label="Pressure1 (Bar)", color='red')
ax1.plot(df["Time"], Press_2_Cal, label="Pressure2 (Bar)", color='orange')
ax1.plot(df["Time"], Press_3_Cal, label="Pressure3 (Bar)", color='green')
ax1.plot(df["Time"], Press_4_Cal, label="Pressure4 (Bar)", color='blue')
ax1.set_xlabel("Time (sample number)")
ax1.set_ylabel("Pressure (Bar)", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True)
ax1.legend(loc="upper left")

# Right Y-axis: Valve Angle
ax2 = ax1.twinx()
ax2.plot(df["Time"], df["ValveAngle_deg"], label="Valve Angle (°)", color='purple', linewidth=2, linestyle='--')
ax2.set_ylabel("Valve Angle (°)", color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# Title and layout
plt.title("Pressure and Valve Angle vs Time")
plt.tight_layout()
plt.show()


# Plot all 4 temperatures on the left Y-axis and valve angle on the right Y-axis
plt.figure(figsize=(12, 6))

# Left Y-axis: Temperatures
ax1 = plt.gca()
ax1.plot(df["Time"], Temp_1_Cal, label="Temp1 (°C)", color='red')
ax1.plot(df["Time"], Temp_2_Cal, label="Temp2 (°C)", color='orange')
ax1.plot(df["Time"], Temp_3_Cal, label="Temp3 (°C)", color='green')
ax1.plot(df["Time"], Temp_4_Cal, label="Temp4 (°C)", color='blue')
ax1.set_xlabel("Time (sample number)")
ax1.set_ylabel("Temperature (°C)", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True)
ax1.legend(loc="upper left")

# Right Y-axis: Valve Angle
ax2 = ax1.twinx()
ax2.plot(df["Time"], df["ValveAngle_deg"], label="Valve Angle (°)", color='purple', linewidth=2, linestyle='--')
ax2.set_ylabel("Valve Angle (°)", color='purple')
ax2.tick_params(axis='y', labelcolor='purple')

# Title and layout
plt.title("Temperatures and Valve Angle vs Time")
plt.tight_layout()
plt.show()



# --- Subset Plots ---

# Torque vs Valve Angle (subset)
plt.figure(figsize=(12, 6))
ax1 = plt.gca()
ax1.plot(df_subset["Time"], Torque_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Torque (Nm)", color='red')
ax1.set_xlabel("Time (sample number)")
ax1.set_ylabel("Torque (Nm)", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True)
ax1.legend(loc="upper left")
ax2 = ax1.twinx()
ax2.plot(df_subset["Time"], df_subset["ValveAngle_deg"], label="Valve Angle (°)", color='purple', linewidth=2, linestyle='--')
ax2.set_ylabel("Valve Angle (°)", color='purple')
ax2.tick_params(axis='y', labelcolor='purple')
plt.title("Torque and Valve Angle vs Time (Subset)")
plt.tight_layout()
plt.show()

# Load vs Valve Angle (subset)
plt.figure(figsize=(12, 6))
ax1 = plt.gca()
ax1.plot(df_subset["Time"], Load_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Mass (kg)", color='red')
ax1.set_xlabel("Time (sample number)")
ax1.set_ylabel("Mass (kg)", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True)
ax1.legend(loc="upper left")
ax2 = ax1.twinx()
ax2.plot(df_subset["Time"], df_subset["ValveAngle_deg"], label="Valve Angle (°)", color='purple', linewidth=2, linestyle='--')
ax2.set_ylabel("Valve Angle (°)", color='purple')
ax2.tick_params(axis='y', labelcolor='purple')
plt.title("Mass and Valve Angle vs Time (Subset)")
plt.tight_layout()
plt.show()

# Pressures vs Valve Angle (subset)
plt.figure(figsize=(12, 6))
ax1 = plt.gca()
ax1.plot(df_subset["Time"], Press_1_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Pressure1 (Bar)", color='red')
ax1.plot(df_subset["Time"], Press_2_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Pressure2 (Bar)", color='orange')
ax1.plot(df_subset["Time"], Press_3_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Pressure3 (Bar)", color='green')
ax1.plot(df_subset["Time"], Press_4_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Pressure4 (Bar)", color='blue')
ax1.set_xlabel("Time (sample number)")
ax1.set_ylabel("Pressure (Bar)", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True)
ax1.legend(loc="upper left")
ax2 = ax1.twinx()
ax2.plot(df_subset["Time"], df_subset["ValveAngle_deg"], label="Valve Angle (°)", color='purple', linewidth=2, linestyle='--')
ax2.set_ylabel("Valve Angle (°)", color='purple')
ax2.tick_params(axis='y', labelcolor='purple')
plt.title("Pressure and Valve Angle vs Time (Subset)")
plt.tight_layout()
plt.show()

# Temperatures vs Valve Angle (subset)
plt.figure(figsize=(12, 6))
ax1 = plt.gca()
ax1.plot(df_subset["Time"], Temp_1_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Temp1 (°C)", color='red')
ax1.plot(df_subset["Time"], Temp_2_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Temp2 (°C)", color='orange')
ax1.plot(df_subset["Time"], Temp_3_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Temp3 (°C)", color='green')
ax1.plot(df_subset["Time"], Temp_4_Cal.iloc[Under_Limit_Range:Upper_Limit_Range], label="Temp4 (°C)", color='blue')
ax1.set_xlabel("Time (sample number)")
ax1.set_ylabel("Temperature (°C)", color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.grid(True)
ax1.legend(loc="upper left")
ax2 = ax1.twinx()
ax2.plot(df_subset["Time"], df_subset["ValveAngle_deg"], label="Valve Angle (°)", color='purple', linewidth=2, linestyle='--')
ax2.set_ylabel("Valve Angle (°)", color='purple')
ax2.tick_params(axis='y', labelcolor='purple')
plt.title("Temperatures and Valve Angle vs Time (Subset)")
plt.tight_layout()
plt.show()


