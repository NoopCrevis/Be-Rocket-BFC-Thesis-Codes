# -*- coding: utf-8 -*-
"""
Created on Tue Jun  3 17:00:08 2025

@author: niels
"""

# Importing needed libs
import CoolProp.CoolProp as CP  # Importing CoolProp
import matplotlib.pyplot as plt  # Importing matplotlib
import numpy as np  # Importing numpy
import pandas as pd #Importing pandas
import math


# --- Tank emptying model --- #

# Input parameters
timestep = 0.001  # [s], Simulation time step
min_temperature_n2o = 182.23  # [K], Restriction from CoolProp Package, beacause it becomes solid at that temperature!
tank_temperature = 270.15  # [K], Tank temperature, typically close to the atmospheric temperature (depends on time between tank fill and tank emptying)
combustion_chamber_pressure = 20e5  # [Pa], Pressure of the combustion chamber
pressure_drop_inj = 20e5  # [Pa], Pressure drop over the injector (Needed for the performance of the injector)
vapourizing_factor_n2o = 1.0575 # [/], Factor that multiplies the N2O amount, so this extra part will vapourize to keep the pressure
Burn_time = 6.535386021327031 # [s], Duration of the burn
mass_flow_N2O = 0.465 # [kg/s], We go for a constant mass flow rate, because the valve will limit this to the wanted flow
Transient_flow_factor = 1.05 # [/], Factor so some more mass is taken into the tank, this is mainly for transient phases, NOT for vapourizing

target_pressure = 80e5  # [Pa], Supercharged pressure, can also be seen as the maximum pressure of the weakest part
molar_mass_n2 = 0.0280134  # [kg/mol], Molar masss of nitrogen
R = 8.314  # [J/(mol*K)], Universal gas constant 

fluid = 'N2O'

# --- Pipe 1 parameter ---
D_inner_pipe_1 = 0.007 #[m], Inner diameter of pipe 1 between the tank and the valve
A_pipe_1 = ((np.pi * D_inner_pipe_1 ** 2) / 4)

# --- Valve Parameters ---
Cd = 0.8                                #[/], Estimate, should be looked into deeper in literature and experimental tests!
d_valve = 0.010                         #[m], Diameter of the valve
A_valve_max = (np.pi * d_valve ** 2) / 4      #[m²], max orifice area                          
r_ball = d_valve / 2                    #[m], radius of the ball valve
r_stem = r_ball / 2                     #[m], Stem of valve
T_break = 0.5                           #[Nm], Torque needed for the moving the valve without pressure
k_f = 0.10                             #[/], Friction coefficient of the valve seat 
W = 2.5                                   #[N], sealing force

# --- Experimental valve parameters ---
Valve_no_P_Torque = 0.06 #[Nm], torque needed to move the valve when no delta P is applied

# --- Pipe 2 parameter ---
D_inner_pipe_2 = 0.007 #[m], Inner diameter of pipe 1 between the tank and the valve
A_pipe_2 = ((np.pi * D_inner_pipe_1 ** 2) / 4)


#Calculating the initial liquid mass
initial_liquid_mass = Burn_time * mass_flow_N2O * Transient_flow_factor # [kg], Liquid mass needed for the whole burn


# Calculate initial liquid N2O properties
liquid_density_n2o = CP.PropsSI('D', 'T', tank_temperature, 'Q', 0, 'N2O') #Calculate the density of liq n2o at given temp and vapor quality 0 (=Saturated liq)

# Determine volume for N2O
liquid_volume_n2o = initial_liquid_mass / liquid_density_n2o
liquid_volume_n2o_vap = liquid_volume_n2o * vapourizing_factor_n2o

# --- Supercharging gas ---
#We first look at how much N2 (nitrogen) is needed to have min pressure before the injector
Pressure_min_in_inj = combustion_chamber_pressure + pressure_drop_inj
#Via the ideal gas law we can calcualte the mass of N2 needed for this
mol_min_in_inj = (Pressure_min_in_inj * liquid_volume_n2o_vap) / (R * tank_temperature)
#Mass of N2
mass_min_in_inj_N2 = mol_min_in_inj * molar_mass_n2
#Now we look at the space needed to store this at the max tank operating pressure, we work back via the ideal gas law
V_tank_N2_min_in_inj = (Pressure_min_in_inj * liquid_volume_n2o_vap) / (target_pressure)
#Now we have the minimum space needed to store the mass of oxidizer needed, and the supercharging gas
V_tank_N2O_N2 = liquid_volume_n2o_vap + V_tank_N2_min_in_inj
#But now that the volume is bigger the minimum pressure when the tank is empty will be lower! But this can be countered by the N2O that will evaporate (That is the 10% extra --> vapourizing_factor_n2o)


# Simulation variables
liquid_mass = initial_liquid_mass  # [kg], liquid mass
vapour_mass_N2 = mass_min_in_inj_N2  # [kg], vapour mass, there is no N2O vapour, because it is above the vapour pressure until a certain point
time = 0  # [s], starting time of model
total_vaporized_mass = 0  # Initialize accumulator

# Arrays for storing results
time_array = []

#Tank store arrays
Tank_v_out = [] #Empty array for storing the speed of the flow that comes out of the tank
tank_pressure_array = []
n2_pressure_new_array = []
n2O_volume_new_array = [] 
n2_volume_new_array = []
n2o_pressure_vapour_array = []
latent_heat_array = []
specific_heat_liquid_array = []
heat_removed_array = []
temperature_drop_array = []
tank_temperature_array = []
n2o_liq_mass_array = []
vaporized_mass_array = []
n2o_liq_density_array = []
Tank_rho_out = []

#Pipe 1 store arrays
Pipe_1_P_out = [] #Empty array for storing the outlet pressure of the first pipe
Pipe_1_T_out = [] #Empty array for storing the outlet temperature of the first pipe
Pipe_1_T_drop = []  #Empty array for storing the temperature drop or increase over the first pipe
Pipe_1_v_out = [] #Empty array for storing the outlet velocity in the first pipe
Pipe_1_v_drop = [] #Empty array for storing the velocity drop over the first pipe
Pipe_1_rho_out = []

#Valve store arrays
Valve_P_out = [] #Empty array for storing the valve outlet pressure
Valve_P_drop = [] #Empty array for storing the pressure drop over the valve
Valve_T_out = [] #Empty array for storing the valve outlet temperature
Valve_T_drop = [] #Empty array for storing the temperature drop over the valve
Valve_v_out = [] #Empty array for storing the valve outlet speed
Valve_Re_out = []
Valve_v_drop = [] #Empty array for storing the velocity drop over the valve
Valve_h_in = []
Valve_h_out = []
Valve_Q_dot = []
Valve_rho_out = []
Valve_Area = []
Valve_torques_theo = []
Valve_torques_exp = []
Valve_torques_hyb = []

#Pipe 2 store arrays
Pipe_2_P_out = [] #Empty array for storing the outlet pressure of the second pipe
Pipe_2_T_out = [] #Empty array for storing the outlet pressure of the second pipe
Pipe_2_T_drop = [] #Empty array for storing the temperature drop over the second pipe
Pipe_2_v_out = [] #Empty array for storing the outlet velocity of the second pipe
Pipe_2_v_drop = [] #Empty array for storing the velocity drop over the second pipe



#While loop for model N2O emptyting
while liquid_mass > 0:
    # Update masses
    delta_mass = mass_flow_N2O * timestep
    liquid_mass -= delta_mass
    if liquid_mass < 0:
        break
    
    #New N2O volume due to the N2O leaving the tank
    Volume_N2O_new = liquid_mass / liquid_density_n2o 
    Volume_N2_new = V_tank_N2O_N2 - Volume_N2O_new
    Pressure_N2_new = (mol_min_in_inj * R * tank_temperature) / (Volume_N2_new)
    
    #Vapour pressure of N2O
    Vapour_pressure_N2O = CP.PropsSI('P', 'T', tank_temperature, 'Q', 1, 'N2O')
    
    latent_heat = CP.PropsSI('H', 'T', tank_temperature, 'Q', 1, 'N2O') - CP.PropsSI('H', 'T', tank_temperature, 'Q', 0, 'N2O')
    specific_heat_liquid = CP.PropsSI('C', 'T', tank_temperature, 'Q', 0, 'N2O')
    
    liquid_density_new = CP.PropsSI('D', 'T', tank_temperature, 'Q', 0, 'N2O') #Calculate the density of liq n2o at given temp and vapor quality 0 (=Saturated liq)


    if Pressure_N2_new > Vapour_pressure_N2O:
        Pressure_tank = Pressure_N2_new
    else:
        # Heat loss from vaporizing nitrous
        heat_removed = delta_mass * latent_heat
        temperature_drop = heat_removed / (liquid_mass * specific_heat_liquid) if liquid_mass > 0 else 0
        tank_temperature = max(tank_temperature - temperature_drop, min_temperature_n2o)
        Pressure_tank = Vapour_pressure_N2O
        
        heat_removed_array.append(heat_removed)
        temperature_drop_array.append(temperature_drop)
        
        mass_vaporized = heat_removed / latent_heat  # This simplifies to delta_mass, but keep it explicit
        total_vaporized_mass += mass_vaporized
        
        vaporized_mass_array.append(mass_vaporized)

    V_tank_out = mass_flow_N2O / (liquid_density_new * A_pipe_1)

    #Pipe 1
    #The pressure will initially drop due to the friction, but once we are in the full N2O phase, it will vapourize a tiny bit of liquid, and the pressure will be restored
    #We could say the temperature will also drop, due to heat needed for this vapourization but this will be a tiny amount,  that the tube walls will counter this
    P_pipe_1_out = Pressure_tank
    #Density of air at outlet of the pipe 1
    Rho_pipe_1_out = CP.PropsSI('D', 'T', tank_temperature, 'Q', 0, 'N2O')
    #Velocity at outlet of pipe 1
    V_pipe_1_out = mass_flow_N2O / (Rho_pipe_1_out * A_pipe_1)
    #Velocity drop over pipe 1
    V_drop_pipe_1 = V_tank_out - V_pipe_1_out

    #Valve
    #Pressure drop & valve area (estimation)
    if P_pipe_1_out >= Pressure_min_in_inj:
        delta_P_valve = P_pipe_1_out - 40e5
        P_valve_out = 40e5
    else:
        delta_P_valve = 0
        P_valve_out = P_pipe_1_out

    
    if P_pipe_1_out >= Pressure_min_in_inj + 1e5: #This is done, otherwise a peak is created beacuse the delta_P_valve gets to small!
        A_valve = (mass_flow_N2O / (Cd * np.sqrt(2 * Rho_pipe_1_out * delta_P_valve)))
    else:
        A_valve = A_valve_max
    
    #Torque calculation valve
    F = delta_P_valve * A_valve #[N], Force needed to open valve with certain delta P and area exposed to this
    torque_valve = (F * d_valve) / 2 #[Nm], Torque needed only due to pressure on area
    torque_friction = k_f * W * r_stem #[Nm], Torque needed due to friction
    torque_total = torque_valve + torque_friction #[Nm], Total torque needed for operating the valve (not from zero movement)

    #Hybrid torque, here the pressure difference is being hold, but the friction part is changed by experimental value for
    torque_hybrid = torque_valve + Valve_no_P_Torque

    #Valve toraue via experimental formula
    torques_experimental = 0.0016 * (P_pipe_1_out / 1e5) + Valve_no_P_Torque
    

    # Get inlet enthalpy
    h_valve_in = CP.PropsSI('H', 'P', P_pipe_1_out, 'T', tank_temperature, fluid)
    # Find outlet temperature assuming isenthalpic expansion
    T_valve_out = CP.PropsSI('T', 'P', P_valve_out, 'H', h_valve_in, fluid)
    # Compute energy per second to reheat from outlet back to inlet temp
    h_valve_out = CP.PropsSI('H', 'P', P_valve_out, 'T', T_valve_out, fluid)
    Q_dot_valve = mass_flow_N2O * (h_valve_in - h_valve_out)
    T_valve_drop = tank_temperature - T_valve_out
    Rho_valve_out = CP.PropsSI('D', 'T', T_valve_out, 'Q', 0, 'N2O')
    V_valve_out = mass_flow_N2O / (Rho_valve_out * A_valve)
    
    
    #Pipe 2
    #Same reasoning as P_pipe_1_out
    #Calculate the pressure at the end of pipe 2
    P_pipe_2_out = P_valve_out
    Rho_pipe_2_out = Rho_valve_out
    #Velocity at outlet of pipe 2
    V_pipe_2_out = mass_flow_N2O / (Rho_pipe_2_out * A_pipe_2)
    #Velocity drop over pipe 1
    V_drop_pipe_2 = V_pipe_2_out - V_pipe_2_out


    # Update time
    time += timestep

    # Store of calculated values

    #torques.append(torque_total)

    #Tank store arrays
    time_array.append(time)
    n2_pressure_new_array.append(Pressure_N2_new)
    tank_pressure_array.append(Pressure_tank)
    n2O_volume_new_array.append(Volume_N2O_new)
    n2_volume_new_array.append(Volume_N2_new)
    n2o_pressure_vapour_array.append(Vapour_pressure_N2O)
    latent_heat_array.append(latent_heat)
    specific_heat_liquid_array.append(specific_heat_liquid)
    n2o_liq_density_array.append(liquid_density_new)
    tank_temperature_array.append(tank_temperature)
    n2o_liq_mass_array.append(liquid_mass)
    Tank_v_out.append(V_tank_out)

    #Pipe 1 store arrays
    
    Pipe_1_P_out.append(P_pipe_1_out)
    #Pipe_1_T_out.append(T_out_pipe_1)
    #Pipe_1_T_drop.append(T_drop_pipe_1)
    Pipe_1_v_out.append(V_pipe_1_out)
    Pipe_1_v_drop.append(V_drop_pipe_1)
    Pipe_1_rho_out.append(Rho_pipe_1_out)

    #Valve store arrays
    
    Valve_P_out.append(P_valve_out)
    Valve_P_drop.append(delta_P_valve)
    Valve_T_out.append(T_valve_out)
    Valve_T_drop.append(T_valve_drop)
    Valve_v_out.append(V_valve_out)
    #Valve_v_drop.append()
    Valve_h_in.append(h_valve_in)
    Valve_h_out.append(h_valve_out)
    Valve_Q_dot.append(Q_dot_valve)
    Valve_rho_out.append(Rho_valve_out)
    Valve_Area.append(A_valve * 1e6)
    Valve_torques_theo.append(torque_total)
    Valve_torques_exp.append(torques_experimental)
    Valve_torques_hyb.append(torque_hybrid)
    
    #Pipe 2 store arrays

    Pipe_2_P_out.append(P_pipe_2_out)
    #Pipe_2_T_out.append(T_out_pipe_2)
    #Pipe_2_T_drop.append(T_drop_pipe_2)
    Pipe_2_v_out.append(V_pipe_2_out)
    Pipe_2_v_drop.append(V_drop_pipe_2)



plt.figure(figsize=(10, 6))
plt.plot(time_array, tank_pressure_array, label='Tank Pressure [Pa]')
plt.plot(time_array, Pipe_1_P_out, label='Pipe 1 Outlet Pressure [Pa]')
plt.plot(time_array, Valve_P_out, label='Valve Outlet Pressure [Pa]')
plt.plot(time_array, Pipe_2_P_out, label='Pipe 2 Outlet Pressure [Pa]')
plt.xlabel('Burn Time [s]')
plt.ylabel('Pressure [Pa]')
plt.title('Pressures vs Burn Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(time_array, Tank_v_out, label='Tank Velocity [m/s]')
plt.plot(time_array, Pipe_1_v_out, label='Pipe 1 Outlet Velocity [m/s]')
plt.plot(time_array, Valve_v_out, label='Valve Outlet Velocity [m/s]')
plt.plot(time_array, Pipe_2_v_out, label='Pipe 2 Outlet Velocity [m/s]')
plt.xlabel('Burn Time [s]')
plt.ylabel('Velocities [m/s]')
plt.title('Velocities vs Burn Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(time_array, Valve_Area, label='Valve Area [mm^2]')
plt.xlabel('Burn Time [s]')
plt.ylabel('Valve Area [mm^2]')
plt.title('Valve Area vs Burn Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(time_array, Valve_torques_theo, label='Valve Torque Theoretical [Nm]')
plt.plot(time_array, Valve_torques_exp, label='Valve Torque Experimental [Nm]')
plt.plot(time_array, Valve_torques_hyb, label='Valve Torque Hybrid [Nm]')
plt.xlabel('Burn Time [s]')
plt.ylabel('Valve torques [Nm]')
plt.title('Ball valve torques vs Burn Time')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()