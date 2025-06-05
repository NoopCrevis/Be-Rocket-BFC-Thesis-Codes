# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 20:18:19 2025

@author: niels
"""


import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
import math

# --- System Parameters ---
fluid = 'Air'                           #The fluid is defined, air or water at the moment
V_tank = 0.006                          #[m³], Volume of the tank
P0 = 40e5                               #[Pa], Tank pressure (Starting temp)
T0 = 273.15 + 16                                #[K], Tank temperature (Starting temp)
P_atm = 101325                             #[Pa], Outlet pressure (Need to look into outlet like injector (40 bar))

# --- Tank Wall Properties ---
T_wall = T0                      # Initial tank wall temperature [K]
h = 100.0                        # Heat transfer coefficient [W/m²·K] (tunable)
A_wall = 0.53                    # Estimated inner tank surface area [m²]
m_wall = 16.5                    # Tank wall mass [kg]
c_wall = 900.0                   # Specific heat capacity of wall [J/kg·K]
adiabatic = False                # Set to True for adiabatic behavior, False for heat transfer from walls (More realistic)

# --- Pipe 1 parameter ---
L_pipe_1 = 0.150 # [m], Pipe 1 length between the tank and the valve
D_inner_pipe_1 = 0.0085 #[m], Inner diameter of pipe 1 between the tank and the valve
epsilon_pipe_1 = 0.00020  # [m], pipe roughness
A_pipe_1 = ((np.pi * D_inner_pipe_1 ** 2) / 4)
T_wall_pipe_1 = T0

# --- Valve Parameters ---
Cd = 0.8                                #[/], Estimate, should be looked into deeper in literature and experimental tests!
d_valve = 0.010                         #[m], Diameter of the valve
A_max = (np.pi * d_valve ** 2) / 4      #[m²], max orifice area 
t_open = 2.0                            #[s], time to fully open valve 
r_ball = d_valve / 2                    #[m], radius of the ball valve
r_stem = r_ball / 2                     #[m], Stem of valve
T_break = 0.06                          #[Nm], Torque needed for the moving the valve without pressure
k_f = 2e-6                              #[/], Friction coefficient of the valve seat
W = 2                                   #[N], sealing force

# --- Experimental valve parameters ---
Valve_no_P_Torque = 0.06 #[Nm], torque needed to move the valve when no delta P is applied

# --- Pipe 2 parameter ---
L_pipe_2 = 0.130 # [m], Pipe 1 length between the tank and the valve
D_inner_pipe_2 = 0.0085 #[m], Inner diameter of pipe 1 between the tank and the valve
epsilon_pipe_2 = 0.00020  # [m], pipe roughness
A_pipe_2 = ((np.pi * D_inner_pipe_1 ** 2) / 4)
T_wall_pipe_2 = T0

# --- Gas Properties ---
R = CP.PropsSI('GAS_CONSTANT', fluid) / CP.PropsSI('MOLAR_MASS', fluid)                                 #[J/(k mol)], Gas constant used for further calculations
gamma = CP.PropsSI('CPMASS', 'T', T0, 'P', P0, fluid) / CP.PropsSI('CVMASS', 'T', T0, 'P', P0, fluid)   #[/], Gamma, depends on T, and P, so it can change with during the expansion of the gas
rho0 = CP.PropsSI('D', 'T', T0, 'P', P0, fluid)                                                         #[kg / m3], Density at time 0 (Starting density), will change with the T and P, so flow rate will be affected
m0 = rho0 * V_tank                                                                                      #[kg], Mass at time 0 (Starting mass), depends on density which depends on the pressure and temperature
mu = 1.9635*1e-5  # Pa.s (dynamic viscosity of oxidizer)

# --- Time Setup ---
dt = 0.01                   #[s], Time step, lower --> Finer graph, higher --> Coarser graph
t_max = 10                  #[s], Maximum simulation time (Make sure it's long enough, so the whole flow is simaluted, it will stop earlier if mass is low enough)
steps = int(t_max / dt)     #[/], Amount of steps for the simulation

# --- Arrays to store results ---
#General store arrays
areas, torques_theo, torques_hyb, torques_exp, T_outlet, angle = [], [], [], [], [], [],
times = [] #Empty array for storing the time steps of the simulation
mdots = [] #Empty array for storing the mass flow during the time steps

#Tank store arrays
Tank_pressure = [] #Empty array for storing the tank pressure
Tank_temperature = [] #Empty array for storing the tank temperature
Tank_mass = [] #Empty array for storing the tank mass (only fluid inside)
Tank_v_out = [] #Empty array for storing the speed of the flow that comes out of the tank
Tank_Re_out = []
Tank_M_out = []
Tank_a_out = []
Tank_rho = []

#Pipe 1 store arrays
Pipe_1_P_out = [] #Empty array for storing the outlet pressure of the first pipe
Pipe_1_P_drop = [] #Empty array for storing the pressure drop over the first pipe
Pipe_1_T_out = [] #Empty array for storing the outlet temperature of the first pipe
Pipe_1_T_drop = []  #Empty array for storing the temperature drop or increase over the first pipe
Pipe_1_v_out = [] #Empty array for storing the outlet velocity in the first pipe
Pipe_1_v_drop = [] #Empty array for storing the velocity drop over the first pipe
Pipe_1_Re_out = []
Pipe_1_rho_out = []

#Valve store arrays
Valve_P_out = [] #Empty array for storing the valve outlet pressure
Valve_P_drop = [] #Empty array for storing the pressure drop over the valve
Valve_T_out = [] #Empty array for storing the valve outlet temperature
Valve_T_drop = [] #Empty array for storing the temperature drop over the valve
Valve_v_out = [] #Empty array for storing the valve outlet speed
Valve_Re_out = []
Valve_v_drop = [] #Empty array for storing the velocity drop over the valve

#Pipe 2 store arrays
Pipe_2_P_out = [] #Empty array for storing the outlet pressure of the second pipe
Pipe_2_P_drop = [] #Empty array for storing the pressure drop over the second pipe
Pipe_2_T_out = [] #Empty array for storing the outlet pressure of the second pipe
Pipe_2_T_drop = [] #Empty array for storing the temperature drop over the second pipe
Pipe_2_v_out = [] #Empty array for storing the outlet velocity of the second pipe
Pipe_2_Re_out = []
Pipe_2_v_drop = [] #Empty array for storing the velocity drop over the second pipe

#As you can see there are no "in" or "inlet" values in these store arrays above, this is becasue the outlet form the previous is equal to the inlet of the following part
#The starting conditions for now are stated via the initial conditions (here below, and in the parameters at the top of the script)

# --- Initial Conditions ---
P, T, m = P0, T0, m0 #This lines sets the inital conditions (P = P0, T = T0 and m = m0)

def valve_area(t): #Definition for the valve area dependend only on the time t
    if t < t_open: #As long as t is smaller than t_opne the valve will keep going more open
        angle = (np.pi / 2) * (t / t_open) #Sine function is used to simulate a ball valve (circle)
        return A_max * np.sin(angle) #Area of the valve at certain t
    else: #If t = t_open --> The valve movement will stop, and A_max will be used for the rest of the simulation time
        return A_max #Usage of A_max for area

def mass_flow_rate(P, T, P_down, A, Cd): #Definition for the mass flow rate with variables: Pressure, Temperature, Pressure downstream, Area valve and discharge coefficient
    critical_pressure = P * (2 / (gamma + 1)) ** (gamma / (gamma - 1)) #Calculating the critical pressure, mainly depends on the P in the tank and the gamma of the medium
    if P_down <= critical_pressure: #If loop, if downstream pressure is < critical pressure --> Chocked (Shockwaves & choked flow)
        return Cd * A * P * np.sqrt(gamma / (R * T) * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))) #Return the mas flow rate for the choked flow
    else: #Else loop for if loop, when P_down > Critical pressure --> no choked flow, no shockwaves
        pressure_ratio = P_down / P #Pressure ratio for further calculation mass flow rate
        term = 2 * gamma / (R * T * (gamma - 1)) #Term to make formula mass flow rate a bit shorter
        return Cd * A * P * np.sqrt(term * (pressure_ratio ** (2 / gamma) - pressure_ratio ** ((gamma + 1) / gamma))) #Returns mass flow rate (non choked)


# --- Time Loop ---
for step in range(steps): #For loop to calculate all the values for certain amount of steps
    time = step * dt #[s], Total time of simulation (Can be shorter or equal to t_max)
    A = valve_area(time) #[m2], Area of valve at certain time step
    mdot = mass_flow_rate(P, T, P_atm, A, Cd) #[kg/s], Mass flow rate via definition
    dm = mdot * dt #[kg], Partial mass for each time step
    m -= dm #[kg], Medium mass in tank after some mass "dm" left the system / tank

    if m <= 1e-5: #If statement
        break #Simulation will stop once the mass is smaller than 0.01 kg, otherwise the simulation will stop only at t_max, with very little change for some values




    #Tank
    #The density depends on the temperature and pressure inside of the tank, so this will change
    Rho_tank = P / (R * T)
    V_tank_out = mdot / (Rho_tank * A_pipe_1)


    if adiabatic: #If statement for adiabatic condition
        T = T0 * (m / m0)**(gamma - 1) #Classic isentropic expansion
    else: #Approximate energy balance with heat transfer
        Q_dot = h * A_wall * (T_wall - T)  # Heat from wall to gas [W]
        dU_gas = -R * T * dm + Q_dot * dt  # Internal energy change: mass loss + heating
        cv = CP.PropsSI('CVMASS', 'T', T, 'P', P, fluid)

        if m > 1e-6 and cv > 0: #Boundaries set so the calculations are still in feasible boundaries
            T += dU_gas / (m * cv) #Adding value so the effect of heat from the walls is taken into account
        else:
            break  #Avoid division by near-zero or bad physics

        #Update wall temperature from energy loss
        T_wall -= (Q_dot * dt) / (m_wall * c_wall) #The wall itself will cooldown, because heat is transfered from the wall to the medium (will maybe be measured)

    if not np.isfinite(T) or T < 1.0: #Statement if the some non real values start to appear
        print(f"Non-physical temperature detected at t={time:.3f}s: T={T:.2f} K")
        break #Simulation will stop, once these values are detected

    P = m / V_tank * R * T #[Pa], With new temperature (More realistic), pressure is calculated via ideal gas law --> Look at other formulas?

    a_tank_out = CP.PropsSI('A', 'T', T, 'P', P, fluid) #[m/s], speed of sound tank
    M_tank_out = V_tank_out / a_tank_out


    #Pipe 1
    #The mass flow through the pipe depends on the mass flow through the valve (so directly dependend on it)
    #The area depends on the diameter of tube, here the flow will be slower than in through the valve
    
    mu_pipe_1 = CP.PropsSI('VISCOSITY', 'T', T, 'P', P, fluid)  # Dynamic viscosity in [Pa·s]

    #Pressure drop through pipe 1
    def reynolds_number(Rho_tank, P_tank_out, D_inner_pipe_1, mu_pipe_1):
        """Calculate Reynolds number"""
        return (Rho_tank * V_tank_out * D_inner_pipe_1) / mu_pipe_1

    def friction_factor(Re_pipe_1, D_inner_pipe_1, epsilon_pipe_1):
        """Calculate friction factor using laminar formula or Colebrook equation for turbulent flow"""
        if Re_pipe_1 < 2300:
            return 64 / Re_pipe_1  # Laminar flow
        else:
            # Colebrook equation (implicit, solved iteratively)
            def colebrook(f):
                arg = epsilon_pipe_1 / (3.7 * D_inner_pipe_1) + 2.51 / (Re_pipe_1 * math.sqrt(f))
                if arg <= 0:
                    raise ValueError("Invalid argument for log10 in Colebrook equation")
                return 1 / math.sqrt(f) + 2.0 * math.log10(arg)
            
            # Initial guess and solving using Newton-Raphson method
            f_guess = 0.05  # Safer initial guess
            for _ in range(20):  # Iterate to refine f
                try:
                    f_new = f_guess - colebrook(f_guess) / (-0.5 / (f_guess ** 1.5))
                    if abs(f_new - f_guess) < 1e-6:
                        break
                    f_guess = f_new
                except ValueError:
                    print("Warning: Could not solve Colebrook equation, returning estimated friction factor.")
                    return 0.02  # Return a default reasonable value if the solver fails
            return f_guess

    def pressure_drop_pipe_1(Rho_tank, V_tank_out, D_inner_pipe_1, L_pipe_1, mu_pipe_1, epsilon_pipe_1):
        """Calculate pressure drop using Darcy-Weisbach equation"""
        Re_pipe_1 = reynolds_number(Rho_tank, V_tank_out, D_inner_pipe_1, mu_pipe_1)
        f = friction_factor(Re_pipe_1, D_inner_pipe_1, epsilon_pipe_1)
        delta_P_pipe_1 = (f * (L_pipe_1 / D_inner_pipe_1) * (Rho_tank * V_tank_out ** 2) / 2)
        return delta_P_pipe_1, Re_pipe_1, f

    # Compute pressure drop in pipe 1
    delta_P_pipe_1, Re_pipe_1, f = pressure_drop_pipe_1(Rho_tank, V_tank_out, D_inner_pipe_1, L_pipe_1, mu_pipe_1, epsilon_pipe_1)
    #Pressure at outlet of the pipe
    P_pipe_1_out = P - delta_P_pipe_1
    #Density of air at outlet of the pipe 1
    Rho_pipe_1_out = P_pipe_1_out / (R * T) #Assume T has not changed yet!!!
    #Velocity at outlet of pipe 1
    V_pipe_1_out = mdot / (Rho_pipe_1_out * A_pipe_1) #[m/s], Speed of flow through pipe 1
    #Velocity drop over pipe 1
    V_drop_pipe_1 = V_tank_out - V_pipe_1_out


    def update_temperature_pipe_1(Rho_pipe_1_out, V_pipe_1_out, D_inner_pipe_1, L_pipe_1, mu, T, T_wall_pipe_1, mdot, fluid):
        # Reynolds number
        Re_pipe_1 = (Rho_pipe_1_out * V_pipe_1_out * D_inner_pipe_1) / mu
    
        # Friction factor
        if Re_pipe_1 < 2300:
            f = 64 / Re_pipe_1  # Laminar flow
        else:
            # Colebrook-White equation solver (approximation)
            epsilon = 0.00007  # m (pipe roughness) <-- you could also pass this in
            relative_roughness = epsilon / D_inner_pipe_1
            f_guess = 0.02
            for _ in range(10):
                f_guess = (-2 * np.log10(relative_roughness/3.7 + 2.51/(Re_pipe_1*np.sqrt(f_guess))))**-2
            f = f_guess
    
        # Pressure drop
        delta_P_pipe_1 = f * (L_pipe_1 / D_inner_pipe_1) * (Rho_pipe_1_out * V_pipe_1_out**2) / 2
    
        # Temperature change (simplified)
        cp_pipe_1 = CP.PropsSI("CPMASS", "T", T, "P", P0, fluid)   # Specific heat capacity at constant pressure
        m_pipe_1 = Rho_pipe_1_out * (np.pi * (D_inner_pipe_1/2)**2) * L_pipe_1  # mass inside pipe section
    
        q_dot_pipe_1 = 2.5  # W/m²·K estimate --> Match with experimental setup !
        
        A_inner = np.pi * D_inner_pipe_1 * L_pipe_1  # Internal surface area
        heat_transfer = q_dot_pipe_1 * A_inner * (T_wall_pipe_1 - T)  # Heat gain or loss
    
        if m_pipe_1 > 0 and cp_pipe_1 > 0:
            dT_pipe_1 = heat_transfer / (m_pipe_1 * cp_pipe_1)
        else:
            dT_pipe_1 = 0.0
    
        T_out_pipe_1 = T + dT_pipe_1
    
        # Nusselt number
        if Re_pipe_1 < 2300:
            Nu_pipe_1 = 3.66  # Constant for fully developed laminar flow
        else:
            Nu_pipe_1 = 0.023 * (Re_pipe_1**0.8) * (0.7**0.3)  # Dittus-Boelter for turbulent flow, Pr ≈ 0.7 for air
    
        return T_out_pipe_1, heat_transfer, q_dot_pipe_1, Re_pipe_1, Nu_pipe_1, dT_pipe_1
    
    T_out_pipe_1, h_pipe_1, q_dot_pipe_1, Re_pipe_1, Nu_pipe_1, T_drop_pipe_1 = update_temperature_pipe_1(Rho_pipe_1_out, V_pipe_1_out, D_inner_pipe_1, L_pipe_1, mu, T, T_wall_pipe_1, mdot, fluid)


    #Valve
    #Pressure drop
    delta_P_valve = mass_flow_rate(P_pipe_1_out, T_out_pipe_1, P_atm, A, Cd) ** 2 / (Cd ** 2 * A ** 2 * 2 * Rho_pipe_1_out) #[Pa], Pressure difference

    if delta_P_valve < 0.005 or np.isnan(delta_P_valve):
        delta_P_valve = 0

    #Valve outlet pressure
    P_valve_out = P_pipe_1_out - delta_P_valve
    
    #Torque calculation valve
    F = delta_P_valve * A #[N], Force needed to open valve with certain delta P and area exposed to this
    torque_valve = (F * d_valve) / 2 #[Nm], Torque needed only due to pressure on area
    torque_friction = k_f * W * r_stem #[Nm], Torque needed due to friction
    torque_total = torque_valve + torque_friction #[Nm], Total torque needed for operating the valve (not from zero movement)

    #Hybrid torque, here the pressure difference is being hold, but the friction part is changed by experimental value for
    torque_hybrid = torque_valve + Valve_no_P_Torque

    #Valve toraue via experimental formula
    torques_experimental = 0.0016 * (P_pipe_1_out / 1e5) + Valve_no_P_Torque

    # --- Real gas outlet temperature using isentropic expansion ---
    try:
        s = CP.PropsSI('S', 'T', T_out_pipe_1, 'P', P_pipe_1_out, fluid)
        T_out = CP.PropsSI('T', 'P', P_valve_out, 'S', s, fluid)
    except ValueError:
        # Fallback to ideal gas approximation if CoolProp fails
        if P_atm <= P * (2 / (gamma + 1)) ** (gamma / (gamma - 1)):
            T_out = T * (2 / (gamma + 1))
        else:
            T_out = T * (P_atm / P_valve_out) ** ((gamma - 1) / gamma)

    T_drop_valve = T_out_pipe_1 - T_out


    #Pipe 2
    mu_pipe_2 = mu
    
    #The density depends on the temperature and pressure inside of the pipe, so this will change
    Rho_pipe_2_out = P_pipe_1_out / (R * T_out)
    V_pipe_2_out = mdot / (Rho_pipe_2_out * A_pipe_2)
    
    #Pressure drop through pipe 2
    def reynolds_number_pipe_2(Rho_pipe_2_out, P_valve_out, D_inner_pipe_2, mu_pipe_2):
        """Calculate Reynolds number"""
        return (Rho_pipe_2_out * V_pipe_2_out * D_inner_pipe_2) / mu_pipe_2

    def friction_factor_pipe_2(Re_pipe_2, D_inner_pipe_2, epsilon_pipe_2):
        """Calculate friction factor using laminar formula or Colebrook equation for turbulent flow"""
        if Re_pipe_2 < 2300:
            return 64 / Re_pipe_2  # Laminar flow
        else:
            # Colebrook equation (implicit, solved iteratively)
            def colebrook_pipe_2(f_2):
                arg_2 = epsilon_pipe_2 / (3.7 * D_inner_pipe_2) + 2.51 / (Re_pipe_2 * math.sqrt(f_2))
                if arg_2 <= 0:
                    raise ValueError("Invalid argument for log10 in Colebrook equation")
                return 1 / math.sqrt(f_2) + 2.0 * math.log10(arg_2)
            
            # Initial guess and solving using Newton-Raphson method
            f_guess_2 = 0.05  # Safer initial guess
            for _ in range(20):  # Iterate to refine f
                try:
                    f_new_2 = f_guess_2 - colebrook_pipe_2(f_guess_2) / (-0.5 / (f_guess_2 ** 1.5))
                    if abs(f_new_2 - f_guess_2) < 1e-6:
                        break
                    f_guess_2 = f_new_2
                except ValueError:
                    print("Warning: Could not solve Colebrook equation, returning estimated friction factor.")
                    return 0.02  # Return a default reasonable value if the solver fails
            return f_guess_2

    def pressure_drop_pipe_2(Rho_pipe_2_out, V_pipe_2_out, D_inner_pipe_2, L_pipe_2, mu_pipe_2, epsilon_pipe_2):
        """Calculate pressure drop using Darcy-Weisbach equation"""
        Re_pipe_2 = reynolds_number_pipe_2(Rho_pipe_2_out, V_pipe_2_out, D_inner_pipe_2, mu_pipe_2)
        f_2 = friction_factor_pipe_2(Re_pipe_2, D_inner_pipe_2, epsilon_pipe_2)
        delta_P_pipe_2 = (f_2 * (L_pipe_2 / D_inner_pipe_2) * (Rho_pipe_2_out * V_pipe_2_out ** 2) / 2)
        return delta_P_pipe_2, Re_pipe_2, f_2

    # Compute pressure drop in pipe 2
    delta_P_pipe_2, Re_pipe_2, f = pressure_drop_pipe_2(Rho_tank, V_tank_out, D_inner_pipe_1, L_pipe_1, mu_pipe_1, epsilon_pipe_1)
    #Pressure at outlet of pipe 2
    P_pipe_2_out = P_valve_out - delta_P_pipe_2
    #Density of air at outlet of the pipe 1
    Rho_pipe_2_out = P_pipe_1_out / (R * T_out) #Assume T has not changed yet!!!
    #Velocity at outlet of pipe 2
    V_pipe_2_out = mdot / (Rho_pipe_2_out * A_pipe_2) #[m/s], Speed of flow through pipe 1
    #Velocity drop over pipe 1
    V_drop_pipe_2 = V_pipe_2_out - V_pipe_2_out


    def update_temperature_pipe_2(Rho_pipe_2_out, V_pipe_2_out, D_inner_pipe_2, L_pipe_2, mu, T_out, T_wall_pipe_2, mdot, fluid):
        # Reynolds number
        Re_pipe_2 = (Rho_pipe_2_out * V_pipe_2_out * D_inner_pipe_2) / mu_pipe_2
    
        # Friction factor
        if Re_pipe_2 < 2300:
            f_pipe_2 = 64 / Re_pipe_2  # Laminar flow
        else:
            # Colebrook-White equation solver (approximation)
            epsilon_pipe_2 = 0.00007  # m (pipe roughness) <-- you could also pass this in
            relative_roughness_pipe_2 = epsilon_pipe_2 / D_inner_pipe_2
            f_guess_pipe_2 = 0.02
            for _ in range(10):
                f_guess_pipe_2 = (-2 * np.log10(relative_roughness_pipe_2/3.7 + 2.51/(Re_pipe_2*np.sqrt(f_guess_pipe_2))))**-2
            f_pipe_2 = f_guess_pipe_2
    
        # Pressure drop
        delta_P_pipe_2 = f_pipe_2 * (L_pipe_2 / D_inner_pipe_2) * (Rho_pipe_2_out * V_pipe_2_out**2) / 2
    
        # Temperature change (simplified)
        cp_pipe_2 = CP.PropsSI("CPMASS", "T", T_out, "P", P_atm, fluid)   # Specific heat capacity at constant pressure
        m_pipe_2 = Rho_pipe_2_out * (np.pi * (D_inner_pipe_2/2)**2) * L_pipe_2  # mass inside pipe section
    
        q_dot_pipe_2 = 1.5  # W/m²·K estimate --> Match with experimental setup !
        A_inner_pipe_2 = np.pi * D_inner_pipe_2 * L_pipe_2  # Internal surface area
        heat_transfer_pipe_2 = q_dot_pipe_2 * A_inner_pipe_2 * (T_wall_pipe_2 - T_out)  # Heat gain or loss
    
        if m_pipe_2 > 0 and cp_pipe_2 > 0:
            dT_pipe_2 = heat_transfer_pipe_2 / (m_pipe_2 * cp_pipe_2)
        else:
            dT_pipe_2 = 0.0
    
        T_out_pipe_2 = T_out + dT_pipe_2
    
        # Nusselt number
        if Re_pipe_2 < 2300:
            Nu_pipe_2 = 3.66  # Constant for fully developed laminar flow
        else:
            Nu_pipe_2 = 0.023 * (Re_pipe_2**0.8) * (0.7**0.3)  # Dittus-Boelter for turbulent flow, Pr ≈ 0.7 for air
    
        return T_out_pipe_2, heat_transfer_pipe_2, q_dot_pipe_2, Re_pipe_2, Nu_pipe_2, dT_pipe_2
    
    
    T_out_pipe_2, h_pipe_2, q_dot_pipe_2, Re_pipe_2, Nu_pipe_2, T_drop_pipe_2 = update_temperature_pipe_2(Rho_pipe_2_out, V_pipe_2_out, D_inner_pipe_2, L_pipe_2, mu, T_out, T_wall_pipe_2, mdot, fluid)




    # Store of calculated values
    #General store arrays
    areas.append(A * 1e6)
    torques_theo.append(torque_total)
    torques_exp.append(torques_experimental)
    torques_hyb.append(torque_hybrid)
    T_outlet.append(T_out)
    angle.append(angle)

    
    times.append(time)
    mdots.append(mdot)
    
    #Tank store arrays
    Tank_pressure.append(P / 1e5)
    Tank_temperature.append(T)
    Tank_mass.append(m)
    Tank_v_out.append(V_tank_out)
    #Tank_Re_out.append()
    Tank_M_out.append(M_tank_out)
    Tank_a_out.append(a_tank_out)
    Tank_rho.append(Rho_tank)

    #Pipe 1 store arrays
    Pipe_1_P_out.append(P_pipe_1_out / 1e5)
    Pipe_1_P_drop.append(delta_P_pipe_1 / 1e5)
    Pipe_1_T_out.append(T_out_pipe_1)
    Pipe_1_T_drop.append(T_drop_pipe_1)
    Pipe_1_v_out.append(V_pipe_1_out)
    Pipe_1_v_drop.append(V_drop_pipe_1)
    Pipe_1_Re_out.append(Re_pipe_1)
    Pipe_1_rho_out.append(Rho_pipe_1_out)

    #Valve store arrays
    Valve_P_out.append(P_valve_out / 1e5)
    Valve_P_drop.append(delta_P_valve / 1e5)
    Valve_T_out.append(T_out)
    Valve_T_drop.append(T_drop_valve)
    #Valve_v_out.append()
    #Valve_v_drop.append()
    #Valve_Re_out.append()
    
    #Pipe 2 store arrays
    Pipe_2_P_out.append(P_pipe_2_out / 1e5)
    Pipe_2_P_drop.append(delta_P_pipe_2 / 1e5)
    Pipe_2_T_out.append(T_out_pipe_2)
    Pipe_2_T_drop.append(T_drop_pipe_2)
    Pipe_2_v_out.append(V_pipe_2_out)
    Pipe_2_v_drop.append(V_drop_pipe_2)
    Pipe_2_Re_out.append(Rho_pipe_2_out)




# --- Plotting ---


plt.figure(figsize=(8, 5))
plt.plot(times, Tank_mass, label="Tank mass")
plt.title("Medium Mass vs Flow Time")
plt.ylabel("Mass [kg]")
plt.xlabel("Time [s]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(times, mdots, label="Mass flow rate")
plt.title("Mass Flow Rate vs Flow Time")
plt.ylabel("Mass flow rate [kg/s]")
plt.xlabel("Time [s]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(times, areas, label="Vavle area")
plt.title("Valve Area vs Flow Time")
plt.ylabel("Area [mm2]")
plt.xlabel("Time [s]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(times, torques_theo, label="Valve Torque Theoretical")
plt.plot(times, torques_exp, label="Valve Torque Experimental")
plt.plot(times, torques_hyb, label="Valve Torque Hybrid")
plt.title("Valve Torques vs Flow Time")
plt.ylabel("Valve Torque [Nm]")
plt.xlabel("Time [s]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(times, Tank_pressure, label="Tank Pressure")
plt.plot(times, Pipe_1_P_out, label="Pipe 1 Outlet Pressure")
plt.plot(times, Valve_P_out, label="Valve Outlet Pressure")
plt.plot(times, Pipe_2_P_out, label="Pipe 2 Outlet Pressure")
plt.title("Pressures vs Flow Time")
plt.ylabel("Pressure [Bar]")
plt.xlabel("Time [s]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(times, Tank_temperature, label="Tank Temperature")
plt.plot(times, Pipe_1_T_out, label="Pipe 1 Outlet Temperature")
plt.plot(times, Valve_T_out, label="Valve Outlet Temperature")
plt.plot(times, Pipe_2_T_out, label="Pipe 2 Outlet Temperature")
plt.title("Temperature vs Flow Time")
plt.ylabel("Temperature [K]")
plt.xlabel("Time [s]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


"""

plt.figure(figsize=(12, 10))

plt.subplot(3, 2, 1)
plt.plot(times, Tank_pressure)
plt.ylabel("Pressure in tank [bar]")
plt.xlabel("Time [s]")
plt.grid(True)

plt.subplot(3, 2, 2)
plt.plot(times, Tank_temperature)
plt.ylabel("Temperature in tank [K]")
plt.xlabel("Time [s]")
plt.grid(True)

plt.subplot(3, 2, 3)
plt.plot(times, Tank_mass)
plt.ylabel("Mass in tank [kg]")
plt.xlabel("Time [s]")
plt.grid(True)

plt.subplot(3, 2, 4)
plt.plot(times, mdots)
plt.ylabel("Mass flow rate [kg/s]")
plt.xlabel("Time [s]")
plt.grid(True)

plt.subplot(3, 2, 5)
plt.plot(times, areas)
plt.ylabel("Valve Area [mm²]")
plt.xlabel("Time [s]")
plt.grid(True)

plt.subplot(3, 2, 6)
plt.plot(times, torques)
plt.ylabel("Valve Torque [Nm]")
plt.xlabel("Time [s]")
plt.grid(True)

plt.tight_layout()
plt.show()


plt.figure(figsize=(6, 4))
plt.plot(times, T_outlet)
plt.title("Gas Temperature After Valve (Real Gas)")
plt.ylabel("Outlet Temperature [K]")
plt.xlabel("Time [s]")
plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(6, 4))
plt.plot(times, Pipe_1_T_out)
plt.title("Gas Temperature outlet pipe 1 (Real Gas)")
plt.ylabel("Outlet Temperature [K]")
plt.xlabel("Time [s]")
plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(6, 4))
plt.plot(times, Tank_v_out)
plt.plot(times, Pipe_1_v_out)
plt.plot(times, Pipe_2_v_out)
plt.title("Gas velocity pipe 2 out (Real Gas)")
plt.ylabel("Gas velocity pipe 2 out [m/s]")
plt.xlabel("Time [s]")
plt.grid(True)
plt.tight_layout()
plt.show()


plt.figure(figsize=(6, 4))
plt.plot(times, Tank_M_out)
plt.title("Mach number flow tank")
plt.ylabel("Mach number [/]")
plt.xlabel("Time [s]")
plt.grid(True)
plt.tight_layout()
plt.show()



# Create a figure and first y-axis
fig, ax1 = plt.subplots()

# Plot mass flow rate on the left y-axis
ax1.plot(times, Tank_M_out, 'b-', label='Mach number [/]')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Mach number [/]', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create a second y-axis sharing the same x-axis
ax2 = ax1.twinx()

# Plot tank pressure on the right y-axis
ax2.plot(times, Tank_pressure, 'r--', label='Tank Pressure [bar]')
ax2.set_ylabel('Tank Pressure [bar]', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Add title and grid
plt.title('Mach number and Tank Pressure Over Time')
fig.tight_layout()
plt.grid(True)
plt.show()


# Create a figure and first y-axis
fig, ax1 = plt.subplots()

# Plot mass flow rate on the left y-axis
ax1.plot(times, Tank_M_out, 'b-', label='Mach number [/]')
ax1.set_xlabel('Time [s]')
ax1.set_ylabel('Mach number [/]', color='b')
ax1.tick_params(axis='y', labelcolor='b')

# Create a second y-axis sharing the same x-axis
ax2 = ax1.twinx()

# Plot tank pressure on the right y-axis
ax2.plot(times, Tank_temperature, 'r--', label='Tank Temperature [K]')
ax2.set_ylabel('Tank Temperature [K]', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Add title and grid
plt.title('Mach number and Tank Temperature Over Time')
fig.tight_layout()
plt.grid(True)
plt.show()


"""
