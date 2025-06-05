# -*- coding: utf-8 -*-
"""
Created on Thu May 22 18:42:05 2025

@author: niels
"""

#This model is based on the technical paper: Modeling the nitrous run tank emptying, from Aspire space
#That model was purely based on a blowdown tank, our is pressure regulated and supercharged, so some workarounds and extra's were needed

# Importing needed libs
import CoolProp.CoolProp as CP  # Importing CoolProp
import matplotlib.pyplot as plt  # Importing matplotlib
import numpy as np  # Importing numpy
import pandas as pd #Importing pandas


# --- Tank emptying model --- #

# Input parameters
timestep = 0.0001  # [s], Simulation time step
min_temperature_n2o = 182.23  # [K], Restriction from CoolProp Package, beacause it becomes solid at that temperature!
tank_temperature = 283.15  # [K], Tank temperature, typically close to the atmospheric temperature (depends on time between tank fill and tank emptying)
combustion_chamber_pressure = 20e5  # [Pa], Pressure of the combustion chamber
pressure_drop_inj = 20e5  # [Pa], Pressure drop over the injector (Needed for the performance of the injector)
vapourizing_factor_n2o = 1.0575 # [/], Factor that multiplies the N2O amount, so this extra part will vapourize to keep the pressure
Burn_time = 6.535386021327031 # [s], Duration of the burn
mass_flow_N2O = 0.465 # [kg/s], We go for a constant mass flow rate, because the valve will limit this to the wanted flow
Transient_flow_factor = 1.05 # [/], Factor so some more mass is taken into the tank, this is mainly for transient phases, NOT for vapourizing

target_pressure = 80e5  # [Pa], Supercharged pressure, can also be seen as the maximum pressure of the weakest part
molar_mass_n2 = 0.0280134  # [kg/mol], Molar masss of nitrogen
R = 8.314  # [J/(mol*K)], Universal gas constant 

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


# --- Liquid emptying model ---


# Simulation variables
liquid_mass = initial_liquid_mass  # [kg], liquid mass
vapour_mass_N2 = mass_min_in_inj_N2  # [kg], vapour mass, there is no N2O vapour, because it is above the vapour pressure until a certain point
time = 0  # [s], starting time of model
total_vaporized_mass = 0  # Initialize accumulator

# Arrays for storing results
time_array = []
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


    # Update time
    time += timestep

    #Storing of the results
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


#Printing some essential values

print(f"Tank Volume: {V_tank_N2O_N2:.5f} m3")
print(f"Amount of N2 gas: {mass_min_in_inj_N2:.4f} kg")
print(f"Total N2O vaporized to maintain pressure: {total_vaporized_mass:.4f} kg")
print(f"Fraction of initial N2O mass that was vaporized: {100 * total_vaporized_mass / initial_liquid_mass:.2f}%")

average_tank_density = np.mean(n2o_liq_density_array)
print(f"Average Tank Pressure during burn: {average_tank_density:.2f} kg/m3")



# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_array, tank_pressure_array, label='Tank Pressure (Pa)')
plt.xlabel('Time (s)')
plt.ylabel('Tank Pressure [Pa]')
plt.title('Nitrous Oxide Tank Dynamics')
plt.grid()
plt.legend()
plt.show()

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_array, n2_pressure_new_array, label='N2 Pressure (Pa)')
plt.xlabel('Time (s)')
plt.ylabel('N2 Pressure [Pa]')
plt.title('Nitrous Oxide Tank Dynamics')
plt.grid()
plt.legend()
plt.show()

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_array, n2o_pressure_vapour_array, label='N2O Pressure (Pa)')
plt.xlabel('Time (s)')
plt.ylabel('N2O Pressure [Pa]')
plt.title('Nitrous Oxide Tank Dynamics')
plt.grid()
plt.legend()
plt.show()


# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_array, tank_temperature_array, label='Tank Temperature (K)')
plt.xlabel('Time (s)')
plt.ylabel('Tank Temperature [K]')
plt.title('Nitrous Oxide Tank Dynamics')
plt.grid()
plt.legend()
plt.show()


# Plotting
plt.figure(figsize=(10, 6))
plt.plot(time_array, n2o_liq_density_array, label='N2O liquid density [kg/m3]')
plt.xlabel('Time (s)')
plt.ylabel('N2O liquid density [kg/m3]')
plt.title('Nitrous Oxide Tank Dynamics')
plt.grid()
plt.legend()
plt.show()



#Data saving in csv
# Create a DataFrame with the relevant data
data = {
    'Time [s]': time_array,
    'Tank Pressure [Pa]': tank_pressure_array
}

df = pd.DataFrame(data)

# Save to CSV
df.to_csv('tank_pressure_data.csv', index=False)

print("CSV file 'tank_pressure_data.csv' saved successfully.")





# === STRUCTURAL DESIGN VERIFICATION ===
#(See thesis section about oxidizer tank for cites)

# Load the Excel file
file_path = "Tank_Design_Database.xlsx"
df_bolts = pd.read_excel(file_path, sheet_name=1, engine="openpyxl")

# Extract and filter bolt diameters (2.9 mm to 14 mm)
bolt_diameters = df_bolts["Diameter_Stress_Cross_Section_[m]"].values
bolt_diameters = bolt_diameters[(bolt_diameters >= 0.0033) & (bolt_diameters <= 0.0096)]

# Tank parameters
D_out = 0.140  # Outer diameter [m]
t_wall = 0.005  # Wall thickness [m]
D_in = D_out - 2 * t_wall  # Inner diameter [m]
MEOP = 80e5  # Maximum Expected Operating Pressure [Pa]
SF = 2  # Safety Factor
P_SF = MEOP * SF  # Pressure with safety factor [Pa]
Tensile_strength_alu = 215 * 1e6  # Tensile strength of aluminum [Pa]
Shear_strength_alu = 130 * 1e6 # Shear strength of aluminum [Pa]
Ullage =  vapourizing_factor_n2o - 1 # Ullage factor (20%)
Oxi_mass = initial_liquid_mass  # Oxidizer mass [kg]
Density_alu = 2700  # Aluminum density [kg/m^3]
Bolt_shear = 800 * 1e6  # Bolt shear strength [Pa]
Density_Bolt_Mat = 7850 # Density of the material from the bolts [kg/m3]

#Bulkhead material specs
#Aluminium EN AW-6082 T6 https://www.aluminiumopmaatgemaakt.be/aluminium-rondstaf.html chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://xometry.eu/wp-content/uploads/2020/09/Datasheet-EN-AW-6082.pdf
Yield_stress_6082_T6 = 260 * 1e6 #[Pa], yield stress of aluminium 6082 T6

#2:1 Elliptical bulkhead
K = 1.0 #[/], shape factor for type of elliptical
E = 1.0 #[/], Welding factor, no welds so = 1.0

#Here you can change the amount of bolt circles (2 is more than enough)
Bolt_Circles = 2
Arbitary_Spacing = 0.015 #[m], see descroption at code lines of E_min

# Storage arrays for results
results = []
results_MC = []

# Loop through bolt sizes and numbers
for Number_of_bolts in range(20, 51, 1):  # Varying bolt numbers (10 to 40)
    for d_bolt in bolt_diameters:
            # Calculate bolt shear stress
            Bolt_Shear_Calc = ((np.pi / 4) * D_in ** 2 * P_SF) / (Number_of_bolts * (np.pi / 4) * d_bolt ** 2)
            
            # Force per bolt
            F_bolt = (((np.pi) / 4) * D_in**2 * P_SF) / Number_of_bolts  
            
            # Edge distance considerations for multiple circles
            if Bolt_Circles == 1:
                E_min = max(F_bolt / (Shear_strength_alu * 2 * t_wall), 2 * d_bolt)  
            else:
                #Overlap of area is still missing in this code!!! (Page 42), so be sure they don't overlap!
                E_min_1 = max(F_bolt / (Shear_strength_alu * 2 * t_wall), 2 * d_bolt)  
                E_min_2 = E_min_1 + Arbitary_Spacing  # Arbitrary spacing between bolt circles, the second row of bolts already has enough material between the bolt, and the edge, so a chosen spacing can be used.
                #You can even say, arbitrary spacing is zero, theoretically, okay! But the bolt heads would collide during assembly!
                E_min = (E_min_1 + E_min_2) / 2  # Average
            
            # Casing Tensile Strength
            if Bolt_Circles == 1:
                Tensile_Strength_Tube = ((np.pi / 4) * D_in**2 * P_SF) / (((D_out - t_wall) * np.pi - Number_of_bolts * d_bolt) * t_wall)
            else:
                Tensile_Strength_Tube = ((np.pi / 4) * D_in**2 * P_SF) / (((D_out - t_wall) * np.pi - Number_of_bolts / Bolt_Circles * d_bolt) * t_wall)
            
            # Bearing Stress
            Bearing_Stress = F_bolt / (d_bolt * t_wall)
            
            # Store results
            if Bolt_Circles == 1:
                results.append([d_bolt, Number_of_bolts, Bolt_Shear_Calc, E_min, Tensile_Strength_Tube, Bearing_Stress])
            else:
                results_MC.append([d_bolt, Number_of_bolts, Bolt_Shear_Calc, E_min, Tensile_Strength_Tube, Bearing_Stress])

# Convert results to DataFrames
df_results = pd.DataFrame(results, columns=["Bolt Diameter (m)", "Number of Bolts", "Bolt Shear Calc (Pa)", "E_min (m)", "Tensile Strength (Pa)", "Bearing Stress (Pa)"])
df_results_MC = pd.DataFrame(results_MC, columns=["Bolt Diameter (m)", "Number of Bolts", "Bolt Shear Calc (Pa)", "E_min (m)", "Tensile Strength (Pa)", "Bearing Stress (Pa)"])

# Save results
df_results.to_excel("Bolt_Influence_Analysis.xlsx", index=False)
df_results_MC.to_excel("Bolt_Influence_Analysis_MC.xlsx", index=False)

# **Visualization**
def plot_results(df, title_suffix, include_threshold=True):
    parameters = ["Bolt Shear Calc (Pa)", "Tensile Strength (Pa)", "Bearing Stress (Pa)", "E_min (m)"]
    y_labels = ["Bolt Shearing (Pa)", "Tensile Strength (Pa)", "Bearing Stress (Pa)", "E_min (m)"]
    titles = ["Influence of Bolt Quantity on Bolt Shearing", "Influence of Bolt Quantity on Tensile Strength", "Influence of Bolt Quantity on Bearing Stress", "Influence of Bolt Quantity on E_min"]
    thresholds = [Bolt_shear, Tensile_strength_alu, Tensile_strength_alu, None]
    
    
    for param, y_label, title, threshold in zip(parameters, y_labels, titles, thresholds):
        plt.figure(figsize=(10, 6))
        for d_bolt in np.unique(df["Bolt Diameter (m)"]):
            subset = df[df["Bolt Diameter (m)"] == d_bolt]
            plt.plot(subset["Number of Bolts"], subset[param], label=f"Diameter {d_bolt*1e3:.1f} mm")
        if include_threshold and threshold is not None:
            plt.axhline(y=threshold, color='r', linestyle='--', label='Threshold')
        plt.xlabel("Number of Bolts")
        plt.ylabel(y_label)
        plt.title(f"{title} {title_suffix}")
        plt.legend()
        plt.grid()
        plt.show()

if Bolt_Circles == 1:
    plot_results(df_results, "(Single Bolt Circle)")
else:
    plot_results(df_results_MC, "(Multiple Bolt Circles)")


#Mass tank (Estimation)
#Mass bulkheads:

#Minimum thickness of flat bulkhead(s)
Thickness_min_bulkhead_flat = np.sqrt((3 * P_SF * (D_in / 2)**2) / (4 * Yield_stress_6082_T6))
Mass_Bulkhead_flat = ((np.pi * D_in ** 2) / 4) * Thickness_min_bulkhead_flat * Density_alu

#Minumun thickness of 2:1 elliptical bulkhead(s)
Thickness_min_bulkhead_elliptical = (P_SF * D_in) / (2 * K * Yield_stress_6082_T6 * E - 0.2 * P_SF)
Volume_Bulkhead_elliptical = ((np.pi * D_in ** 3) / 24) - ((np.pi * (D_in - Thickness_min_bulkhead_elliptical) ** 3) / 24)
Mass_Bulkhead_elliptical = Volume_Bulkhead_elliptical * Density_alu

#Minimum thickness of hemispherical bulkhead(s)
Thickness_min_bulkhead_hemispherical = (P_SF * (D_in / 2)) / (2 * Yield_stress_6082_T6)
Volume_Bulkhead_hemispherical = ((2 / 3) * np.pi * (D_in / 2) ** 3) - ((2 / 3) * np.pi * ((D_in - Thickness_min_bulkhead_hemispherical) / 2) ** 3)
Mass_Bulkhead_hemispherical = Volume_Bulkhead_hemispherical * Density_alu

#Tube for flat bulkhead
#Length of the tube (only for oxidzer, so inner dimension, between inner faces of the bulkheads)
L_tube_Flat = V_tank_N2O_N2 / ((np.pi * D_in**2) / 4)
#Mass of the tube
Mass_Tube_Flat = ((((np.pi * D_out ** 2) / 4) - ((np.pi * D_in ** 2) / 4)) * L_tube_Flat) * Density_alu # [kg], Mass of tube only

#Tube for elliptical bulkhead
#Length of the tube
V_tank_Elliptical = V_tank_N2O_N2 - ((np.pi * (D_in - Thickness_min_bulkhead_elliptical) ** 3) / 24)
L_tube_Elliptical = V_tank_Elliptical / ((np.pi * D_in**2) / 4)
#Mass of the tube
Mass_Tube_Elliptical = ((((np.pi * D_out ** 2) / 4) - ((np.pi * D_in ** 2) / 4)) * L_tube_Elliptical) * Density_alu # [kg], Mass of tube only

#Tube for hemispherical bulkhead
#Length of tube
V_tank_Hemispherical = V_tank_N2O_N2 - ((2 / 3) * np.pi * ((D_in - Thickness_min_bulkhead_hemispherical) / 2) ** 3)
L_tube_Hemispherical = V_tank_Hemispherical / ((np.pi * D_in**2) / 4)
#Mass of the tube
Mass_Tube_Hemispherical = ((((np.pi * D_out ** 2) / 4) - ((np.pi * D_in ** 2) / 4)) * L_tube_Hemispherical) * Density_alu # [kg], Mass of tube only


#Bolt mass estimatiom calculation:
#This part will calculate the bolt mass with a chosen bolt type
#For now fill in the bolt diameter yourself after investigating the approtiate bolt from the curves
Bolt_Diameter_Chosen = 0.006 # Bolt Diameter [m]
Bolt_Length_Chosen = 0.010 # Bolt length chosen [m]
Bolt_Head_Diameter_Chosen = 0.010 # Bolt head diameter chosen [m]
Bolt_Head_Height_Chosen = 0.006 # Bolt head height chosen [m]
Amount_Bolts_Chosen = 80 # This highly depends on the configuration you are using!
Total_Mass_Bolt = (((np.pi * Bolt_Diameter_Chosen ** 2) / 4) * Bolt_Length_Chosen) * Density_Bolt_Mat + (((np.pi * Bolt_Head_Diameter_Chosen ** 2) / 4) * Bolt_Head_Height_Chosen) * Density_Bolt_Mat
All_Bolts_Mass = Total_Mass_Bolt * Amount_Bolts_Chosen

#Total estimation mass tank
Total_Mass_Tank_bulkhead_flat = Mass_Tube_Flat + 2 * Mass_Bulkhead_flat + All_Bolts_Mass #[kg], Estimation of total mass of run tank
Total_Mass_Tank_bulkhead_elliptical = Mass_Tube_Elliptical + 2 * Mass_Bulkhead_elliptical + All_Bolts_Mass #[kg], Estimation of total mass of run tank
Total_Mass_Tank_bulkhead_hemispherical = Mass_Tube_Hemispherical + 2 * Mass_Bulkhead_hemispherical + All_Bolts_Mass #[kg], Estimation of total mass of run tank

#Printing solutions with descriptions
print(f"Minimum Flat Bulkhead Thickness: {Thickness_min_bulkhead_flat:.5f} m")
print(f"Minimum Elliptical Bulkhead Thickness: {Thickness_min_bulkhead_elliptical:.5f} m")
print(f"Minimum Hemispherical Bulkhead Thickness: {Thickness_min_bulkhead_hemispherical:.5f} m")
print(f"Minimum length inside of the tank flat bulkheads: {L_tube_Flat:.4f} m")
print(f"Minimum length inside of the tank elliptical bulkheads: {L_tube_Elliptical:.4f} m")
print(f"Minimum length inside of the tank hemispherical bulkheads: {L_tube_Hemispherical:.4f} m")
print(f"Mass of Flat Bulkhead: {Mass_Bulkhead_flat:.3f} kg")
print(f"Mass of Elliptical Bulkhead: {Mass_Bulkhead_elliptical:.3f} kg")
print(f"Mass of Hemispherical Bulkhead: {Mass_Bulkhead_hemispherical:.3f} kg")
print(f"Mass of Tube: {Mass_Tube_Flat:.2f} kg")
print(f"Mass of all bolts : {All_Bolts_Mass:2f} kg")
print(f"Total Estimated Tank Mass + flat bulkheads: {Total_Mass_Tank_bulkhead_flat:.2f} kg")
print(f"Total Estimated Tank Mass + elliptical bulkheads: {Total_Mass_Tank_bulkhead_elliptical:.2f} kg")
print(f"Total Estimated Tank Mass + hemispherical bulkheads: {Total_Mass_Tank_bulkhead_hemispherical:.2f} kg")
print(f"Be aware that the masses of the bulkheads are just the pressure holding thickness, so there is no material for the bolts and O-rings taken into account!")


#Minimum wall thickness check of hoop stress and axial stress
#Hoop Stress
Hoop_Stress = (P_SF * (D_in / 2)) / t_wall # [Pa], Hoop stress of tank
if Hoop_Stress < Tensile_strength_alu:
    print(f"Hoop stress is OKI: {Hoop_Stress:.2f} Pa")
else:
    print(f"Hoop stress is to high!: {Hoop_Stress:.2f} Pa > {Tensile_strength_alu} Pa")
#Axial Stress
Axial_Stress = (P_SF * (D_in / 2)) / (2 * t_wall)
if Axial_Stress < Tensile_strength_alu:
    print(f"Axial stress is OKI: {Axial_Stress:.2f} Pa")
else:
    print(f"Axial stress is to high!: {Axial_Stress:.2f} Pa > {Tensile_strength_alu} Pa")


