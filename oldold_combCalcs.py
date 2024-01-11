import numpy as np
import matplotlib.pyplot as plt


def plot_flame_temperature(thrust, cp_air=1000, cp_gas=1150, cp_fuel=2000, ox_mass_fraction=0.2):
    mdot_air, press_in, temp_in, temp_out, mdot_fuel = get_conditions(thrust)
    
    mdot_ox = ox_mass_fraction * mdot_air
    fuel_burnt_frac = 0
    for i, station in enumerate(["Inlet", "PZ", "SZ", "DZ"]):
        # calculate how much fuel we burn based on o2 and fuel contents
        # calculate O/F ratio
        # calculate heat released
        # reinterpolate cp of gas based on current fuel_burnt_frac
        # calculate new mix temperature with added air flow
        
def compute_fuel_burn(mdot_ox, mdot_air, mdot_fuel, fuel_burnt_frac):
    #calvalue= 48*HfCO2 + 46*HfH2O - 4*HfC12H23 - 81*HfO2
    # oxygen omitted as it's formation enthalpy is zero
    released_per_fuel_mole = 48/4 * form_enthalpy["CO2"] + 46/4 * form_enthalpy["H2O"] - form_enthalpy["fuel"]
    calvalue = released_per_fuel_mole / molar_mass["fuel"]
    moles_fuel = mdot_fuel * fuel_burnt_frac / molar_mass["fuel"]
    moles_ox = mdot_ox / molar_mass["O2"]
    o2_per_fuel = 81/4  # mol
    
    
    if moles_ox => moles_fuel * o2_per_fuel:
        # fuel lean
        fuel_burnt = moles_ox / o2_per_fuel * molar_mass["fuel"]
        energy_released = -fuel_burnt * calvalue # minus to make it positve
        fuel_burnt_frac += fuel_burnt / mdot_fuel
        print("Ox rich")
        
    else:
        # fuel rich
        fuel_burnt = 
        print("Fuel rich")
    
    if energy_released < 0:
        raise ValueError(f"reaction is endothermic: \nburnt fuel: {fuel_burnt :.3f} kg \ncalval: {calvalue :.3f} j/kg"
    # figure out whether we run out of o2 or fuel first
    return energy_released, mdot_ox, fuel_burnt_frac, equivalence ratio


def get_conditions(thrust, spec=eng_spec)
    idx = spec["thrust"].index(thrust)
    return spec["mdot_comb_in"][idx], spec["press_comb_in"][idx], spec["temp_comb_in"][idx], spec["temp_comb_out"][idx], spec["mdot_fuel"][idx]


eng_spec = {
    "thrust" : [20, 40, 80, 120],   # kn
    "mdot_comb_in" : [10,16,29,39],  #kg/s
    "press_comb_in" : [700e3,1100e3,2200e3,3200e3],  #Pa
    "temp_comb_in" : [540,610,725,820],  #K
    "temp_comb_out" : [1150,1230,1410,1630],  #K
    "mdot_fuel" : [0.24329509704489904, 0.39565366601399976, 0.7923044430814622, 1.259949527155469]  # kg/s
    }

# all in kg/mol
molar_mass = {
    "CO2": (12.011 + 15.99*2) * 1e-3,
    "H2O": (1.0078*2 + 15.99) * 1e-3,
    "fuel": (12*12.011 + 23 * 1.0078) * 1e-3,  # c12h23
    "O2": (15.99*2) * 1e-3
    }

# all in J/mol
form_enthalpy = {
    "CO2": -393.5e3,
    "H2O": -285.8e3,
    "fuel": -61.7*4184 ,  # c12h23
    "O2": 0
    }



injection_start = 0.05
pz_start = 0.07  # m
sz_start = 0.09  # m
dz_start = 0.15  # m
combuster_length = 0.25  # m

x1 = 0.07
x2 = 0.08
x3 = 0.42 - x1 - x2 - 0.1
station_lengths = [pz_start - injection_start, 
                sz_start - pz_start, 
                dz_start - sz_start, 
                combuster_length - dz_start]


# these are the air fractions entering at the injector, PZ, SZ, and DZ
# the 0.2 entering at the end is the cooling air initially injected in the SZ
# can't really cool if it's on fire right?
zonal_admittance = [0.16, 2*x1, 2*x2, 2*x3 + 0.2] 


if __name__ == "__main__":        
    
    for thrust in eng_spec["thrust"]:
        plot_flame_temperature(thrust)
