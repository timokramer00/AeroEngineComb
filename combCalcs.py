import numpy as np
import matplotlib.pyplot as plt
from polyLib import CO2, water, N2, fuel
from dataclasses import dataclass


class Massflows:
    O2: float
    CO2: float
    fuel: float
    N2: float
    water: float


class Combustor:
    def __init__(self, form_enthalpy, eng_spec, molar_mass, volume=0.25, x1 = 0.07, x2 = 0.08, combustion_eff = 0.995):
        self.cp_gas = 1150  # J/kg
        self.cp_air = 1000  # J/kg
        self.FAR_stoich = 0.682
        self.fuel_temp = 300  # K
        self.calvalue = (48/4 * form_enthalpy["CO2"] + 46/4 * form_enthalpy["H2O"] - form_enthalpy["fuel"]) / molar_mass["fuel"]  # J/kg
        self.eng_spec = eng_spec
        self.molar_mass = molar_mass
        self.volume = volume  # m3
        self.ox_frac = 0.2
        self.comb_eff = combustion_eff
        self.mdots = Massflows(0, 0, 0, 0, 0)
        
        x3 = 0.42 - x1 - x2 - 0.1
        self.zonal_admittance = [0.16, 2*x1, 2*x2, 2*x3 + 0.2] 
        self.stations_positions = [0.05, 0.07, 0.09, 0.15, 0.25]
        self.stations = ["Injector", "Primary Zone", "Secondary Zone", "Dilution Zone"]
        
    def plot_all(self):
        # TODO modify to add equiv ratio to plots
        plt.grid()
        plt.xlabel("Axial position [m]")
        plt.ylabel("Temperature [K]")
        for thrust in self.eng_spec["thrust"]:
            temps = self.compute_temp_total_equiv(thrust)
            print("temps", temps)
            print("stations", self.stations)
            self.add_to_plot(thrust, temps, 69)
        plt.legend()
        plt.show()
        
    def compute_temp_total_equiv(self, thrust):
        air_in, press_in, temp_in, temp_out, self.mdots.fuel = self.get_conditions(thrust)
        
        temp = temp_in
        temperatures = [temp]
        
        for airfrac in self.zonal_admittance:
            # get added air
            self.mdots.N2 += (1 - self.ox_frac) * airfrac * air_in
            self.mdot.O2 = self.ox_frac airfrac * air_in
            
            temp = self.compute_burn(temp)
            
            temperatures.append(temp)
        #TODO implement equivalence ratio
        return temperatures
    
    def compute_burn(self, mdot_O2, mdot_fuel, temp_in):
        moles_fuel = mdot_fuel / self.molar_mass["fuel"]
        moles_ox = mdot_O2 / self.molar_mass["O2"]
        # mol ratios
        o2_per_fuel = 81/4
        CO2_per_fuel = 0
        water_per_fuel = 0
        
        
        # TODO implement iterative solver
        # fuel rich
        if moles_ox/moles_fuel <= o2_per_fuel:
            # compute new mole values
            pass
        
        # ox rich
        else:
            pass
                    
                             
        return mdot_fuel, new_co2, new_water, mdot_O2, product_temp
    
    def calculate_temperature(self, temperature, current_air_mass, new_air, energy_in, temp_in):
        
        return
    
    def add_to_plot(self, thrust, temperature, total_equivalence_ratio):
        # TODO modify to add equiv ratio
        plt.plot(self.stations_positions, temperature, label=f"{thrust} kN")
        plt.scatter(self.stations_positions, temperature)
        
    
    def get_conditions(self, thrust):
        spec = self.eng_spec
        idx = spec["thrust"].index(thrust)
        return spec["mdot_comb_in"][idx], spec["press_comb_in"][idx], spec["temp_comb_in"][idx], spec["temp_comb_out"][idx], spec["mdot_fuel"][idx]
        
        
if __name__ == "__main__":
        
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
        "O2": (15.99*2) * 1e-3,
        "N2": 28.0134*1e-3
        }

    # all in J/mol
    form_enthalpy = {
        "CO2": -393.5e3,
        "H2O": -285.8e3,
        "fuel": -61.7*4184 ,  # c12h23
        "O2": 0
        }
    comb = Combustor(form_enthalpy, eng_spec, molar_mass)
    comb.plot_temp_distr_all_thrust()
