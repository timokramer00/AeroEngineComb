import numpy as np
import matplotlib.pyplot as plt
from polyLib import CO2, water, N2, fuel, O2
from dataclasses import dataclass


class Massflows:
    O2 = 0
    CO2 = 0
    fuel = 0
    N2 = 0
    water = 0
    
    def __str__(self):
        r = f"O2 {self.O2 :.3f}, \nCO2 {self.CO2 :.3f}, \nfuel {self.fuel :.3f}, \nN2 {self.N2 :.3f}, \nwater {self.water :.3f} \n \n"
        return r
              
              

    molar_mass = {
        "CO2": (12.011 + 15.99*2) * 1e-3,
        "H2O": (1.0078*2 + 15.99) * 1e-3,
        "fuel": (12*12.011 + 23 * 1.0078) * 1e-3,  # c12h23
        "O2": (15.99*2) * 1e-3,
        "N2": 28.0134*1e-3
        }


def massflows_to_moles_funcs(mdots: Massflows):
    moles = np.array([mdots.O2, mdots.CO2, mdots.fuel, mdots.N2, mdots.water])
    funcs = [O2, CO2, fuel, N2, water]
    mole_masses = np.array([molar_mass["O2"], molar_mass["CO2"], molar_mass["fuel"], molar_mass["N2"], molar_mass["H2O"]])
    moles /= mole_masses
    return moles, funcs


def compute_enthalpy(massflows, temp, ref_temp=298):
    component_moles, component_funcs = massflows_to_moles_funcs(massflows)
    total_enth = 0
    for c_moles, c_func in zip(component_moles, component_funcs):
        total_enth += c_moles * (c_func.enthalpy(temp) - c_func.enthalpy(ref_temp))
    return total_enth
        


class Combustor:
    def __init__(self, form_enthalpy, eng_spec, molar_mass, volume=0.25, x1 = 0.07, x2 = 0.08, combustion_eff = 0.995):
        self.cp_gas = 1150  # J/kg
        self.cp_air = 1000  # J/kg
        self.FAR_stoich = 0.682
        self.fuel_temp = 300  # K
        self.calvalue = 81/4 * form_enthalpy["CO2"] + 48/4 * form_enthalpy["H2O"] - form_enthalpy["fuel"]  # J/mol
        self.eng_spec = eng_spec
        self.molar_mass = molar_mass
        self.volume = volume  # m3
        self.ox_frac = 0.21
        self.n_iter = 2500
        self.comb_eff = combustion_eff
        self.mdots = None
        
        x3 = 0.42 - x1 - x2 - 0.1
        self.zonal_admittance = [0.16, 2*x1, 2*x2, 2*x3 + 0.2] 
        self.stations_positions = [0.05, 0.07, 0.09, 0.15, 0.25]
        self.stations = ["Injector", "Primary Zone", "Secondary Zone", "Dilution Zone"]
        
    def plot_all(self):
        # TODO modify to add equiv ratio to plots
        plt.grid()
        plt.xlabel("Axial position [m]")
        plt.ylabel("Temperature [K]")
        for thrust in self.eng_spec["thrust"][:1]:
            temps = self.compute_temp_total_equiv(thrust)
            print("temps", temps)
            print("stations", self.stations)
            self.add_to_plot(thrust, temps, 69)
        plt.legend()
        plt.show()
        
    def compute_temp_total_equiv(self, thrust):
        self.mdots = Massflows()
        air_in, press_in, temp_in, temp_out, self.mdots.fuel = self.get_conditions(thrust)
        
        temp = temp_in
        temperatures = [temp]
        
        for i, airfrac in enumerate(self.zonal_admittance):
            # get added air
            self.mdots.N2 += (1 - self.ox_frac) * airfrac * air_in
            self.mdots.O2 += self.ox_frac * airfrac * air_in
            
            temp = self.compute_burn_equiv_ratio(temp)
            
            temperatures.append(temp)
        return temperatures
    
    def compute_burn_equiv_ratio(self, reactant_temp):
        moles_fuel = self.mdots.fuel / self.molar_mass["fuel"]
        moles_ox = self.mdots.O2 / self.molar_mass["O2"]
        # mol ratios
        o2_per_fuel = 81/4
        CO2_per_fuel = 48/4
        water_per_fuel = 81/4
        
        products = Massflows()  # inits to zero
        # figure out released energy and products
        print(self.mdots)
        
        # fuel rich
        if moles_ox/moles_fuel <= o2_per_fuel:
            print("fuel rich")
            # all O2 gone, so stays zero kg
            fuel_burn = moles_ox / o2_per_fuel
            products.CO2 = fuel_burn * CO2_per_fuel * self.molar_mass["CO2"]
            products.water = fuel_burn * water_per_fuel * self.molar_mass["H2O"]
            products.N2 = self.mdots.N2
            products.fuel = self.mdots.fuel - fuel_burn * self.molar_mass["fuel"]
            
            energy_released = fuel_burn * self.calvalue
            
        
        # ox rich
        else:
            print("ox rich")
            fuel_burn = moles_fuel
            products.O2 = fuel_burn * o2_per_fuel * self.molar_mass["O2"]
            products.CO2 = fuel_burn * CO2_per_fuel * self.molar_mass["CO2"]
            products.water = fuel_burn * water_per_fuel * self.molar_mass["H2O"]
            products.N2 = self.mdots.N2
            products.fuel = self.mdots.fuel - fuel_burn * self.molar_mass["fuel"]
            
            energy_released = fuel_burn * self.calvalue
        
        
        # determine temperature iteratively
        print(products)
        reactant_enthalpy = compute_enthalpy(self.mdots, reactant_temp)
        target = reactant_enthalpy - energy_released
        temp_prod = reactant_temp
    
        running = True
        n_iter = 0
        while running:
            product_enthalpy = compute_enthalpy(products, temp_prod)
            
            rel = product_enthalpy / target
            
            if abs(rel - 1) < 1e-3:
                running = False
            elif rel < 1:
                temp_prod += 1
            elif rel > 1:
                temp_prod -= 1
            #print(rel)
            #print(temp_prod)
            n_iter += 1
            if n_iter > self.n_iter:
                raise ValueError("Max number of iterations exceeded!")
        
        print()
        self.mdots = products
        return temp_prod
    
    
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
    comb.plot_all()
