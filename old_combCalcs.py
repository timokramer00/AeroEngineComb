import numpy as np
import matplotlib.pyplot as plt

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
        
        x3 = 0.42 - x1 - x2 - 0.1
        self.zonal_admittance = [0.16, 2*x1, 2*x2, 2*x3 + 0.2] 
        self.stations_positions = [0.05, 0.07, 0.09, 0.15, 0.25]
        self.stations = ["Injector", "Primary Zone", "Secondary Zone", "Dilution Zone"]
        self.mdot_fuel = None
        
    def plot_temp_distr_all_thrust(self):
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
        air_in, press_in, temp_in, temp_out, self.mdot_fuel = self.get_conditions(thrust)
        
        temp = temp_in
        temperatures = [temp]
        
        current_air = 0
        ox_content = 0
        fuel_burnt_frac = 0
        for airfrac in self.zonal_admittance:
            ox_content += airfrac * air_in * self.ox_frac
            
            fuel_burnt_frac, energy_released, ox_content = self.compute_burn_and_energy_released(fuel_burnt_frac, ox_content)
            
            temp, current_air = self.calculate_temperature(temp, current_air, airfrac*air_in, energy_released, temp_in)
            
            temperatures.append(temp)
        #TODO implement equivalence ratio
        return temperatures
    
    def compute_burn_and_energy_released(self, fuel_burnt_frac, ox_content):
        moles_fuel = self.mdot_fuel * (1-fuel_burnt_frac) / self.molar_mass["fuel"]
        moles_ox = ox_content / self.molar_mass["O2"]
        o2_per_fuel = 81/4  # mol
        
        e = moles_fuel* self.molar_mass["fuel"]
        print(f"cur fuel: {e :.3f}")
    
        if moles_ox >= moles_fuel * o2_per_fuel:
            # fuel limiting
            fuel_burnt = moles_fuel * self.molar_mass["fuel"]
            energy_released = -fuel_burnt * self.calvalue # minus to make it positve
            fuel_burnt_frac += fuel_burnt / self.mdot_fuel
            ox_content -= moles_fuel * o2_per_fuel * self.molar_mass["O2"]
            print("Ox rich")
            print(f"Fuel burnt: {fuel_burnt :.3f}")
            print(f"Fuel burnt frac: {fuel_burnt_frac :.3f}")
            print(f"Ox content: {ox_content :.3f}")
            
        else:
            # fuel rich
            fuel_burnt = moles_ox / o2_per_fuel * self.molar_mass["fuel"]
            energy_released = -fuel_burnt * self.calvalue
            fuel_burnt_frac += fuel_burnt / self.mdot_fuel
            ox_content = 0
            print("Fuel rich")
            print(f"Fuel burnt: {fuel_burnt :.3f}")
            print(f"Fuel burnt frac: {fuel_burnt_frac :.3f}")
            print(f"Ox content: {ox_content :.3f}")
            
        if energy_released < 0:
            raise ValueError(f"reaction is endothermic: \nburnt fuel: {fuel_burnt :.3f} kg \ncalval: {self.calvalue :.3f} j/kg")
                             
        return fuel_burnt_frac, energy_released*self.comb_eff, ox_content
    
    def calculate_temperature(self, temperature, current_air_mass, new_air, energy_in, temp_in):
        # total energy
        energy_convected_in = temperature * current_air_mass*self.cp_gas + new_air * temp_in*self.cp_air
        air_out = current_air_mass + new_air
        print(f"air out: {air_out :.3f}")
        print(f"energy in {energy_convected_in :.3f}")
        print(f"comb in: {energy_in:.3f} \n")
        return (energy_convected_in + energy_in) / air_out / self.cp_gas, air_out
        # divide by heat cap
    
    def add_to_plot(self, thrust, temperature, total_equivalence_ratio):
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
