import numpy as np

class MultiShomate:
    def __init__(self, shomate_coeff:np.ndarray, temp_range: np.ndarray):
        if shomate_coeff.shape[0] != temp_range.size - 1:
            raise ValueError(f"number of shomate coeffs and temperature range don't agree: {shomate_coeff.shape[0]} sets of coeffs provided but temp range only supports {temp_range.size - 1}.")
        self.coeffs = shomate_coeff
        self.temp_range = temp_range
    
    def enthalpy(self, temperature):
        # figure out which temp range we're in
        poly_num = int(np.where(temperature < self.temp_range)[0][0]) - 1
        c = self.coeffs[poly_num]
        # 0 1 2 3 4 5 6 7
        # a b c d e f g h
        return c[0] * temperature + c[1] * temperature**2 / 2 + c[2] * temperature ** 3 / 3 + c[3] * temperature**4 / 4 + c[5] - c[7]
    
    def cp(self, temperature):
        # figure out which temp range we're in
        poly_num = int(np.where(temperature < self.temp_range)[0])
        c = self.coeffs[poly_num]
        # 0 1 2 3 4 5 6 7
        # a b c d e f g h
        return c[0] + c[1] * temperature + c[2] * temperature ** 2 + c[3] * temperature ** 3 + c[4] / temperature ** 2
    

# obtained from https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1&Type=JANAFG&Table=on#JANAFG
CO2_ShoCoeff = np.array([[24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224],
                         [58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393.5224]])
CO2_TempRange = np.array([298, 1200, 6000])


co2 = MultiShomate(CO2_ShoCoeff, CO2_TempRange)

if __name__ == "__main__":
        print(co2.enthalpy(300)/1e3)
        print(co2.enthalpy(800)/1e3)
        print(co2.enthalpy(1199)/1e3)
        print(co2.enthalpy(1201)/1e3)
