import numpy as np

class MultiShomate:
    def __init__(self, shomate_coeff:np.ndarray, temp_range: np.ndarray):
        if shomate_coeff.shape[0] != temp_range.size - 1:
            raise ValueError(f"number of shomate coeffs and temperature range don't agree: {shomate_coeff.shape[0]} sets of coeffs provided but temp range only supports {temp_range.size - 1}.")
        self.coeffs = shomate_coeff
        self.temp_range = temp_range
    
    def enthalpy(self, temperature):
        """
        Takes temperature in Kelvin, returns enthalpy in J/mol
        """
        # figure out which temp range we're in
        poly_num = int(np.where(temperature < self.temp_range)[0][-1]) - 1
        c = self.coeffs[poly_num]
        t = temperature/1000
        # 0 1 2 3 4 5 6 7
        # a b c d e f g h
        return (c[0] * t + c[1] * t**2 / 2 + c[2] * t ** 3 / 3 + c[3] * t**4 / 4 - c[4]/t + c[5] - c[7]) * 1e3
    
    def cp(self, temperature):
        """
        Takes temperature in Kelvin, returns Cp in J/mol/K
        """
        # figure out which temp range we're in
        poly_num = int(np.where(temperature < self.temp_range)[0][0]) - 1
        c = self.coeffs[poly_num]
        t = temperature/1000
        # 0 1 2 3 4 5 6 7
        # a b c d e f g h
        return c[0] + c[1] * t + c[2] * t ** 2 + c[3] * t ** 3 + c[4] / t ** 2
    
class NASAPoly:
    def __init__(self, polyCoeff:np.ndarray, temp_range:np.ndarray):
        if polyCoeff.size != 14:
            raise ValueError(f"PolyCoeff array is too small: 14 entries are needed, it has {polyCoeff.size}")
        if temp_range.size != 3:
            raise ValueError(f"temp_range array is too small: 3 entries are needed, it has {temp_range.size}")
        self.c_hi = polyCoeff[:7]
        self.c_lo = polyCoeff[7:]
        self.temp_range = temp_range
        self.R = 8.314  # universal gas constant
        
    def cp(self, temperature):
        """
        Takes temperature in Kelvin, returns Cp in J/mol/K
        """
        t = temperature
        c = self._get_poly(temperature)
        return self.R * (c[0] + c[1] * t + c[2] * t ** 2 + c[3] * t ** 3 + c[4] * t ** 4)
    
    def enthalpy(self, temperature):
        """
        returns enthalpy in J/mol/K
        """
        t = temperature
        c = self._get_poly(temperature)
        return self.R * t * (c[0] + c[1] * t / 2 + c[2] * t ** 2 / 3 + c[3] * t **3 /4 + c[4] * t **4 / 5 + c[5] / t)
        #pass
        
    def _get_poly(self, temp):         
        poly_num = int(np.where(temp < self.temp_range)[0][0]) - 1
        if poly_num == 0: # low range
            c = self.c_lo
        elif poly_num == 1: # high range
            c = self.c_hi
        else:
            print(f"que? {poly_num}")
        return c
        

# obtained from https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1&Type=JANAFG&Table=on#JANAFG
CO2_ShoCoeff = np.array([[24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431, -393.5224],[58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125, -393.5224]])
CO2_TempRange = np.array([298, 1200, 6000])

O2_ShoCoeff = np.array([[31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0], 
                        [30.03235, 8.772972, -3.988133, 0.788313, -0.741599, -11.32468, 236.1663, 0], 
                        [20.91111, 10.72071, -2.020498, 0.146449, 9.245722, 5.337651, 237.6185, 0]])
O2_TempRange = np.array([100, 700, 2000, 6000])

H2O_ShoCoeff = np.array([[30.092, 6.832514, 6.793435, -2.534480, 0.082139, -250.8810, 223.3967, -241.8264], [41.96426, 8.622053, -1.499780, 0.098119, -11.15764, -272.1797, 219.7809, -241.8264]])
H2O_TempRange = np.array([500, 1700, 6000])

N2_ShoCoeff = np.array([[19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900, 0], 
                        [35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 224.9810, 0]])
N2_TempRange = np.array([500, 2000, 6000])

# Poly from: http://combustion.berkeley.edu/gri_mech/data/nasa_plnm.html
# Values from: https://web.stanford.edu/group/haiwanglab/HyChem/approach/Report_Jet_Fuel_Thermochemical_Properties_v6.pdf
fuel_NASACoeff = np.array([0.20115053E+02, 0.69181815E-01, -0.23171997E-04, 0.19519781E-08, 0.27041864E-12, -0.42549684E+05, -0.79951393E+02, 0.40264606E-01, 0.88749342E-01, 0.42724278E-04,
-0.11716354E-06, 0.53996904E-10, -0.35190535E+05, 0.32478203E+02])
fuel_TempRange = np.array([298, 1000, 9000])

CO2 = MultiShomate(CO2_ShoCoeff, CO2_TempRange)
O2 = MultiShomate(O2_ShoCoeff, O2_TempRange)
water = MultiShomate(H2O_ShoCoeff, H2O_TempRange)
N2 = MultiShomate(N2_ShoCoeff, N2_TempRange)
fuel = NASAPoly(fuel_NASACoeff, fuel_TempRange)

if __name__ == "__main__":
        #print(CO2.enthalpy(300))
        #print(CO2.enthalpy(800))
        #print(CO2.enthalpy(1199))
        #print(CO2.enthalpy(1201))
        #print(CO2.enthalpy(1200))
        
        
        molar_mass = {  # kg/mol
            "CO2": (12.011 + 15.99*2) * 1e-3,
            "H2O": (1.0078*2 + 15.99) * 1e-3,
            "fuel": (12*12.011 + 23 * 1.0078) * 1e-3,  # c12h23
            "O2": (15.99*2) * 1e-3,
            "N2": 28.0134*1e-3
        }
        
        import matplotlib.pyplot as plt
        
        def vectorize(func, arg):
            out = []
            for a in arg:
                out.append(func(a)/molar_mass["fuel"])
            return out
        
        
        t1 = np.linspace(300, 2999, 500)
        t2 = np.linspace(1000, 2999, 500) 
        
        cp1 = vectorize(fuel.enthalpy, t1)
        cp2 = vectorize(fuel.enthalpy, t2)
        plt.plot(t1, cp1)
        plt.plot(t2, cp2)
        plt.grid()
        plt.show()
        
