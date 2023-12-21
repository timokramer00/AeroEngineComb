import numpy as np
import matplotlib.pyplot as pt 

#4C12H23 + 81O2 -> 48CO2 + 46H2O
HfC12H23=-61.7*4184 #J/mole
HfO2=0 #J/mole
HfCO2=-393.5e3 #J/mole
HfH2O=-241.8e3 #J/mole

calvalue= 48*HfCO2 + 46*HfH2O - 4*HfC12H23 - 81*HfO2