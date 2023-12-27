import numpy as np
import matplotlib.pyplot as pt 

#4C12H23 + 81O2 -> 48CO2 + 46H2O
HfC12H23=-61.7*4184 #J/mole
HfO2=0 #J/mole
HfCO2=-393.5e3 #J/mole
HfH2O=-241.8e3 #J/mole

calvalue= 48*HfCO2 + 46*HfH2O - 4*HfC12H23 - 81*HfO2

#Values at different thrust levels
#20 kN, 40 kN, 80 kN, 120 kN 

mdot3=[10,16,29,39] #kg/s
P3=[700e3,1100e3,2200e3,3200e3] #Pa
T3=[540,610,725,820] #K
T4=[1150,1230,1410,1630] #K

cceff=0.995
cp_gas=1150

mdot=[]
for i in range(4):
    mdot.append((mdot3[i]*cp_gas*(T4[i]-T3[i]))/(cceff*abs(calvalue)))


#Equivalence Ratio Calcs
FAR_stoich_ker=0.0682 #READER
ER=[]

for i in range(4):
    FAR_actual=mdot[i]/mdot3[i]
    ER.append(FAR_actual/FAR_stoich_ker)

#Equivalence Ratio Calcs (In zones)
