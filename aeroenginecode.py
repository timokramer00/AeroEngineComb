import numpy as np
import matplotlib.pyplot as pt 

#4C12H23 + 81O2 -> 48CO2 + 46H2O
HfC12H23=-61.7*4184 #J/mole
HfO2=0 #J/mole
HfCO2=-393.5e3 #J/mole
HfH2O=-241.8e3 #J/mole

calvalue= (12*HfCO2 + 46/4*HfH2O - HfC12H23)/0.1673114

#Values at different thrust levels
#20 kN, 40 kN, 80 kN, 120 kN 

mdot3=[10,16,29,39] #kg/s
P3=[700e3,1100e3,2200e3,3200e3] #Pa
T3=[540,610,725,820] #K
T4=[1150,1230,1410,1630] #K
perc_air_zones=[0.16+0.14,0.2+0.16]

cceff=0.995
cp_gas=1150

mdot=[]
for i in range(4):
    mdot.append((mdot3[i]*cp_gas*(T4[i]-T3[i]))/(cceff*abs(calvalue)))


#Equivalence Ratio Calcs
FAR_stoich_ker=0.0682 #READER
ER=[]
ER_PZ=[]
ER_SZ=[]
rho3=[]
R=287
A=2*np.pi*0.255*0.25
V=[]
res_time=[]


for i in range(4):
    FAR_actual=mdot[i]/mdot3[i]
    ER.append(FAR_actual/FAR_stoich_ker)
    rho3.append(P3[i]/(R*T3[i]))
    V.append(mdot3[i]/(rho3[i]*A))
    res_time.append(0.25/V[i])

mdot_air_PZ=[]
mdot_air_SZ=[]

#Equivalence Ratio Calcs (In zones)

for i in range(4):
    ER_PZ.append(ER[i]/perc_air_zones[0])
    ER_SZ.append(ER[i]/(perc_air_zones[0]+perc_air_zones[1]))

print(res_time)
print(ER)

