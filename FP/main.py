import timeit
import pandas as pd 

from pipeprop import Rohrelemente as Rohr
from controllvolume import Prop_in
from controllvolume import controllvolume as cv

# Gib die Parameter der Rohrelemente
rohr = Rohr(R = 1.5e-3, H = 2, T_w = 205, n = 1000)

# Gib die Eintrittseigenschaften 
prop_0 = Prop_in(fluid='CO2', m_g_in=0.05e-3, m_s_in=0.05e-3, p_in=1.01325e5, 
                 h_g_in=4.23e5,# trendprop('co2', input_code = 'PSUBV+', prop1 = 3e5,prop2 = 1, err_capt = True, unit_p = 'Pa', unit_T = 'C').Prop('H', unit_T = 'C',unit_m = 'kg'),
                 G_in=0.9)  

# initialisieren den sublimierte Massestrom und controllvolumen sowie Propertyliste
# m_sub = 0
propslist = [prop_0]
cvlist = []
df = pd.DataFrame(propslist[0], index=[0]).T 

# Starten den Timer
start = timeit.default_timer()

# zählen 0-9
for i in range(rohr['n']):  
    # initialisieren das Kontrollvolumen 
    cvlist.append(cv(propslist[i], rohr)) 
    # Simulationsvorgan des Kontrollvolumens
    cvlist[i].solve_out()
    # die Ausgangsdaten
    propslist.append(cvlist[i].prop_out())
    # m_sub += propslist[i+1]['m_sub'] 
    # DataFrame + concat 
    df_i = pd.DataFrame(propslist[i+1], index=[i+1]).T
    df = pd.concat([df, df_i], axis=1)

print('Time: ', timeit.default_timer() - start)
