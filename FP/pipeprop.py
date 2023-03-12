import timeit
import math
import numpy as np
from Korrelation.trend_func import trendprop
from scipy.optimize import root
from scipy.optimize import minimize

def Trend_p_min(fluid):
    return 1

def Trend_TripSI(fluid):               
        triple = trendprop(fluid).triple()
        T_trip = triple[0]
        p_trip = triple[1]
        return T_trip, p_trip
     
def check_props(props, prop_ls = None):
    prop_ls_default = [
        'p' ,
        'T',
        'Q' ,
        'cp' ,
        'cv' ,
        'd' ,
        'h' ,
        's',
        'vs',
        'L' ,
        'betas', 
        'ks' ,
        'Ds/Dp|d' ,
        'Ds/Dp|E' ,
        'Dp/DT|d' ,
        'Dp/Dd|T' ,
        'Dh/Dp|T' ,
        'Dh/Dp|E' ,
        'DT/Dp|E' ,
    ]
    if prop_ls is None:
        prop_ls = prop_ls_default
    else:
        if prop_ls < prop_ls_default:
            pass
    res = True
    if isinstance(props, dict):        
        for key, val in props.items():
            if key not in prop_ls:
                print('invalid property: {0} = {1}'.format(key, val))
                res = False        
    else:
        res = False
    return res

def Trend_props_E(fluid, prop1, prop2, in1 = None, in2 = None, 
                 inkl = [], exkl = []):
    """
    equilibrium phase property
    ink = list()
    exk = list()
    unit system K, Pa, kg
    """
    
    prop_codes = {
        'p' : 'P',          # Pa
        'T' : 'T',          # K
        'Q' : 'Q',          # kg/kg
        'cp' : 'CP',        # J/kg/K
        'cv' : 'CV',        # J/kg/K
        'betas' : 'betas',
        'd' : 'D',          # kg/m3
        'h' : 'H',          # J/kg
        's' : 'S',          # J/kgK
        'vs': 'VS',         # m2/s
        'L' : 'L',          # W/mK
        'ks' : 'ks',
        'Dp/DT|d' : 'Dp/DT|d',
        'Dp/Dd|T' : 'Dp/Dd|T',
    }
    
    props = {}
    if inkl:
        prop_ls = inkl
    else:
        prop_ls = list(set(prop_codes.keys() - set(exkl)))
    
    for prop in prop_ls:
        
        try:
            props.update(
                {
                    prop : trendprop(
                            fluid, 
                            input_code = in1 + in2, 
                            prop1 = prop1, prop2 = prop2, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                prop_codes[prop], 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                            )
                }
            )
        except ValueError as e:
            if 'iteration' in str(e):
                if in2 == 'Q' and (in1 == 'T' or in1 == 'p'):
                    h0 = trendprop(
                            fluid, 
                            input_code = in1 + 'liq', 
                            prop1 = prop1, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    h1 = trendprop(
                            fluid, 
                            input_code = in1 + 'vap', 
                            prop1 = prop1, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                            )
                    h = h0 + (h1 - h0) * prop2
                    props.update(
                        {
                            prop : trendprop(
                                    fluid, 
                                    input_code = in1 + 'h', 
                                    prop1 = prop1, prop2 = h, 
                                    unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                    err_capt = True
                                ).Prop(
                                        prop_codes[prop], 
                                        unit_T = 'K', 
                                        unit_p = 'Pa', 
                                        unit_m = 'kg'
                                    )
                        }
                    )
                elif in1 == 'Q' and (in2 == 'T' or in2 == 'p'):
                    h0 = trendprop(
                            fluid, 
                            input_code = in2 + 'liq', 
                            prop1 = prop2, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    h1 = trendprop(
                            fluid, 
                            input_code = in2 + 'vap', 
                            prop1 = prop2, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                            )
                    h = h0 + (h1 - h0) * prop1
                    props.update(
                        {
                            prop : trendprop(
                                    fluid, 
                                    input_code = in2 + 'h', 
                                    prop1 = prop2, prop2 = h, 
                                    unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                    err_capt = True
                                ).Prop(
                                        prop_codes[prop], 
                                        unit_T = 'K', 
                                        unit_p = 'Pa', 
                                        unit_m = 'kg'
                                    )
                        }
                    )

                else:
                    raise ValueError(e)
            elif str(e) == 'Call of Flash_pure with press <= ptp or press >= pc '\
            or str(e) == 'Call of Flash_pure with Temp <= ttp or Temp >= tc '\
            or str(e) == 'Subroutine PHASEDET_PURE: Flash calculation did not converge '\
            or str(e) == 'Temperature <= Tmin of transport equations '\
            or str(e) == 'Temperature <= Tminfluid ':
                if prop == 'Q' and (in1 == 'Q' or in2 == 'Q'):
                    if in1 == 'Q':
                        Q = prop1
                    else:
                        Q = prop2
                    props.update({prop : Q})
                elif prop == 'Q' and (in1 == 'T' or in1 == 'p'):
                    h = trendprop(
                            fluid, 
                            input_code = in1 + in2 + '+', 
                            prop1 = prop1, prop2 = prop2, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    h0 = trendprop(
                            fluid, 
                            input_code = in1 + 'subs+', 
                            prop1 = prop1, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    h1 = trendprop(
                            fluid, 
                            input_code = in1 + 'subv+', 
                            prop1 = prop1, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                            )
                    x = (h - h0) / (h1 - h0)
                    props.update({prop : x})
                elif prop == 'Q' and (in2 == 'T' or in2 == 'p'):
                    h = trendprop(
                            fluid, 
                            input_code = in1 + in2 + '+', 
                            prop1 = prop1, prop2 = prop2, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    h0 = trendprop(
                            fluid, 
                            input_code = in2 + 'subs+', 
                            prop1 = prop1, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    h1 = trendprop(
                            fluid, 
                            input_code = in2 + 'subv+', 
                            prop1 = prop2, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                            )
                    x = (h - h0) / (h1 - h0) 
                    props.update({prop : x})
                elif in2 == 'Q' and (in1 == 'T' or in1 == 'p'):
                    h0 = trendprop(
                            fluid, 
                            input_code = in1 + 'subs+', 
                            prop1 = prop1, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    h1 = trendprop(
                            fluid, 
                            input_code = in1 + 'subv+', 
                            prop1 = prop1, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                            )
                    h = h0 + (h1 - h0) * prop2
                    p = trendprop(
                            fluid, 
                            input_code = in1 + 'subs+', 
                            prop1 = prop1, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'p', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    props.update(
                        {
                            prop : trendprop(
                                    fluid, 
                                    input_code = 'ph+', 
                                    prop1 = p, prop2 = h, 
                                    unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                    err_capt = True
                                ).Prop(
                                        prop_codes[prop], 
                                        unit_T = 'K', 
                                        unit_p = 'Pa', 
                                        unit_m = 'kg'
                                    )
                        }
                    )
                elif in1 == 'Q' and (in2 == 'T' or in2 == 'p'):
                    h0 = trendprop(
                            fluid, 
                            input_code = in2 + 'subs+', 
                            prop1 = prop2, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
                    h1 = trendprop(
                            fluid, 
                            input_code = in2 + 'subv+', 
                            prop1 = prop2, prop2 = 0, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'h', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                            )
                    h = h0 + (h1 - h0) * prop1
                    props.update(
                        {
                            prop : trendprop(
                                    fluid, 
                                    input_code = in2 + 'h+', 
                                    prop1 = prop2, prop2 = h, 
                                    unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                    err_capt = True
                                ).Prop(
                                        prop_codes[prop], 
                                        unit_T = 'K', 
                                        unit_p = 'Pa', 
                                        unit_m = 'kg'
                                    )
                        }
                    )
                else:

                    props.update(
                        {
                            prop : trendprop(
                                    fluid, 
                                    input_code = in1 + in2 + '+', 
                                    prop1 = prop1, prop2 = prop2, 
                                    unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                    err_capt = True
                                ).Prop(
                                        prop_codes[prop], 
                                        unit_T = 'K', 
                                        unit_p = 'Pa', 
                                        unit_m = 'kg'
                                    )
                        }
                    )
            elif str(e) == 'Property in homogeneous region undefined ':
                props.update({prop : 'NA'})
            else:
                print(in1, in2)
                print(prop1, prop2)
                raise ValueError(e)

    return props

"""
        l:solid length 
        R:rohr radius
        H:rohr length
        A:area
        U:velocity
        Q:transferred heat flow
        G:gas area ratio
        
        g:gas
        s:solid
        w:wall
        
        1:inout
        2:output
        m:average
        
"""
# Gib die Parameter der Rohrelemente
def Rohrelemente(R, H, T_w, n):
    rohr = {
        'R' : R,
        'H' : H,
        'T_w' : T_w,
        'n' : n
        }
    return rohr


# definiert die Festphase 
class Solid():
    
    # Initialisieren die Eigenschaften der Feststoffphase 
    def __init__(self, fluid, m_s1, p_1, G_1, R, H):
        
        # Initialisieren die verwendete Fluid
        self.fluid = fluid
        # Initialisieren den Eingangmassestrome von Gas-und Festphase  
        self.m_s1 = m_s1
        # Initialisieren den Eingangsdruck
        self.p_1 = p_1
        # Initialisieren relative Flaeche
        self.G_1 = G_1
        # Initialisieren den Radius und Rohrlaenge
        self.R = R
        self.H = H
        # Mit Trend-Funktion andere Eigenschaften berechnen 
        self.h_s1 = self.prop1()['h']
        self.T_s1 = self.prop1()['T']
        self.d_s1 = self.prop1()['d']
        # definieren die gesamte Querschnittsflaeche 
        self.A = math.pi * R ** 2
        
    # Mit Trend-Funktion die Eingangseigenschaften der Festphase berechnen 
    def prop1(self):
        # Mit Eingangsdruck die andere Eigenschaften auf der Sublimationskurve berechnen 
        props_1 = Trend_props_E(self.fluid, self.p_1, 1e-8,
                                 in1 = 'p', in2 = 'Q',
                                 inkl = ['T', 'd', 'h']
                                 )
        return props_1
    
    # Mit Trend-Funktion die Eigenschaften der festphase am Ausgang 2 berechnen
    def prop2(self, p_2):
        # Mit Ausgangdruck die andere Eigenschaften auf der Sublimationskurve berechnen 
        props_2 = Trend_props_E(self.fluid, p_2, 1e-8,
                                 in1 = 'p', in2 = 'Q',
                                 inkl = ['T','d','h']
                                 )   
        return props_2
    
    # Mit Trend-funktion die durchschnittlichen Eigenschaften der Festphase berechnen
    def propm(self, p_m):
        # Mit durchschnittliche Druck die andere Eigenschaften auf der Sublimationskurve berechnen 
        props_m = Trend_props_E(self.fluid, p_m, 1e-8,
                                 in1 = 'p', in2 = 'Q',
                                 inkl = ['T','d','cp']
                                 )
        return props_m
    
    # Gib die durchschnittliche relative Flaeche der Gasphase
    def Gm(self, G_2):
        return 0.5 * (self.G_1 + G_2)
    
    # Gib die mittele Massestrom der Festphase  
    def solidmass(self, m_s2):
        return 0.5 * (self.m_s1 + m_s2)
    
    # Gib die mittlere Feststofftemperatur 
    def solidtemperatur(self, T_s2):
        return 0.5 * (self.T_s1 + T_s2)
    
    # Gib die mittele Solidgeschwindigkeit
    def solidvelocity(self, U_s1, U_s2):
        return 0.5 * (U_s1 + U_s2)
      
    # Gib die Flaesche von Feststoff-Stroemung 
    def solidarea(self, G):
        return (1 - G) * self.A
    
    # Gib die Flaesche von Feststoff-Wand-surface 
    def solidwallarea(self, G):
        # den Winkel berechnen
        res = root(fun = lambda x: (x - np.sin(x) - 2 * math.pi * G),
                   x0 = 0.001)
        a = res.x[0]
        
        A_ws = (2 * math.pi - a) * self.R * self.H
        return A_ws
        
    # Gib die von Wand zur Solidphase uebertragener Waermestrom
    def Wallsolidheattransfer(self, k, T_w, T_s, A_ws):
        Q_ws = k * A_ws * (T_w - T_s)
        return Q_ws
    
    
# definiert die Gasphase   
class Gas():
    
    # Initialisieren die Eigenschaften der Gasphase 
    def __init__(self, fluid, m_g1, h_g1, p_1, G_1, R, H):
        # Initialisieren die verwendete Fluid
        self.fluid = fluid
        # Initialisieren den Eingangmassestrome von Gas-und Festphase  
        self.m_g1 = m_g1
        self.h_g1 = h_g1
        # Initialisieren den Eingangsdruck
        self.p_1 = p_1
        # Initialisieren relative Flaeche
        self.G_1 = G_1
        # Initialisieren den Radius und Rohrlaenge
        self.R = R
        self.H = H
        # Mit Trend-Funktion andere Eigenschaften berechnen
        self.T_g1 = self.prop1()['T']
        self.d_g1 = self.prop1()['d']
        # definieren die gesamte Querschnittsflaeche 
        self.A = math.pi * R ** 2    
    
    # Mit Trend-Funktion die Eigenschaften der Gasphase am Eingang 1 berechnen
    def prop1(self):
        # Mit Eingangsdruck- und Enthalpie die andere Eigenschaften berechnen
        props_1 = Trend_props_E(self.fluid, self.p_1, self.h_g1,
                                 in1 = 'p', in2 = 'h',
                                 inkl = ['T', 'd',]
                                 )
        return props_1 
    
    # Mit Trend-Funktion die Eigenschaften der Gasphase am Ausgang 2 berechnen
    def prop2(self, p_2, h_g2):
        # Mit Ausgangsdruck- und Enthalpie die andere Eigenschaften berechnen
        props_2 = Trend_props_E(self.fluid, p_2, h_g2,
                                 in1 = 'p', in2 = 'h',
                                 inkl = ['T', 'd']
                                 )    
        return props_2
    
    # Mit Trend-Funktion die durchschnittliche Eigenschaften der Gasphase berechnen
    def propm(self, p_m, T_m):
        # Mit durchschnittliche Druck- und Temperatur die andere Eigenschaften berechnen
        props_m = Trend_props_E(self.fluid, p_m, T_m,
                                 in1 = 'p', in2 = 'T',
                                 inkl = ['h', 'd', 'L', 'cp', 'vs'])
        return props_m
    
    # Gib die mittele Gasdruck
    def pressurem(self, p_2):
        return 0.5 * (self.p_1 + p_2)
    
    # Gib die mittele Gasdruck
    def gasenthalpiem(self, h_g2):
        return 0.5 * (self.h_g1 + h_g2)
    
    # Gib die mittele Gasmassestrom
    def gasmass(self, m_g2):
        return 0.5 * (self.m_g1 + m_g2)
    
    # Gib die mittele Gastemperatur 
    def gastemperatur(self, T_g2):
        return 0.5 * (self.T_g1 + T_g2)
    
    # Gib die mittele Gasgeschwindigkeit
    def gasvelocity(self, U_g1, U_g2):
        return 0.5 * (U_g1 + U_g2)
    
    # Gib die Flaesche von Gas-Stroemung 
    def gasarea(self, G):
        return G * self.A
    
    # Aequivalente Laenge der Gasphase
    def equradius(self,G_m):
        # den Winkel berechnen
        res = root(fun = lambda x: (x - np.sin(x) - 2 * math.pi * G_m),
                   x0 = 0.001)
        a = res.x
        # Gas-Wand-Oberflaeche berechnen  
        A_g = self.gasarea(G_m)
        # Gas-Wand-Umfang berechnen
        u_g = a * self.R
        # aequivalente Laenge berechnen
        d_h = 4 * A_g / u_g
        return d_h
    
    # Gib die Flaesche von Gas-Wand-surface 
    def gaswallarea(self, G):
        # den Winkel berechnen
        res = root(fun = lambda x: (x - np.sin(x) - 2 * math.pi * G),
                   x0 = 0.001)
        a = res.x[0]
        
        A_wg = a * self.R * self.H
        return A_wg
    
    # Gib die Flaesche von Gas-Feststoff-Querschnitt
    def gassolidarea(self, G):
        # den Winkel berechnen
        res = root(fun = lambda x: (x - np.sin(x) - 2 * math.pi * G),
                   x0 = 0.001)
        a = res.x[0]
        
        A_Q = 2 * self.R * np.sin(0.5 * a) * self.H
        return A_Q    
    
    def LmTd(self, T_in1, T_in2, T_out1, T_out2):
        return ((T_out1 - T_out2) - (T_in1 - T_in2)) / (math.log((T_out1 - T_out2) / (T_in1 - T_in2)))
    
    # Gib den von Wand zur Gasphase uebertragenen Waermestrom
    def Wallgasheattransfer(self, alpha, T_w, T_g2, A_wg):
        # DTm = self.LmTd(T_w, self.T_g1, T_w, T_g2)
        Q_wg = alpha * A_wg * (T_w - T_g2)
        return Q_wg
    
    # Gib den von Gasphase zur Festphase uebertragenen Waermestrom
    def Gassolidheattransfer(self, k, T_s, T_g2, A_Q):
        # DTm = self.LmTd(self.T_g1, T_s, T_g2, T_s)
        Q_gs = k * A_Q * (T_g2 - T_s)
        return Q_gs

if __name__ == '__main__':
    start = timeit.default_timer()

    h_subs = trendprop('co2', input_code = 'PSUBS+', prop1 = 1.01325e5, prop2 = 1, err_capt = True, 
                       unit_p = 'Pa', unit_T = 'C').Prop('H', unit_T = 'C',unit_m = 'kg')
    h_subv = trendprop('co2', input_code = 'PSUBV+', prop1 = 1.01325e5, prop2 = 1, err_capt = True, 
                       unit_p = 'Pa', unit_T = 'C').Prop('H', unit_T = 'C',unit_m = 'kg') 
    
    
    prop_gas = Trend_props_E(fluid='co2', prop1=1.01e5, prop2=5.14e5,
                                    in1 = 'p', in2 = 'h',
                                    inkl = ['T','h','d','s','cp','vs','L'])
    
    prop_solid = Trend_props_E(fluid='co2', prop1=1.01325e5, prop2=-2e5,
                                    in1 = 'p', in2 = 'h',
                                    inkl = ['T','h','d','s',])

    stop = timeit.default_timer()
    print('Time: ', stop - start)     
    
def Trend_props_E_ext(fluid, prop1, prop2, in1 = None, in2 = None, 
                 inkl = [], exkl = []):

    prop_codes = {
        'p' : 'P',      # Pa
        'T' : 'T',      # K
        'Q' : 'Q',      # kg/kg  Dampfqualität
        'cp' : 'CP',    # J/kg/K
        'cv' : 'CV',    # J/kg/K
        'betas' : 'betas',
        'd' : 'D',      # kg/m3
        'h' : 'H',      # J/kg
        's' : 'S',      # J/kgK
        'vs': 'VS',     # m2/s
        'L' : 'L',      # W/m/K
        'ks' : 'ks',
        'Ds/Dp|E' : '',
        'Dp/DT|d' : 'Dp/DT|d',
        'Dp/Dd|T' : 'Dp/Dd|T',
        'Dh/Dp|T' : '',
        'Dh/Dp|E' : '',
        'Ds/Dp|E' : '',
        'DT/Dp|E' : '',
    }
    
    props = {}
    post = False
    prop_ls_post = []
    exkl_post = []
    if inkl:
        prop_ls = inkl
    else:
        prop_ls = list(set(prop_codes.keys()  - set(exkl)))

    for prop in prop_ls:
        if prop == 'Dh/Dp|T':
            prop_ls1 = ['d', 'T', 'Dp/DT|d', 'Dp/Dd|T']
            for prop0 in prop_ls1:
                if prop0 not in prop_ls:
                    prop_ls.append(prop0)
                    exkl_post.append(prop0)
                post = True
            prop_ls_post.append(prop)
            exkl.append(prop)

        if prop == 'Ds/Dp|E':
            if in1 == 'p' or in1 == 'T':
                props.update(
                    {
                        'Ds/Dp|E' : ds_dp_lE(fluid, prop1, in1 = in1)
                    }
                )
            elif in2 == 'p' or in2 == 'T':
                props.update(
                    {
                        'Ds/Dp|E' : ds_dp_lE(fluid, prop2, in1 = in2)
                    }
                )
            exkl.append(prop) 
        if prop == 'Dh/Dp|E':
            props.update({'Dh/Dp|E' : dh_dp_E(fluid, prop1, prop2, 
                                              in1 = in1, in2 = in2)})   
            exkl.append(prop)  
        if prop == 'DT/Dp|E':
            props.update({'DT/Dp|E' : dT_dp_E(fluid, prop1, prop2, 
                                              in1 = in1, in2 = in2)})               
            exkl.append(prop) 
  
    props.update(Trend_props_E(fluid, prop1, prop2, in1 = in1, in2 = in2, 
                     exkl = list(set(prop_codes.keys()) - set(prop_ls))))
     
    
    if post is True:
        for prop in prop_ls_post:
            if prop == 'Dh/Dp|T':
                res = (- props['d'] * props['Dp/Dd|T'] + props['T'] * props['Dp/DT|d']) / ( - props['d'] ** 2 * props['Dp/Dd|T'])
                props.update({prop : res})
            
        for prop in exkl_post:
            del props[prop]
                
    return props 

def ds_dp_lE(fluid, prop1, in1 = 'p'):
    sp0 = trendprop(
                            fluid, 
                            input_code = in1 + 'liq', 
                            prop1 = prop1, prop2 = 1, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        )
    try:
        s0 = sp0.Prop('s', unit_T = 'K', unit_p = 'Pa', unit_m = 'kg')
        sp1 = trendprop(
                    fluid, 
                    input_code = in1 + 'liq', 
                    prop1 = prop1 * 1.0001, prop2 = 1, 
                    unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                    err_capt = True
                )
    except ValueError as e:
        if str(e) == 'Call of Flash_pure with press <= ptp or press >= pc ':
            sp0 = trendprop(
                    fluid, 
                    input_code = in1 + 'subs+', 
                    prop1 = prop1, prop2 = 1, 
                    unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                    err_capt = False
                )
            s0 = sp0.Prop('s', unit_T = 'K', unit_p = 'Pa', unit_m = 'kg')
            sp1 = trendprop(
                fluid, 
                input_code = in1 + 'subs+', 
                prop1 = prop1 * 0.9999, prop2 = 1, 
                unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                err_capt = False
            )
            
    p0 = sp0.Prop('p', unit_T = 'K', unit_p = 'Pa', unit_m = 'kg')

    s1 = sp1.Prop('s', unit_T = 'K', unit_p = 'Pa', unit_m = 'kg')
    p1 = sp1.Prop('p', unit_T = 'K', unit_p = 'Pa', unit_m = 'kg')
    return (s1 - s0) / (p1 - p0)     

def dh_dp_E(fluid, prop1, prop2, in1 = 'p', in2 = 'Q', prop_E = Trend_props_E):
    sp0 = prop_E(fluid, prop1, prop2, in1 = in1, in2 = in2, inkl = ['p', 'h'])
    sp1 = prop_E(fluid, prop1 * 0.9999, prop2, in1 = in1, in2 = in2, inkl = ['p', 'h'])        
    return (sp1['h'] - sp0['h']) / (sp1['p'] - sp0['p'])     

def dT_dp_E(fluid, prop1, prop2, in1 = 'p', in2 = 'Q'):
    sp0 = Trend_props_E(fluid, prop1, prop2, in1 = in1, in2 = in2, inkl = ['p', 'T'])
    sp1 = Trend_props_E(fluid, prop1 * 0.9999, prop2, in1 = in1, in2 = in2, inkl = ['p', 'T'])
    print(sp0, sp1)
    if sp1['T'] == 0:
        sp1 = Trend_props_E(fluid, prop1 * 1.0001, prop2, in1 = in1, in2 = in2, inkl = ['p', 'T']) 
    return (sp1['T'] - sp0['T']) / (sp1['p'] - sp0['p'])  

def props(fluid, p, s):
    props_E = Trend_props_E(fluid, p, s, in1 = 'p', in2 = 's', 
                 inkl = ['cp', 'cv' , 'd'])

    return props_E

def props_lMeta(fluid, p, s):
    props_E = Trend_props_lMeta(fluid, p, s, in1 = 'p', in2 = 's', 
                 inkl = ['cp', 'cv' , 'd'])

    return props_E

def integ_vdp_s(fluid, p1, p2, s1, step = 100):    
    dp = (p2 - p1) / step    
    p = np.arange(p1, p2, dp)
    
    v = np.zeros(step-1)
    gam = np.zeros(step-1)
    vdp = np.zeros(step-1)
    for i in np.arange(step - 1):
        p_m = (p[i] + p[i + 1]) / 2
        gam[i] = props(fluid, p_m, s1)['cp'] / props(fluid, p_m, s1)['cv']
        v[i] = 1 / props(fluid, p_m, s1)['d']
                
        #pdV[i] = gam[i] /  (gam[i] - 1) * (p[i+1] * v[i+1] - p[i] * v[i])
        vdp[i] = v[i] * dp

    res = sum(vdp)
    return res    

def integ_vdp_s_lMeta(fluid, p1, p2, s1, step = 100):    
    dp = (p2 - p1) / step    
    p = np.arange(p1, p2, dp)
    
    v = np.zeros(step-1)
    gam = np.zeros(step-1)
    vdp = np.zeros(step-1)
    for i in np.arange(step - 1):
        p_m = (p[i] + p[i + 1]) / 2
        gam[i] = props_lMeta(fluid, p_m, s1)['cp'] / props_lMeta(fluid, p_m, s1)['cv']
        v[i] = 1 / props_lMeta(fluid, p_m, s1)['d']
        
        #pdV[i] = gam[i] / (gam[i] - 1) * (p[i+1] * v[i+1] - p[i] * v[i])
        vdp[i] = v[i] * dp
    res = sum(vdp)
    return res

def Trend_props_lMeta(fluid, prop1, prop2, in1 = 'p', in2 = 's', 
                 inkl = [], exkl = []):    
    """
    in1 = 'p' or 'T'
    in2
    """
    prop_codes = {
        'p' : 'P', # Pa
        'T' : 'T', # K
        'cp' : 'CP',    # J/kg/K
        'cv' : 'CV',    # J/kg/K
        'd' : 'D',     # kg/m3
        'h' : 'H', # J/kg
        's': 'S',     # J/kgK
    }
    
    props = {}
    
    if inkl:
        prop_ls = inkl
    else:
        prop_ls = list(set(prop_codes.keys() - set([in1, in2]) - set(exkl)))
    
    if (in1  == 'T' and in2 == 'd') \
    or (in1  == 'p' and in2 == 'T') \
    or (in1 == 'T' and in2 == 'p'):
        prophom = trendprop(
                            fluid, 
                            input_code = in1 + in2 +'&' , 
                            prop1 = prop1, prop2 = prop2, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        )
        prophom.allprop_hom()
        for prop in prop_ls:
            try: 
                props.update(
                            {
                                prop : prophom.Prop(
                                            prop_codes[prop], 
                                            unit_T = 'K', 
                                            unit_p = 'Pa', 
                                            unit_m = 'kg'
                                        )
                            }
                        )
            except ValueError as e:
                print(e)
    elif in1 == 'T':
        # calculate saturated liquid density 
        d0 = trendprop(
                            fluid, 
                            input_code = in1 + 'liq' , 
                            prop1 = prop1, prop2 = prop2, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'D', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
        def func_prop2_zero(d):
            prophom = trendprop(
                                fluid, 
                                input_code = 'TD&', 
                                prop1 = prop1, prop2 = d, 
                                unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                err_capt = False
                            )
            return abs( prop2 - prophom.Prop_hom(in2, 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg')
                )
            
        res = minimize(func_prop2_zero, 
                       d0,     
                       method = 'Nelder-Mead'
        )
        d = res.x[0]
        #print(res)
        prophom = trendprop(
                                fluid, 
                                input_code = 'TD&', 
                                prop1 = prop1, prop2 = d, 
                                unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                err_capt = True
                            )
        for prop in prop_ls:
            try: 
                props.update(
                            {
                                prop : prophom.Prop(
                                            prop_codes[prop], 
                                            unit_T = 'K', 
                                            unit_p = 'Pa', 
                                            unit_m = 'kg'
                                        )
                            }
                        )
            except ValueError as e:
                print(e)
                
    elif in1 == 'p':
        # calculate saturated liquid temperature
        d0 = trendprop(
                            fluid, 
                            input_code = in1 + 'liq' , 
                            prop1 = prop1, prop2 = prop2, 
                            unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                            err_capt = True
                        ).Prop(
                                'D', 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg'
                                )
        def func_prop2_zero(T):
            prophom = trendprop(
                                fluid, 
                                input_code = 'PT&', 
                                prop1 = prop1, prop2 = T, 
                                unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                err_capt = False
                            )
            return abs( prop2 - prophom.Prop_hom(in2, 
                                unit_T = 'K', 
                                unit_p = 'Pa', 
                                unit_m = 'kg')
                )
            
        res = minimize(func_prop2_zero, 
                       d0,     
                       method = 'Nelder-Mead'
        )
        T = res.x[0]
        #print(res)
        prophom = trendprop(
                                fluid, 
                                input_code = 'PT&', 
                                prop1 = prop1, prop2 = T, 
                                unit_T = 'K', unit_p = 'Pa', unit_m = 'kg',
                                err_capt = True
                            )
        for prop in prop_ls:
            try: 
                props.update(
                            {
                                prop : prophom.Prop(
                                            prop_codes[prop], 
                                            unit_T = 'K', 
                                            unit_p = 'Pa', 
                                            unit_m = 'kg'
                                        )
                            }
                        )
            except ValueError as e:
                print(e)
    return props
