import timeit

from scipy.optimize import root
from scipy.optimize import fsolve
from scipy.optimize import minimize

from pipeprop import Trend_TripSI as Trip
from pipeprop import Solid,Gas
from pipeprop import Rohrelemente

from Korrelation.trend_func import trendprop
from Korrelation.velocity import velocity2 as vl
from Korrelation.pressuredrop import fg
from Korrelation.pressuredrop import fp
from Korrelation.force import wallgasfrictionforce as wgff
from Korrelation.force import wallsolidfrictionforce as wsff
from Korrelation.force import shearforce as sf
from Korrelation.purenumber import Nusseltnumber as nu
from Korrelation.purenumber import Prandtlnumber as pr
from Korrelation.purenumber import Renouldsnumber as re
from Korrelation.heattransfer import alpha_g as alg
from Korrelation.heattransfer import finalform2 as als
from Korrelation.heattransfer import thermaltransmittance as ttm

"""     
    Einheit:                          
        'm' : kg/s             
        'U' : m/s 
        'p' : Pa
        'T' : K
        'Q' : J/s = W
        'cp' : J/kg/K
        'cv' : J/kg/K
        'betas' :
        'd' : kg/m3
        'h' : J/kg
        's' : J/kgK
        'vs': m2/s
        'L' : W/mK
        alpha : W/m2/K
        k : W/m2/K
"""

# Gib die Eintrittseigenschaften 
def Prop_in(fluid, m_g_in, m_s_in, p_in, G_in, h_g_in):
    prop_in = {
            'fluid' : fluid,
            'G' : G_in,
            'T_g' : 0,
            'T_s' : 0,
            'h_g' : h_g_in, 
            'p' : p_in,
            'm_g' : m_g_in,
            'm_s' : m_s_in,
            'm_sub' : 0,
            'Q_gs' : 0,
            'Q_wg' : 0,
            'Q_ws' : 0,
            'Q_sub' : 0,
            'F_Q' : 0,
            'FR_wg': 0,
            'FR_ws': 0,
            }
    return prop_in


# definieren die Kontrollvolumen  
class controllvolume():
    
    # initializieren die Klasse Kontrollvolumen
    def __init__(self, prop_1, rohr, e_s = 0.01):
        # Initialisieren die verwendete Fluid
        self.fluid = prop_1['fluid']
        # Initialisieren die Massestrome von Gas-und Festphase am Eingang 1 
        self.m_g1, self.m_s1, self.m_sub_1 = prop_1['m_g'], prop_1['m_s'], prop_1['m_sub']
        # Initialisieren den Druck
        self.p_1 = prop_1['p']
        # Initialisieren die relative Flaeche der Gasphase 
        self.G_1 = prop_1['G']
        # Initialisieren die Enthalpie der Gasphase 
        self.h_g1 = prop_1['h_g']
        # Initialisieren den Radius des Rohrs 
        self.R = rohr['R']
        # Initialisieren die laenge der Kontrollvolumen 
        self.H = rohr['H']/rohr['n']
        # Initialisieren die Wandtemperatur 
        self.T_w = rohr['T_w']
        # Initialisieren die Porositaet
        self.e_s = e_s
        # erstellen ein Objekt für die Klasse Gasphase erstellen
        self.gas = Gas(self.fluid, self.m_g1, self.h_g1, self.p_1, self.G_1, self.R, self.H)
        # erstellen ein Objekt für die Klasse Festphase erstellen
        self.solid = Solid(self.fluid, self.m_s1, self.p_1, self.G_1, self.R, self.H)
        
        # Gib die Temperatur und Dichte der Gas- und Festphase am Eingang 1 
        self.T_g1, self.T_s1 = self.gas.T_g1, self.solid.T_s1
        # trendprop('co2', input_code = 'PSUBS+', prop1 = self.p_1, prop2 = 1, err_capt = True, unit_p = 'Pa', unit_T = 'C').Prop('T', unit_T = 'K',unit_m = 'kg')
        self.d_g1, self.d_s1 = self.gas.d_g1, self.solid.d_s1
        # Gib die Querschnittsflaeche der Gas- und Festphaseam Eingang 1
        self.A_g1, self.A_s1 = self.gas.gasarea(self.G_1), self.solid.solidarea(self.G_1)
        # Gib die Enthalpie der Festphase am Eingang 1  
        self.h_s1 = self.solid.h_s1
        # Gib die Geschwindigkeiten der Gas- und Festphasen am Eingang 1  
        vl_1 = vl(self.m_g1, self.m_s1, self.A_g1, self.A_s1, self.d_g1, self.d_s1, self.e_s)
        self.U_g1, self.U_s1 = vl_1[0], vl_1[1]
        
        # Gib die Schaetzwerte des kontrollvolumens 
        self.T_g2, self.T_s2 = self.T_g1+1e-6, self.T_s1
        self.d_g2, self.d_s2, self.d_g, self.d_s = self.d_g1, self.d_s1, self.d_g1, self.d_s1              
        self.vs = self.gas.propm(self.p_1, self.T_g1+1)['vs']
        self.cp_g = self.gas.propm(self.p_1, self.T_g1+1)['cp']
        self.L = self.gas.propm(self.p_1, self.T_g1+1)['L']
        self.cp_s = self.solid.propm(self.p_1)['cp']
        self.p_2, self.h_g2 = self.p_1, self.h_g1  
        
        self.h_subs = trendprop('co2', input_code = 'PSUBS+', prop1 = self.p_1, prop2 = 1, err_capt = True, 
                           unit_p = 'Pa', unit_T = 'K').Prop('H', unit_T = 'K',unit_m = 'kg')
        # Gib die AusgangsEnthalpie von Gasphase beim Wandtemperatur  
        try:
            self.h_w = trendprop('co2', input_code = 'TP', prop1 = self.T_w, prop2 = self.p_1, err_capt = True, 
                                 unit_p = 'Pa', unit_T = 'K').Prop('H', unit_T = 'K',unit_m = 'kg')
        except:
            self.h_w = trendprop('co2', input_code = 'TP', prop1 = 217, prop2 = self.p_1, err_capt = True, 
                                 unit_p = 'Pa', unit_T = 'K').Prop('H', unit_T = 'K',unit_m = 'kg')

    # Gib die berechnete Groesse an 
    def prop_out(self):
        # Gib die Massestrome der Gas- und Festphase am Ausgang 2 
        self.m_g2, self.m_s2 = self.m_g1 + self.m_sub, self.m_s1 - self.m_sub

        # Gib die Temperatur der Gas- und Festphase am Austritt 2
        self.T_g2, self.T_s2 = self.gas.prop2(self.p_2, self.h_g2)['T'], self.solid.prop2(self.p_2)['T']
        # Gib die Dichte der Gas- und Festphase am Austritt 2
        self.d_g2, self.d_s2 = self.gas.prop2(self.p_2, self.h_g2)['d'], self.solid.prop2(self.p_2)['d']
        # Gib die Querschnittsflaeche der Gas- und Festphase am Austritt 2 
        self.A_g2, self.A_s2 = self.gas.gasarea(self.G_2), self.solid.solidarea(self.G_2)
        
        # Gib die Geschwindigkeiten der Gas- und Festphasen am Ausgang 2 
        vl_2 = vl(self.m_g2, self.m_s2, self.A_g2, self.A_s2, self.d_g2, self.d_s2, self.e_s)
        self.U_g2, self.U_s2 = vl_2[0], vl_2[1]
        
        # Gib die durchschnittliche relative Flaeche der Gasphase
        self.G_m = self.solid.Gm(self.G_2)
        # Gib die durchschnittliche Querschnittsflaeche Flaeche der Gas- und Festphase
        self.A_g, self.A_s = self.gas.gasarea(self.G_m), self.solid.solidarea(self.G_m)
        # Gib die durchschnittliche Geschwindigkeiten der Gas- und Festphase
        U_g, U_s = self.gas.gasvelocity(self.U_g1, self.U_g2), self.solid.solidvelocity(self.U_s1, self.U_s2)     
        # Gib die durchschnittliche Druck 
        self.p_m = self.gas.pressurem(self.p_2)
        # Gib die durchschnittliche Temperatur  der Gas- und Festphase
        self.T_m, self.T_s = self.gas.gastemperatur(self.T_g2), self.solid.solidtemperatur(self.T_s2) 
        # Gib die durchschnittliche Dichte der Gas- und Festphase
        self.d_g, self.d_s = self.gas.propm(self.p_m, self.T_m)['d'], self.solid.propm(self.p_m)['d'] 
        # Gib die Viskositaet 
        self.vs = self.gas.propm(self.p_m, self.T_m)['vs']
        # Gib die isobare Waermekapazitaet
        self.cp_g = self.gas.propm(self.p_m, self.T_m)['cp']
        # Gib die Waermeleitfaehigkeit
        self.L = self.gas.propm(self.p_m, self.T_m)['L']
        # Gib die isobare Waermekapazitaet
        self.cp_s = self.solid.propm(self.p_m)['cp']
        # Gib die delta
        delta = self.cp_s / self.cp_g
        # Gib die eta 
        eta = 0.01 # m_s / m_g
        
        # Gib die Wand-Gas-Oberflaeche
        self.A_wg = self.gas.gaswallarea(self.G_2)
        # Gib die Querschnittsfläche von Feststoff- und Gasphasen
        self.A_Q = self.gas.gassolidarea(self.G_2)
        # Gib die Wand-Feststoff-Oberflaeche
        self.A_ws = self.solid.solidwallarea(self.G_2)
        # Gib die aequivalentslaenge
        self.d_h = self.gas.equradius(self.G_m)
        
        # Gib die dimensionlose Prandtl-Zahl
        self.Pr = pr(self.d_g, self.vs, self.cp_g, self.L)
        # Gib die dimensionlose Reynolds Zahl
        self.Re = re(U_g, self.d_h, self.vs)[0]
        # Gib die dimensionlose Nusselt Zahl
        self.Nu = nu(self.Re, self.Pr, self.R, self.H)
        # Gib die Waermeuebergangskoeffizient in der Gasphase
        self.alpha_g = alg(self.Nu, self.L, self.d_h)[0]
        # Gib die Waermeuebergangskoeffizient in suspension alpha_s
        self.alpha_s = als(self.alpha_g, self.Re, eta, delta)
        # Gig die Effektive Wärmeleitfähigkeit des Solidphase
        self.k = ttm()
        
        # Gib die Fanning friction factor for pure gas f_g 
        self.f_g = fg(self.Re)
        # Gib die friction factor for solid f_p 
        f_p = fp()
        
        # Gib die Querkraft und Reibungskraft
        self.F_Q = sf(self.G_m, self.d_g, self.d_s, U_g, U_s, self.e_s, self.R, self.H)
        self.FR_wg = wgff(self.A_wg, self.d_g, self.f_g, U_g, self.R, self.H)
        self.FR_ws = wsff(self.A_ws, self.d_s, self.d_g, self.e_s, self.H, f_p)
        # Gib die Enthalpie auf der sublimationskurve
        self.h_s2 = self.solid.prop2(self.p_2)['h']
        self.h_subs = trendprop('co2', input_code = 'PSUBS+', prop1 = self.p_2, prop2 = 1, err_capt = True, 
                           unit_p = 'Pa', unit_T = 'K').Prop('H', unit_T = 'K',unit_m = 'kg')
        # Gib die Enthalpie auf der Resublimationskurve
        self.h_subv = trendprop('co2', input_code = 'PSUBV+', prop1 = self.p_2, prop2 = 1, err_capt = True, 
                           unit_p = 'Pa', unit_T = 'K').Prop('H', unit_T = 'K',unit_m = 'kg')
        # Gib die erforderliche Enthalpieänderung für feste Sublimation 
        self.Dh_sub = self.h_subv - self.h_s1
        
        # Gib den von Wand auf Festphase uebertragenen Waermestrom
        self.Q_ws = self.solid.Wallsolidheattransfer(self.k, self.T_w, self.T_s, self.A_ws)
        # Gib die erforderliche Sublimationswärme für Feststoffe 
        self.Q_sub = self.m_sub * self.Dh_sub
        # Gib den von Wand auf Gasphase uebertragenen Waermestrom
        self.Q_wg = self.gas.Wallgasheattransfer(self.alpha_s, self.T_w, self.T_m, self.A_wg)
        # Gib den von Gasphase auf Festphase uebertragenen Waermestrom
        self.Q_gs = self.gas.Gassolidheattransfer(self.k, self.T_s, self.T_m, self.A_Q)
            
        prop_out = {
            'fluid' : self.fluid,
            'G' : self.G_2,
            'T_g' : self.T_g2,
            'T_s' : self.T_s2,
            'h_g' : self.h_g2, 
            'p' : self.p_2,
            'm_g' : self.m_g2,
            'm_s' : self.m_s2,
            'm_sub' : self.m_sub,
            'Q_gs' : self.Q_gs,
            'Q_wg' : self.Q_wg,
            'Q_ws' : self.Q_ws,
            'Q_sub' : self.Q_sub,
            'F_Q' : self.F_Q,
            'FR_wg': self.FR_wg,
            'FR_ws': self.FR_ws,
            }
        
        return prop_out 
    
    # Funktion fuer die Iteration
    def func_1(self, out):
        
        # Gib die zu berechnende Parameter 
        [m_sub, p_2, G_2, h_g2, h_s2] = out

        # Gib die Massestrome der Gas- und Festphase am Ausgang 2 
        m_g2 = self.m_g1 + m_sub
        m_s2 = self.m_s1 - m_sub
        
        # Gib die durchschnittliche Massestrome der Gas- und Festphase 
        # m_g = self.gas.gasmass(m_g2)
        # m_s = self.gas.gasmass(m_s2)

        # Gib die Temperatur der Gas- und Festphase am Austritt 2
        try:
            T_g2 = self.gas.prop2(p_2, h_g2)['T']
            self.T_g2 = T_g2
        except:
            T_g2 = self.T_g2
            
        T_s2 = self.T_s2
        try:
            T_s2 = self.solid.prop2(p_2)['T']
            self.T_s2 = T_s2
        except:
            T_s2 = self.T_s2
            
        # Gib die Dichte der Gas- und Festphase am Austritt 2
        try:
            d_g2 = self.gas.prop2(p_2, h_g2)['d']
            self.d_g2 = d_g2 
        except:
            d_g2 = self.d_g2
        try:
            d_s2 = self.solid.prop2(self.p_2)['d']
            self.d_s2 = d_s2
        except:
            d_s2 = self.d_s2
        # Gib die Querschnittsflaeche der Gas- und Festphase am Austritt 2 
        A_g2 = self.gas.gasarea(G_2)
        A_s2 = self.solid.solidarea(G_2)
        # Gib die Geschwindigkeiten der Gas- und Festphasen am Ausgang 2 
        vl_2 = vl(m_g2, m_s2, A_g2, A_s2, d_g2, d_s2, self.e_s)
        U_g2 = vl_2[0]
        U_s2 = vl_2[1]
        # Gib die durchschnittliche relative Flaeche der Gasphase
        G_m = self.solid.Gm(G_2)
        # Gib die durchschnittliche Querschnittsflaeche Flaeche der Gas- und Festphase
        # A_g = self.gas.gasarea(G_m)
        # A_s = self.solid.solidarea(G_m)
        # Gib die durchschnittliche Geschwindigkeiten der Gas- und Festphase
        U_g = self.gas.gasvelocity(self.U_g1, U_g2)      
        U_s = self.solid.solidvelocity(self.U_s1, U_s2)     
        # Gib die durchschnittliche Druck 
        p_m = self.gas.pressurem(p_2)
        # Gib die durchschnittliche Temperatur  der Gas- und Festphase
        T_m = self.gas.gastemperatur(T_g2)
        try:
            T_s = self.solid.solidtemperatur(T_s2) 
        except:
            T_s = self.T_s1
        # Gib die durchschnittliche Dichte der Gas- und Festphase
        try:
            d_g = self.gas.propm(p_m, T_m)['d']
            self.d_g = d_g
        except:
            d_g = self.d_g            
        try:
            d_s = self.solid.propm(p_m)['d']
            self.d_s = d_s
        except:
            d_s = self.d_s
    
        # Gib die Viskositaet 
        try:
            vs = self.gas.propm(p_m, T_m)['vs']
            self.vs = vs
        except:
            vs = self.vs
            
        # Gib die isobare Waermekapazitaet
        try:
            cp_g = self.gas.propm(p_m, T_m)['cp']
            self.cp_g = cp_g
        except:
            cp_g = self.cp_g
            
        # Gib die Waermeleitfaehigkeit
        try:
            L = self.gas.propm(p_m, T_m)['L']
            self.L = L
        except:
            L = self.L
            
        # Gib die isobare Waermekapazitaet
        try:
            cp_s = self.solid.propm(p_m)['cp']
            self.cp_s = cp_s
        except:
            cp_s = self.cp_s
        # Gib die delta
        delta = cp_s / cp_g
        # Gib die eta 
        eta = 0.01 # m_s / m_g
        
        # Gib die Wand-Gas-Oberflaeche
        A_wg = self.gas.gaswallarea(G_2)
        # Gib die Querschnittsfläche von Feststoff- und Gasphasen
        A_Q = self.gas.gassolidarea(G_2)
        # Gib die Wand-Feststoff-Oberflaeche
        A_ws = self.solid.solidwallarea(G_2)
                        
        # Gib die aequivalentslaenge
        d_h = self.gas.equradius(G_m)
        # Gib die dimensionlose Prandtl-Zahl
        Pr = pr(d_g, vs, cp_g, L)
        # Gib die dimensionlose Reynolds Zahl
        Re = re(U_g, d_h, vs)[0]
        # Gib die dimensionlose Nusselt Zahl
        Nu = nu(Re, Pr, self.R, self.H)
        # Gib die Waermeuebergangskoeffizient in der Gasphase
        alpha_g = alg(Nu, L, d_h)[0]
        # Gib die Waermeuebergangskoeffizient in suspension alpha_s
        try:
            alpha_s = als(alpha_g, Re, eta, delta)    
        except:
            alpha_s = alpha_g + 1e-1
        # Gig die Effektive Wärmeleitfähigkeit des Solidphase
        k = ttm()
        
        # Gib die Fanning friction factor for pure gas f_g 
        f_g = fg(Re)
        # Gib die friction factor for solid f_p 
        f_p = fp()
        
        # Gib die Querkraft und Reibungskraft
        F_Q = sf(G_m, d_g, d_s, U_g, U_s, self.e_s, self.R, self.H)
        FR_wg = wgff(A_wg, d_g, f_g, U_g, self.R, self.H)
        FR_ws = wsff(A_ws, d_s, d_g, self.e_s, self.H, f_p)
            
        # Gib die Enthalpie auf der Sublimationskurve
        try: 
            h_subs = trendprop('co2', input_code = 'PSUBS+', prop1 = p_2, prop2 = 1, err_capt = True, 
                           unit_p = 'Pa', unit_T = 'C').Prop('H', unit_T = 'C',unit_m = 'kg')
            self.h_subs = h_subs
        except:
            h_subs = self.h_subs
            
        # print(h_subs)
        # Gib die Enthalpie auf der Resublimationskurve
        try:
            h_subv = trendprop('co2', input_code = 'PSUBV+', prop1 = p_2, prop2 = 1, err_capt = True, 
                           unit_p = 'Pa', unit_T = 'C').Prop('H', unit_T = 'C',unit_m = 'kg')
            self.h_subv = h_subv
        except:
            h_subv = self.h_g1
        # Gib die erforderliche Enthalpieänderung für feste Sublimation 
        Dh_sub = h_subv - self.h_s1

        # Gib den von Wand auf Festphase uebertragenen Waermestrom
        Q_ws = self.solid.Wallsolidheattransfer(k, self.T_w, T_s, A_ws)
        
        # Randbedingung des Modells
        h_s2 = h_subs
        #m_sub = -(d_s2 * A_s2 * U_s2 - self.d_s1 * self.A_s1 * self.U_s1) * (1 - self.e_s)
        
        # Gib den von Wand auf Gasphase uebertragenen Waermestrom
        Q_wg = self.gas.Wallgasheattransfer(alpha_s, self.T_w, T_m, A_wg)     
        # Gib den von Gasphase auf Festphase uebertragenen Waermestrom
        Q_gs = self.gas.Gassolidheattransfer(k, T_s, T_m, A_Q)
        
        try:
            # Bilanzen eintragen 
            func = [# Impulsbilanz der Gasphase 
                m_g2 * U_g2 - self.m_g1 * self.U_g1 + p_2 * A_g2 - self.p_1 * self.A_g1 + F_Q + FR_wg,
                # Impulsbilanz der Festphase
                m_s2 * U_s2 - self.m_s1 * self.U_s1 + p_2 * A_s2 - self.p_1 * self.A_s1 - F_Q + FR_ws,
                # Energiebilanz der Gasphase 
                m_g2 * h_g2 - self.m_g1 * self.h_g1 - m_sub * h_g2 + Q_gs - Q_wg,
                # m_g2 * h_g2 - self.m_g1 * self.h_g1 + Q_gs - Q_wg,
                # Energiebilanz der Festphase
                m_s2 * h_subs - self.m_s1 * self.h_s1 + m_sub * h_g2 - Q_ws - Q_gs,
                # m_s2 * h_subs - self.m_s1 * self.h_s1 - Q_ws - Q_gs,
                # Massenbilanzen der Gas-Festphase
                # h_s2 - h_subs
                m_sub + (d_s2 * A_s2 * U_s2 - self.d_s1 * self.A_s1 * self.U_s1) * (1 - self.e_s)
                ]
        except:
            func = [# Impulsbilanz der Gasphase 
                m_g2 * U_g2 - self.m_g1 * self.U_g1 + p_2 * A_g2 - self.p_1 * self.A_g1 + F_Q + FR_wg,
                # Impulsbilanz der Festphase
                m_s2 * U_s2 - self.m_s1 * self.U_s1 + p_2 * A_s2 - self.p_1 * self.A_s1 - F_Q + FR_ws,
                # Energiebilanz der Gasphase 
                m_g2 * h_g2 - self.m_g1 * self.h_g1 - m_sub * h_g2 + Q_gs - Q_wg,
                # Energiebilanz der Festphase
                m_s2 * h_subs - self.m_s1 * self.h_s1 + m_sub * h_g2 - Q_ws - Q_gs,
                # Massenbilanzen der Gas-Festphase
                h_s2 - self.h_subs
                # m_sub + (d_s2 * A_s2 * U_s2 - self.d_s1 * self.A_s1 * self.U_s1) * (1 - self.e_s)
                ]
            
        # print(p_2, T_g2, T_s2, h_g2, h_s2)
        
        return func

    # Funktion fuer die Iteration
    def func_2(self, out):
        
        # Gib die zu berechnende Parameter 
        [m_sub, p_2, G_2, h_g2] = out

        # Gib die Massestrome der Gas- und Festphase am Ausgang 2 
        m_g2 = self.m_g1 + m_sub
        m_s2 = self.m_s1 - m_sub
        # Gib die Querschnittsflaeche der Gas- und Festphase am Austritt 2 
        A_g2 = self.gas.gasarea(G_2)
        A_s2 = self.solid.solidarea(G_2)
        
        # Gib die Temperatur der Gas- und Festphase am Austritt 2
        try:
            T_g2 = self.gas.prop2(p_2, h_g2)['T']
            self.T_g2 = T_g2
        except:
            T_g2 = self.T_g2
            
        T_s2 = self.T_s2
        try:
            T_s2 = self.solid.prop2(p_2)['T']
            self.T_s2 = T_s2
        except:
            T_s2 = self.T_s2
            
        # Gib die Dichte der Gas- und Festphase am Austritt 2
        try:
            d_g2 = self.gas.prop2(p_2, h_g2)['d']
            self.d_g2 = d_g2 
        except:
            d_g2 = self.d_g2
        try:
            d_s2 = self.solid.prop2(self.p_2)['d']
            self.d_s2 = d_s2
        except:
            d_s2 = self.d_s2

        # Gib die Geschwindigkeiten der Gas- und Festphasen am Ausgang 2 
        vl_2 = vl(m_g2, m_s2, A_g2, A_s2, d_g2, d_s2, self.e_s)
        U_g2 = vl_2[0]
        U_s2 = vl_2[1]
        # Gib die durchschnittliche relative Flaeche der Gasphase
        G_m = self.solid.Gm(G_2)
 
        # Gib die durchschnittliche Geschwindigkeiten der Gas- und Festphase
        U_g = self.gas.gasvelocity(self.U_g1, U_g2)      
        U_s = self.solid.solidvelocity(self.U_s1, U_s2)     
        # Gib die durchschnittliche Druck 
        p_m = self.gas.pressurem(p_2)
        # Gib die durchschnittliche Temperatur  der Gas- und Festphase
        T_m = self.gas.gastemperatur(T_g2)
        try:
            T_s = self.solid.solidtemperatur(T_s2) 
        except:
            T_s = self.T_s1
        # Gib die durchschnittliche Dichte der Gas- und Festphase
        try:
            d_g = self.gas.propm(p_m, T_m)['d']
            self.d_g = d_g
        except:
            d_g = self.d_g            
        try:
            d_s = self.solid.propm(p_m)['d']
            self.d_s = d_s
        except:
            d_s = self.d_s
    
        # Gib die Viskositaet 
        try:
            vs = self.gas.propm(p_m, T_m)['vs']
            self.vs = vs
        except:
            vs = self.vs
            
        # Gib die isobare Waermekapazitaet
        try:
            cp_g = self.gas.propm(p_m, T_m)['cp']
            self.cp_g = cp_g
        except:
            cp_g = self.cp_g
            
        # Gib die Waermeleitfaehigkeit
        try:
            L = self.gas.propm(p_m, T_m)['L']
            self.L = L
        except:
            L = self.L
            
        # Gib die isobare Waermekapazitaet
        try:
            cp_s = self.solid.propm(p_m)['cp']
            self.cp_s = cp_s
        except:
            cp_s = self.cp_s
        # Gib die delta
        delta = cp_s / cp_g
        # Gib die eta 
        eta = 0.01 # m_s / m_g
        
        # Gib die Wand-Gas-Oberflaeche
        A_wg = self.gas.gaswallarea(G_2)
        # Gib die Querschnittsfläche von Feststoff- und Gasphasen
        A_Q = self.gas.gassolidarea(G_2)
        # Gib die Wand-Feststoff-Oberflaeche
        A_ws = self.solid.solidwallarea(G_2)
                        
        # Gib die aequivalentslaenge
        d_h = self.gas.equradius(G_m)
        # Gib die dimensionlose Prandtl-Zahl
        Pr = pr(d_g, vs, cp_g, L)
        # Gib die dimensionlose Reynolds Zahl
        Re = re(U_g, d_h, vs)[0]
        # Gib die dimensionlose Nusselt Zahl
        Nu = nu(Re, Pr, self.R, self.H)
        # Gib die Waermeuebergangskoeffizient in der Gasphase
        alpha_g = alg(Nu, L, d_h)[0]
        # Gib die Waermeuebergangskoeffizient in suspension alpha_s
        try:
            alpha_s = als(alpha_g, Re, eta, delta)    
        except:
            alpha_s = alpha_g + 1e-1
        # Gig die Effektive Wärmeleitfähigkeit des Solidphase
        k = ttm()
        
        # Gib die Fanning friction factor for pure gas f_g 
        f_g = fg(Re)
        # Gib die friction factor for solid f_p 
        f_p = fp()
        
        # Gib die Querkraft und Reibungskraft
        F_Q = sf(G_m, d_g, d_s, U_g, U_s, self.e_s, self.R, self.H)
        FR_wg = wgff(A_wg, d_g, f_g, U_g, self.R, self.H)
        FR_ws = wsff(A_ws, d_s, d_g, self.e_s, self.H, f_p)
            
        # Gib die Enthalpie auf der Sublimationskurve
        try: 
            h_subs = trendprop('co2', input_code = 'PSUBS+', prop1 = p_2, prop2 = 1, err_capt = True, 
                           unit_p = 'Pa', unit_T = 'C').Prop('H', unit_T = 'C',unit_m = 'kg')
            self.h_subs = h_subs
        except:
            h_subs = self.h_subs
            
        # print(h_subs)
        # Gib die Enthalpie auf der Resublimationskurve
        try:
            h_subv = trendprop('co2', input_code = 'PSUBV+', prop1 = p_2, prop2 = 1, err_capt = True, 
                           unit_p = 'Pa', unit_T = 'C').Prop('H', unit_T = 'C',unit_m = 'kg')
            self.h_subv = h_subv
        except:
            h_subv = self.h_g1
        # Gib die erforderliche Enthalpieänderung für feste Sublimation 
        Dh_sub = h_subv - self.h_s1

        # Gib den von Wand auf Festphase uebertragenen Waermestrom
        Q_ws = self.solid.Wallsolidheattransfer(k, self.T_w, T_s, A_ws)
        
        # Randbedingung des Modells
        h_s2 = self.solid.prop2(p_2)['h']
        m_sub = -(d_s2 * A_s2 * U_s2 - self.d_s1 * self.A_s1 * self.U_s1) * (1 - self.e_s)
        
        # Gib den von Wand auf Gasphase uebertragenen Waermestrom
        Q_wg = self.gas.Wallgasheattransfer(alpha_s, self.T_w, T_m, A_wg)     
        # Gib den von Gasphase auf Festphase uebertragenen Waermestrom
        Q_gs = self.gas.Gassolidheattransfer(k, T_s, T_m, A_Q)
        
        print(p_2, T_g2)
        
        # Bilanzen eintragen 
        func = [# Impulsbilanz der Gasphase 
                m_g2 * U_g2 - self.m_g1 * self.U_g1 + p_2 * A_g2 - self.p_1 * self.A_g1 + F_Q + FR_wg,
                # Impulsbilanz der Festphase
                m_s2 * U_s2 - self.m_s1 * self.U_s1 + p_2 * A_s2 - self.p_1 * self.A_s1 - F_Q + FR_ws,
                # Energiebilanz der Gasphase 
                m_g2 * h_g2 - self.m_g1 * self.h_g1 - m_sub * h_g2 + Q_gs - Q_wg,
                # m_g2 * h_g2 - self.m_g1 * self.h_g1 + Q_gs - Q_wg,
                # Energiebilanz der Festphase
                m_s2 * h_subs - self.m_s1 * self.h_s1 + m_sub * h_g2 - Q_ws - Q_gs,
                # m_s2 * h_subs - self.m_s1 * self.h_s1 - Q_ws - Q_gs,
                ]
        
        return func

    # Lösen die Gleichungen
    def solve_out(self):
        
        # Gib den Anfangswert des Unbekannten an
        m_sub_0 = self.m_sub_1 #+ 0.151e-7
        p_2_0 = self.p_1 - 205
        G_2_0 = self.G_1 
        h_g2_0 = self.h_g1 
        # h_s2_0 = self.h_s1
        
        # Ausnahmebehandlung
        try:
            # Die root-Funktion benutzen
            self.resr = root(self.func_2, [m_sub_0, p_2_0, G_2_0, h_g2_0],
                       method='hybr',
                       tol=1e-10)
                       # options={'maxiter':2}) #hybr
            
            # Die fsolve-Funktion benutzen
            # resf = fsolve(self.func_3, [m_sub_0, p_2_0, G_2_0, h_g2_0],factor=100)
                 
            # Die Ergebnisse dem entsprechenden Parameter zuordnen 
            self.m_sub, self.p_2, self.G_2, self.h_g2 = self.resr.x
            
            # Kontrolliren den Ausgangswert
            self.prop_out()
            self.check_out()
            
        except ValueError:
            raise #("Problem auftritt")
        else:
            print("Die Berechnung wird einmal durchgeführt")
        
    # Überprüfen den Ausgangswert
    def check_out(self):
        # Kontrollieren den Eingangstemperatur 
        if self.T_g1 > self.T_w:
            raise ValueError("T_g1 cannot exceed the temperature T_w")      
        # Kontrollieren die Ausgangstemperatur 
        if self.T_g2 > self.T_w:
            raise ValueError("T_g2 cannot exceed the temperature T_w")            
        # Kontrollieren den Eingangsdruck
        p_trip = 1e6 * Trip(self.fluid)[1]
        if self.p_1 > p_trip:
            raise ValueError("pressure p cannot exceed the triplepressure")        
        # Kontrollieren den Ausgangsdruck 
        if self.p_2 > self.p_1:
            raise ValueError("There is only loss of pressure, not increase")
        # Kontrollieren die relative Fläche am Austritt 2 
        if (self.G_2 > 1 or self.G_2 < 0):
            raise ValueError("G can only be a number between 0 and 1 ")
        # Kontrollieren den Ausgangsmassestrom von Festphase 
        if self.m_s2 < 0:
            raise ValueError("m_s cannot be less than 0")
        # Kontrollieren den sublimierten Feststoff
        if self.m_sub < 0 :
            raise ValueError("m_sub cannot be less than 0")

    
if __name__ == '__main__':
    
    start = timeit.default_timer()
    
    rohr = Rohrelemente(R=1.5e-3, H=2, T_w=213, n=1000)
    prop_in = Prop_in(fluid='co2', m_g_in=0.05e-3, m_s_in=0.05e-3, p_in=1.01325e5, 
                      h_g_in= trendprop('co2', input_code = 'PSUBV+', prop1 = 1.01325e5,
                                      prop2 = 1, err_capt = True, unit_p = 'Pa', unit_T = 'C')
                                      .Prop('H', unit_T = 'C',unit_m = 'kg'),
                      G_in=0.9) 
    
    # cvm = controllvolume_minimize(prop_in, rohr)
    
    cv = controllvolume(prop_in, rohr)
    cv.solve_out()
    prop_out = cv.prop_out()
    
    stop = timeit.default_timer()
    print('Time: ', stop - start)
    