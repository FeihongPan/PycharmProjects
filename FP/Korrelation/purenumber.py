# Gib die dimensionlose Groesse
import math

# Berechnung der Renoulds-Zahl Re
"""
    U:Geschwindigkeit                                                                [m/s]
    d_h:charakteristische Laenge der aequivalente bzw. der hydraulische Durchmesser  [m]
    vs:kinematische Viscocity                                                        []   
        
"""
def Renouldsnumber(U, d_h, vs):
    Re = U * d_h / vs
    return Re


# Berechnung der Prandtl-Zahl Pr
"""
    d:density                                                                        [kg/m3]
    vs:kinematische Viscocity                                                        []
    cp:isobaric_heat_capacity(isobare Waermekapazitaet)                              [J/kg/K]
    L:thermal_conductivity(Waermeleitkoeffizient)                                    [W/m/K]   
        
"""
def Prandtlnumber(d, vs, cp, L):
    a = L / (d * cp)
    Pr = vs / a
    return Pr


# Berechnung der Nusselt-Zahl Nu
"""
    GL(1)-(3):
    Verein Deutscher Ingenieure VDI (Hrsg.): VDI-Wärmeatlas. 8. Auflage auf CD-ROM,
    Springer Verlag, Berlin Heidelberg 1997.
            
    GL(4):
    Wagner, W.: Wärmeübertragung. Vogelverlag, 1981.
                 
"""     
def Nusseltnumber(Re, Pr, R, H):
    D = 2 * R
    # laminar thermisch beliebig, hydraulisch ausgebildet GL(1)
    if Re <= 2300:
        Nu = (3.66 ** 3 + (1.61 ** 3 * Re * Pr) * D / H) ** (1/3)        
        
    # turbulent konstante Wandtemperatur,hydraulisch ausgebildet GL(3)       
    elif 1e4 <= Re <= 1e6:
            eta = (1.8 * math.log(Re, 10) - 1.5) ** -2
            Nu = (eta/8) * Re * Pr * (1 + (D/H) ** (2/3)) / (1 + 12.7 * (eta/8) ** 0.5 * (Pr ** (2/3) - 1))
    
    # Stroemung im Uebergangsbereich
    else:
        eta = (0.79 * math.log(Re) - 1.64) ** -2
        Nu_t = (eta/8) * (Re - 1000) * Pr * (1 + (D/H) ** (2/3)) / (1 + 12.7 * (eta/8) ** 0.5 * (Pr ** (2/3) - 1))
        Nu_l = (3.66 ** 3 + (1.61 ** 3 * Re * Pr) * D / H) ** (1/3)  # 0.664 * (Re * D * Pr ** 0.33 / H) ** 0.5
        x = (Re - 2300) / (1e4 - 2300)
        
        Nu = ((1 - x) * Nu_l + x * Nu_t)
    return Nu    
        
        
        
if __name__ == '__main__':         
    Nu = Nusseltnumber(Re=1000, Pr=2.2808, R=1.5e-3, H=2)
        
        
        
        
        
