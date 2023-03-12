"""
    1966_Pfeffer_Analysis and correlation of and 
    friction factor data for heat-transfer coefficient 
    dilute gas-solid suspension
    
"""
# Gib die Fanning friction factor for pure gas f_g 
"""
    f_g:Fanning friction factor for pure gas        
    Selte 28 Gleichung (45)

    Churchill, S.W. (1977). "Friction factor equation spans all fluid-flow regimes". 
    Chemical engineering. 84 (24): 91–92

    Lightfoot, Edwin N.; Stewart, Warren E. (2007). Transport phenomena. Wiley

    Colebrook, C. F.; White, C. M. (3 August 1937). 
    "Experiments with Fluid Friction in Roughened Pipes". 
    Proceedings of the Royal Society of London. 
    Series A, Mathematical and Physical Sciences. 161 (906): 367–381.
"""
def fg(Re):
    if 1e4 <= Re <= 1e6:
        f_g = 0.046 * Re ** (-0.2)
    elif Re <= 2300:
        f_g = 16 / Re
    else:
        f_g = 0.0791 * Re ** (-0.25)
    return f_g

# Gib die friction factor for solid f_p 
def fp():
    f_p = 0.1
    return f_p

# Gib die Druckverlust Dp_s
"""
    gc:gravitational constant = 6.67408 × 10e-11    [m3/kg/s2]
    U_g:gas velocity                                [m/s]
    H:rohr length                                   [m]
    R:rohr radius                                   [m]
    f_s:Fanning friction factor based on velocity 
        of suspension and density of gas            []
    d_g:density of gas                              [kg/m3]
    Dp_s:pressure drop caused by suspension flow    [Pa]
    DP_g:pressure drop caused by gas                [Pa(N/m2)(kgm-1s-2)]
    
"""    
# Gib die f_s:Fanning friction factor based on velocity of suspension and density of gas f_s
"""
    eta:loading ratio = m_s/m_g
    delta:specific heat ratio
    
"""
# Seite 44 Gleichung (69b)
def JulianandDukler(f_g, eta):
    f_s = f_g * (1 + eta) ** 0.3
    return f_s


# Seite 44 Gleichung (70)
def Gorbis1(f_g, Re, eta):
    f_s = f_g * (1 + 4.0 * Re ** -0.32 * eta)
    return f_s


# Seite 44 Gleichung (71)
def Gorbis2(f_g, Re, eta, delta):
    f_s = f_g * (1 + 1/delta * (7.6 * (delta * eta) ** 0.45 * Re ** -0.21 -1))
    return f_s


# Seite 49 Gleichung (70)
def form1(f_g, Re, eta):
    f_s = f_g * (1 + 4.0 * eta * Re ** -0.32 )
    return f_s


# Seite 49 Gleichung (69b)
def form2(f_g, eta):
    f_s = f_g * (1 + eta) ** 0.3
    return f_s



if __name__ == '__main__':
    f_g = fg(6600)
    f_s = form1(f_g, 6600, 0.01)






       
    
    
    