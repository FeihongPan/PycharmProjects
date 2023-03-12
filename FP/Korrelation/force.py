"""
    f_p:solid friction factor    
    f_g:gas friction factor 
    d_s:density of solid                [kg/m3]
    d_g:density of gas                  [kg/m3]
    e_s:voigate of solid-layer          []
    A_s:solid area                      [m2]
    A_g:gas area                        [m2]
    H:rohr length                       [m]
    g:Gravitational acceleration        [m/s2]
"""
# Gib die Wand-Festphase-Reibungskraft FR_ws
def wallgasfrictionforce(A_wg, d_g, f_g, U_g, R, H):
    """
    1966_Pfeffer_Analysis and correlation of and 
    friction factor data for heat-transfer coefficient 
    dilute gas-solid suspension
    
    """
    FR_wg = f_g * H / R * d_g * A_wg * U_g ** 2 
    
    return FR_wg

# Gib die Wand-Festphase-Reibungskraft FR_ws
def wallsolidfrictionforce(A_ws, d_s, d_g, e_s, H, f_p):
    """
        'HANDBOOK OF CONVEYING AND HANDLING OF PARTICU LATE SOLIDS'
            Edited by
            A. LEVY
                Ben-Gurion University of the Negev, Department of Mechanical
                Engineering, Beer Sheva 84105, Israel
                
    """
    """
            Operating limits of low-velocity pneumatic conveying
            J. Yi and P.W. Wypych
            Seite 346 Gleichung (4)
    """  
    FR_ws = f_p * (d_s - d_g) * (1 - e_s) * A_ws * H 
    
    return FR_ws

# Gib die Wand-Festphase-Reibungskraft FR_ws
def wallsolidfrictionforce2(m_s, H, f_p, U_s, g = 9.80665):
    """
    Muschelknautz E, Krambrock W (1969) Chem-Ing-Techn 41:1164–1172
    """  
    FR_ws = f_p * m_s * g * H / U_s
    
    return FR_ws


# Gib die Querkraft
"""
    L_h:momentum transfer factor
    
    Operating limits of low-velocity pneumatic conveying
    J. Yi and P.W. Wypych
    Seite 348 Gleichung (15)
"""

def shearforce0(A_Q, U_g, vs, d_g, R):
    F_Q = A_Q * vs * d_g * U_g / R
    return F_Q

def shearforce(G, d_g, d_s, U_g, U_s, e_s, R, H, L_h = 0.0826):
    F_Q = L_h * d_g * U_g ** 2 * 2 * R * H * (1 - U_s/U_g) * (1 - d_g * e_s / (d_s * (1 - e_s))) * ( (4 * G * (1 - G)) **(1/3) / G ** 2)
    return F_Q



if __name__ == '__main__':
    FR_wg = wallgasfrictionforce(A_wg=0.002328, d_g=3, f_g=0.008776, U_g=5, R=0.005, H=1)
    FR_ws = wallsolidfrictionforce(A_ws=0.0008134, d_s=1570, d_g=3, e_s= 0.1, H=0.1, f_p = 0.1) 
    FR_ws2 = wallsolidfrictionforce2(m_s=1e-3, H=0.1, f_p=0.8, U_s=0.1)
    FQ_0 = shearforce0(A_Q=0.00073, U_g=5, vs=1.5e-5, d_g=3, R=0.005)
    FQ = shearforce(G=0.9, d_g=3, d_s=1570, U_g=5, U_s=0.1, e_s=0.1, R=0.005, H=0.1)    
    
    print(FR_wg, FR_ws, FR_ws2, FQ_0, FQ,) 
    
    
    
    
    
    
