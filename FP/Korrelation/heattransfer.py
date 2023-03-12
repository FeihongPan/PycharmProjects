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
# Gib die Waermeuebergangskoeffizient in pure gas alpha_g (heat transfer coefficient)
def alpha_g(Nu, L, d_h):
    alpha_g = L * Nu / d_h
    return alpha_g

# Gig die Effektive Wärmeleitfähigkeit des Solidphase
"""
Dvornitsyn et al, 2006, Experimental investigation of a bottle-sublimation cooler
"""
def thermaltransmittance():
    k = 200
    return k

# Gib die Waermeuebergangskoeffizient in suspension alpha_s (heat transfer coefficient)

"""
    p:partikel                            
    g:gas phase
    k:functions of other variables 
    alpha:heat transfer coefficient       [W/m2/K]
    eta:loading ratio
    delta:specific heat ratio
    Re:Reynolds number
    Einheit von Partikeldurchmesser:      [1e-6m]
"""

# Einfluss auf Loading ratio eta und specific heat ratio delta         
"""
    1966_Pfeffer_Analysis and correlation of and 
    friction factor data for heat-transfer coefficient 
    dilute gas-solid suspension
    
"""
     
# Seite 20 Gleichung (38) (15)
def farbarandmorley(alpha_g, Re, eta, delta):
    if eta < 10:
        alpha_s = alpha_g * 6.4 * (Re ** -0.2) * ((eta * delta) ** 0.45)
    else:
        alpha_s = alpha_g * 6.8 * (Re ** -0.2) * (eta ** 0.45)
    return alpha_s
 
    
# Seite 21 Gleichung (39) (17)
def Danziger(alpha_g, Re, eta, delta):
    if eta < 10:
        alpha_s = alpha_g * 4.0 * (Re ** -0.14) * ((eta * delta ** 0.45))
    else:
        alpha_s = alpha_g * 3.7 * (Re ** -0.14) * (eta ** 0.45)
    return alpha_s
    
    
# Seite 21 Gleichung (19)
def Schluderberg(alpha_g, eta, delta):
    alpha_s = alpha_g * 0.78 * ((1 + eta * delta ** 0.45))
    return alpha_s
    

# Seite 21 Gleichung (25)
def FranklinInstitute(alpha_g, Re, eta, delta):
    alpha_s = alpha_g * 16.9 * (Re ** -0.3) * ((1 + eta * delta) ** 0.45)
    return alpha_s

    
# Seite 21 Gleichung (34)
def GorbisandBakhtiozin(alpha_g, Re_g, Re_p, eta, delta):
    alpha_s = alpha_g * (1 + 6.3 * (Re_g ** -0.3) * (Re_p ** -0.33))
    return alpha_s

    
# k-Form Gleichung,davon k ist die Funktion von der andere Variable 
# Seite 23 Gleichung (40)
def k1form(eta, delta, alpha_g, k_1):
    if 2 < eta < 10:
        alpha_s = alpha_g * k_1 * ((eta * delta) ** 0.45)
    else:
        raise ValueError()
    return alpha_s


# Seite 23 Gleichung (41)
def k4form(eta, delta, alpha_g, k_4):
    if 2 < eta < 10:
        alpha_s = alpha_g * (1 + k_4 * eta * delta)
    else:
        raise ValueError()
    return alpha_s
    

# Finale Form
# Seite 26 Gleichung (43)
def finalform1(alpha_g, Re, eta, delta):
    if Re < 1e6 and (eta * delta) < 10:
        alpha_s =  alpha_g * 7.6 * (Re ** -0.21) * ((eta * delta) ** 0.45) 
    else:
        raise ValueError()
    return alpha_s
    
# Seite 26 Gleichung (44)
def finalform2(alpha_g, Re, eta, delta):
    if Re < 1e6 and (eta * delta) < 10:
        alpha_s =  alpha_g * (1 + 4.0 * (Re ** -0.32) * eta * delta)   
    else:
        raise ValueError()
    return alpha_s    
    



if __name__ == '__main__':
    alpha_g = alpha_g(26.5, 2, 2)
    alpha_s = finalform2(alpha_g, 6000, 0.5, 1.8)
 

    
    
"""    

    # Einfluss auf den Duechmesser von Partikel  
    def EOparticlediameter:
        pass 

    # Einfluss auf Gas Renoulds Zahl 
    def EOgasrenouldsnumber:
        pass   
    
"""    