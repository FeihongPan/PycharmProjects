# Gib die Geschwindigkeit der Solid- und Gasphase
from scipy.optimize import root

"""
    m_s:solid mass flow             [kg/s]
    d_s:solid density               [kg/m3]
    A_s:solid area                  [m2]
    e_s:voigate of solid-layer      
    m_g:gas mass flow               [kg/s]
    d_g:gas density                 [kg/m3]
    A_g:gas area                    [m2]
    x:quality 
    a:void fraction                 
    S:slip rate
""" 
def velocity1(m_g, m_s, A_g, A_s, d_g, d_s, x_g ,e_s):
    def func(out):
        U_g = out[0]
        U_s = out[1]
        
        fun = [m_s - d_s * A_s * U_s * (1 - e_s) ,# - d_s * A_g * U_g * (1 - x_g)
               m_g - d_g * A_g * U_g * x_g - d_g * A_s * e_s * U_s]
        return fun
    
    U_g0 = 1
    U_s0 = 1
    
    res = root(func, [(U_g0, U_s0)])
    
    U_g, U_s = res.x
    return(U_g, U_s)


def velocity2(m_g, m_s, A_g, A_s, d_g, d_s, e_s):
    U_s = m_s / (d_s * A_s * (1 - e_s))
    U_g = (m_g - d_g * A_s * U_s * e_s) / (d_g * A_g)
    return (U_g, U_s)

if __name__ == '__main__':
    velocity1 = velocity1(1e-3, 1e-3, 7.06e-5, 7.85e-6, 3, 1571, 0.99 ,0.1)
    velocity2 = velocity2(1e-3, 1e-3, 7.06e-5, 7.85e-6, 3, 1571, 0.1)
    
