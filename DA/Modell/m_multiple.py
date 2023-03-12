"""
Frequency Response Modell
+ Parameter:
    + delta: charakristische Funktion
    + eta: Parameter in der Frequency Response zur Substitution
    + omega: Winkelfrequenz [1/s]
    + f: Frequenz [1/s]
    + t_h, t_R: Zeitkonstante [s]
    + a: Kugelradius [m]
    + m: Probenmasse [kg]
    + L: Dicke von der ebenen Platte [m]
    + D: Diffusionskoeffizient [m^2/s]
    + R: Gaskonstante = 8.314 J/molK
    + T_0: Versuchstemperatur [K]
    + V_0: Ausgangsvolumen [m^3]

+ Index:
    + 0: Ausgangsdaten
    + 1: Platte, der Index „1“ steht für die 1-dimensionale Betrachtung
    + 3: Kugel, der Index „1“ steht für die 1-dimensionale Betrachtung
    + c: in-phase
    + s: out-of-phase
"""

import numpy as np
from Modell.m_single import Kugel, Platte


# Klasse für das Multiple Diffusionsmodell
class m_multiple():
    '''
        delta_c_i: charakteristische in-phase Komponente
        delta_s_i: charakteristische out-of-phase Komponente
        K_i: Frequency Response Parameter
    '''

    def __init__(self, list_delta_c, list_delta_s, list_K):
        self.list_delta_c = list_delta_c
        self.list_delta_s = list_delta_s
        self.list_K = list_K

    def func_chara_c(self):
        chara_c = 0
        n = len(self.list_K)
        for i in range(n):
            chara_c += (self.list_K[i] * self.list_delta_c[i])
        return chara_c

    def func_chara_s(self):
        chara_s = 0
        n = len(self.list_K)
        for i in range(n):
            chara_s += (self.list_K[i] * self.list_delta_s[i])
        return chara_s


if __name__ == '__main__':
    np_f = np.linspace(start=0.001, stop=10, num=10000)
    m_1 = Kugel(a=1e-5, D=1e-11, f=np_f, K=1)
    m_2 = Platte(L=1e-5, D=1e-11, f=np_f, K=1)
    m_multi = m_multiple(list_delta_c=None, list_delta_s=[m_1.func_chara_3s(), m_2.func_chara_1s()], list_K=[1, 1])
    chara_s = m_multi.func_chara_s()
    dict = dict(zip(np_f, chara_s))
    print(dict)
