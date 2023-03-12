"""
Frequency Response Modell
+ Parameter:
    + delta: charakristische Funktion
    + eta: Parameter in der Frequency Response zur Substitution
    + omega: Winkelfrequenz [1/s]
    + f: Frequenz [1/s]
    + t_h, t_R: Zeitkonstante [s]
    + a: Partikelradius [m]
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
import matplotlib.pyplot as plt
from Modell.m_single import Kugel, Platte


# Klasse für das Nichtisothermes Diffusionsmodell
class m_noniso():
    '''
        t_h: Zeitkonstante für den Wärmeaustausch zwischen dem Adsorptionsmittel und seiner Umgebung
        gamma: ein Maß für die Nichtisothermie des Adsorbat-Adsorbens-Systems
        delta_c_i: charakteristische in-phase Komponente
        delta_s_i: charakteristische out-of-phase Komponente
        K_i: Frequency Response Parameter
    '''

    def __init__(self, f, list_delta_c, list_delta_s, list_K, t_h, gamma):
        self.f = f
        self.list_delta_c = list_delta_c
        self.list_delta_s = list_delta_s
        self.list_K = list_K
        self.t_h = t_h
        self.gamma = gamma

    # Gibt die Kreisfrequenz der Volumenmodulation vom Modell
    @property
    def func_omega(self):
        return (2 * np.pi * self.f)

    # Charakteristiche Funktion für das nichtisothermrmes Verhalten(in-phase)
    def func_chara_noniso_c(self):
        ans = 0
        n = len(self.list_K)

        def Zaehler(delta_c, delta_s):
            return ((delta_c * (1 + self.func_omega ** 2 * self.t_h ** 2) + self.gamma * (
                    delta_c ** 2 + delta_s ** 2) * self.func_omega ** 2 * self.t_h ** 2))

        def Nenner(delta_c, delta_s):
            return ((1 + self.gamma * delta_s * self.func_omega * self.t_h) ** 2 + (
                    1 + self.gamma * delta_c) ** 2 * self.func_omega ** 2 * self.t_h ** 2)

        for i in range(n):
            ans += (self.list_K[i] * Zaehler(self.list_delta_c[i], self.list_delta_s[i])
                    / Nenner(self.list_delta_c[i], self.list_delta_s[i]))

        return ans

    # Charakteristiche Funktion für das nichtisothermrmes Verhalten(out-of-phase)
    def func_chara_noniso_s(self):
        ans = 0
        n = len(self.list_K)

        def Nenner(delta_c, delta_s):
            return ((1 + self.gamma * delta_s * self.func_omega * self.t_h) ** 2 + (
                    1 + self.gamma * delta_c) ** 2 * self.func_omega ** 2 * self.t_h ** 2)

        def Zaehler(delta_c, delta_s):
            return (delta_s * (1 + self.func_omega ** 2 * self.t_h ** 2) + self.gamma * (
                    delta_s ** 2 + delta_c ** 2) * self.func_omega * self.t_h)

        for i in range(n):
            ans += (self.list_K[i] * Zaehler(self.list_delta_c[i], self.list_delta_s[i])
                    / Nenner(self.list_delta_c[i], self.list_delta_s[i]))

        return ans


if __name__ == '__main__':
    np_f = np.linspace(start=0.001, stop=10, num=10000)
    m_1 = Kugel(a=1e-5, D=1e-11, f=np_f, K=1)
    m_2 = Platte(L=1e-5, D=1e-11, f=np_f, K=1)
    m_noniso = m_noniso(f=np_f, list_delta_c=[m_1.func_delta_3c(), m_2.func_delta_1c()],
                        list_delta_s=[m_1.func_delta_3s(), m_2.func_delta_1s()], list_K=[1, 1], t_h=1, gamma=0.5)
    chara_s = m_noniso.func_chara_noniso_s()
    #m = m_noniso(f=np_f, list_delta_c=[m_1.func_delta_3c()], list_delta_s=[m_1.func_delta_3s()], list_K=[1], t_h=1, gamma=0.5)
    dict = dict(zip(np_f, chara_s))
    print(dict)