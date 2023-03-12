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


# Klasse von dem Einzeldiffusionsprozessmodell
class m_single():

    # Initialisieren das Modell von der Diffusion in der Kugel
    # enthält 3 Attribute: K, D, f (als Eingabe)
    def __init__(self, D, f, K):
        self.D, self.f, self.K = D, f, K

    # Gibt die Kreisfrequenz der Volumenmodulation vom Modell
    @property
    def func_omega(self):
        return (2 * np.pi * self.f)

    # dB_0_dP_0:  der Anstieg der Isothermen im Gleichgewicht.
    def func_K(self, T_0, V_0, dB_0_dP_0, R=8.314):
        self.K = R * T_0 / V_0 * (dB_0_dP_0)
        return self.K


class Kugel(m_single):
    """
    Geometrie von Kugel:
        a: Kugelradius [m]
        D: Diffusionskoeffizient [m^2/s]
        f: Frequenz [Hz]
    """

    def __init__(self, a, D, f, K):
        super().__init__(D, f, K)
        self.a = a

    # Gibt die Parameter in der Frequency Response zur Substitution von der Kugel
    @property
    def func_eta_3(self):
        return np.sqrt(2 * self.a ** 2 * self.func_omega / self.D)

    # Charakteristische Koeffizient in-phase von der Kugel
    def func_delta_3c(self):
        delta_3c = 3 * (np.sinh(self.func_eta_3) - np.sin(self.func_eta_3)) / (
                self.func_eta_3 * (np.cosh(self.func_eta_3) - np.cos(self.func_eta_3)))
        return delta_3c

    # Charakteristische Koeffizient out-of-phase von der Kugel
    def func_delta_3s(self):
        delta_3s = 3 * (np.sinh(self.func_eta_3) + np.sin(self.func_eta_3)) / (
                self.func_eta_3 * (np.cosh(self.func_eta_3) - np.cos(self.func_eta_3))) - 6 / (self.func_eta_3 ** 2)
        return delta_3s

    # Charakteristische Funktion in-phase von der Kugel
    def func_chara_3c(self):
        return self.K * self.func_delta_3c()

    # Charakteristische Funktion out-of-phase von der Kugel
    def func_chara_3s(self):
        return self.K * self.func_delta_3s()


class Platte(m_single):
    """
        Geometrie von Ebene Platte:
            L: Dicke von der Ebene Platte [m]
            D: Diffusionskoeffizient [m^2/s]
            f: Frequenz [Hz]
    """

    def __init__(self, L, D, f, K):
        super().__init__(D, f, K)
        self.L = L

    # Gibt die Parameter in der Frequency Response zur Substitution von der Ebene Platte
    @property
    def func_eta_1(self):
        return np.sqrt((self.func_omega * self.L ** 2) / (2 * self.D))

    # Charakteristische Koeffizient in-phase von der Ebene Platte
    def func_delta_1c(self):
        delta_1c = (np.sinh(self.func_eta_1) + np.sin(self.func_eta_1)) / (
                self.func_eta_1 * (np.cosh(self.func_eta_1) + np.cos(self.func_eta_1)))
        return delta_1c

    # Charakteristische Koeffizient out-of-phase von der Ebene Platte
    def func_delta_1s(self):
        delta_1s = (np.sinh(self.func_eta_1) - np.sin(self.func_eta_1)) / (
                self.func_eta_1 * (np.cosh(self.func_eta_1) + np.cos(self.func_eta_1)))
        return delta_1s

    # Charakteristische Funktion in-phase von der Ebene Platte
    def func_chara_1c(self):
        return self.K * self.func_delta_1c()

    # Charakteristische Funktion out-of-phase von der Ebene Platte
    def func_chara_1s(self):
        return self.K * self.func_delta_1s()


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


class m_hierarchisch():
    '''
        t_R: Zeitkonstante für den Massenaustausch zwischen den beiden Porensorten
        K_omega: relative Gleichgewichtskonstante der Konzentration in den Speicherporen im Verhältnis zu den Transportporen
        delta_c_i: charakteristische in-phase Komponente
        delta_s_i: charakteristische out-of-phase Komponente
        K_i: Frequency Response Parameter
    '''

    def __init__(self, f, list_delta_c, list_delta_s, list_K, t_R, K_omega):
        self.f = f
        self.list_delta_c = list_delta_c
        self.list_delta_s = list_delta_s
        self.list_K = list_K
        self.t_R = t_R
        self.K_omega = K_omega

    # Gibt die Kreisfrequenz der Volumenmodulation vom Modell
    @property
    def func_omega(self):
        return (2 * np.pi * self.f)

    # Charakteristiche Funktion für das hierarchischen Porensystemen(in-phase)
    def func_chara_hierar_c(self):
        ans = 0
        n = len(self.list_K)

        def func_i(delta_c, delta_s):
            chara_c_i = (delta_c + self.K_omega * ((delta_c - delta_s * self.func_omega * self.t_R)
                                                   / (1 + self.func_omega ** 2 * self.t_R ** 2)))
            return chara_c_i

        for i in range(n):
            ans += self.list_K[i] * func_i(self.list_delta_c[i], self.list_delta_s[i])
            return ans

    # Charakteristiche Funktion für das hierarchischen Porensystemen(out-of-phase)
    def func_chara_hierar_s(self):
        ans = 0
        n = len(self.list_K)

        def func_i(delta_c, delta_s):
            chara_s_i = (delta_s + self.K_omega * ((delta_s + delta_c * self.func_omega * self.t_R)
                                                   / (1 + self.func_omega ** 2 * self.t_R ** 2)))
            return chara_s_i

        for i in range(n):
            ans += self.list_K[i] * func_i(self.list_delta_c[i], self.list_delta_s[i])
        return ans


if __name__ == '__main__':
    np_f = np.linspace(start=0.001, stop=10, num=10000)
    m_1 = Kugel(a=1e-5, D=1e-11, f=np_f, K=1)
    m_2 = Platte(L=1e-5, D=1e-11, f=np_f, K=1)
    m_hier = m_hierarchisch(f=np_f, list_delta_c=[m_1.func_chara_3c(), m_2.func_chara_1c()],
                            list_delta_s=[m_1.func_chara_3s(), m_2.func_chara_1s()], list_K=[1, 1], t_R=2.5,
                            K_omega=0.4)
    chara_s = m_hier.func_chara_hierar_c()
    dict = dict(zip(np_f, chara_s))
    print(dict)
