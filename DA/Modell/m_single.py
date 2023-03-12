import numpy as np

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


if __name__ == '__main__':
    np_f = np.linspace(start=0.001, stop=10, num=10000)
    th_in_out_1 = Kugel(a=1e-5, D=1e-11, f=np_f, K=1)
    chara_3c = th_in_out_1.func_chara_3c()
    dict = dict(zip(np_f, chara_3c))
    print(dict)
