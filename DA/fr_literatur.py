import numpy as np
import sys
sys.path.append('../Modell/')

'''
from Modell.m_single import Kugel,Platte
from Modell.m_multiple import m_multiple
from Modell.m_noniso import m_noniso
from Modell.m_hierarchisch import m_hierarchisch
'''
from Modell.modell import *


np_f = np.linspace(start=0.001, stop=10, num=10000)
m_1 = Kugel(a=1e-5, D=1e-11, f=np_f, K=1)
m_2 = Platte(L=1e-5, D=1e-11, f=np_f, K=1)
m_hier = m_hierarchisch(f=np_f, list_delta_c=[m_1.func_chara_3c(), m_2.func_chara_1c()],
                        list_delta_s=[m_1.func_chara_3s(), m_2.func_chara_1s()], list_K=[1, 1], t_R=2.5,
                        K_omega=0.4)
chara_s = m_hier.func_chara_hierar_c()
dict = dict(zip(np_f, chara_s))
print(dict)
