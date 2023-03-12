import timeit
import numpy as np

from scipy.optimize import minimize


from prop_setting import proplib
#print(proplib)
if proplib.upper() == 'TREND':
    from prop import Trend_props_E_ext as props_E
    from prop import Trend_props_lMeta as props_lMeta
    from prop import Trend_p_min as p_min
elif proplib.upper() == 'COOLPROP':
    from prop import CP_props_E as props_E
    from prop import CP_p_min as p_min
 
class crit_flow_model():
    """
    input
        p_0: stagnation pressure ~ static pressure and
        h_0: stagnation enthalpy or
        x_0: stagnation quality        
    sub
        t: throat, 0: stagnation
        g: gas phase
        l: liquid phase
    """
    def __init__(self, fluid, p_0, h_0 = None, x_0 = None):
        self.fluid = fluid
        self.p_0 = p_0
        self.Gc = 0
        self.p_t = 0
        if h_0 is None and x_0 is None:
            self.x_0 = 0
            self.props_0 = self.props0() 
        elif h_0 is not None:
            self.h_0 = h_0
            self.props_0 = self.props0()  
            self.s_0 = self.props_0['s']
            self.props_g0 = self.props_g(self.p_0)
            self.props_l0 = self.props_l(self.p_0)
            h_l0 = self.props_l0['h']
            h_g0 = self.props_g0['h']
            self.x_0 = (self.h_0 - h_l0) / (h_g0 - h_l0)
        elif x_0 is not None:
            self.x_0 = x_0             
            self.props_g0 = self.props_g(self.p_0)
            self.props_l0 = self.props_l(self.p_0)
            h_l0 = self.props_l0['h']
            h_g0 = self.props_g0['h']
            self.h_0 = h_l0 + self.x_0 * (h_g0 - h_l0)
            self.props_0 = self.props0()
            self.s_0 = self.props_0['s']
        else:
            raise ValueError('input either h_0 or x_0') 
            
        self.p_sat_T0 = self.p_sat(self.props_0['T'])
     
    
    
    def props0(self):
        props_0 = props_E(self.fluid, self.p_0, self.h_0, 
                                in1 = 'p', in2 = 'h',
                                inkl = ['T', 'd', 'h', 's']
                                ) 
        
        return props_0    

    def p_sat(self, T): 
        props_sat = props_E(self.fluid, T, 0, 
                              	in1 = 'T', in2 = 'Q',
                                inkl = ['p']) 
        return props_sat['p']    

    def props_g(self, p):
        return props_E(self.fluid, p, 1, in1 = 'p', in2 = 'Q',
                                inkl = ['cp', 'cv', 'd', 'h', 's'])

    def props_1(self, p):
        return props_E(self.fluid, p, 0, in1 = 'p', in2 = 'Q',
                                inkl = ['cp', 'd', 'h', 's'])    
    
    def check_props(self):
        props_ls  = ['T', 'd', 'h', 's']
        for prop in props_ls:
            try: 
                self.props_0[prop]
            except AttributeError:
                raise ValueError()

        
class Henry_Fauske(crit_flow_model):        
    """
    HENRY-FAUSKE MODEL (F-S-MODEL) 1971
    R.HENRY and H.Fauske
    The two-phase critical flow of one-component mixtures in nozzles, 
    orifices, and short tubes. ASME J. Heat Transfer, 179–187  
        sub: E:  equilibrium; NE: Non-equilibrium
    """
    
    def __init__(self, fluid, p_0, h_0 = None, x_0 = None, 
                 C_NE = None, run = True, 
                 **solver_par):
        
        try:
            super().__init__(fluid, p_0, h_0 = h_0, x_0 = x_0)
        except ValueError:
            return
        
        if C_NE is None:
            C_NE = 0.14 
            # partial phase change rate, 
            # 0.14 by default, given by Henry&Fauske 1971, 
            # according to exp.res. from Starkman et al. 1964
        self.C_NE = C_NE
        #self.p_t = 3.9e5
        self.solver_par = self.__default_solver_par()
        if solver_par:
            self.solver_par.update(solver_par)
        self.rel_dev = 1
        if run is True:
            self.main()
        
    @property
    def props_gtE(self):
        return props_E(self.fluid, self.p_t, 1, in1 = 'p', in2 = 'Q',
                            inkl = ['cp', 'cv', 'd', 's'])
    @property
    def props_ltE(self):
        return props_E(self.fluid, self.p_t, 0, in1 = 'p', in2 = 'Q',
                                inkl = ['cp', 'Ds/Dp|E', 's'])     

    def Kappa(self, p):
        return props_E(self.fluid, self.p_0, 1, in1 = 'p', in2 = 'Q',
                                inkl = ['ks'])['ks']

    @property
    def n(self):
        #### thermal equilibrium polytropic exponent (tangren et al. 1949)
        cp_l = self.props_ltE['cp']
        cp_g = self.props_gtE['cp']
        Kappa = self.Kappa(self.p_t)
        #cv_g = self.props_g0['cv']
        return ((1 - self.x_0) * cp_l / cp_g + 1) / \
                ((1 - self.x_0) * cp_l / cp_g + 1 / Kappa ) 
    
    @property
    def __max_dev(self):  
        return 0.001
    
    def N_HF(self, x_tE, C_NE):
        """
        Thermal nonequilibreium factor
        experimental parameter for partial phase change rate
        definition:
            dx/dp = N * (dx/dp)Equilibrium
            applied only to the liquid phase
        """
        if x_tE < C_NE:
                N = x_tE / C_NE
        else:
            N = 1  
        return N
    
    def v_gt_is(self, v_g0, p_t, p_0, Kappa):
        return v_g0 * (p_t / p_0) ** (-1 / Kappa)
            
    def GC2_2PH(self, p_t):
        # @stagnation 0
        p_0 = self.p_0
        x_0 = self.x_0
        
        v_g0 = 1 / self.props_g0['d']
        s_g0 = self.props_g0['s']
        
        v_l0 = 1 / self.props_l0['d']
        s_l0 = self.props_l0['s']
        s_0 = self.props_0['s']

        # 0 -> t
        self.p_t = p_t
        Kappa = self.Kappa(self.p_0) # isentropic exponent        
        x_t = x_0 # assume: no mass transfer during expansion  
        s_lt = s_l0
        s_gt = s_g0 # assume: each phase expands isentropically
        v_gt = self.v_gt_is(v_g0, p_t, p_0, Kappa)
        #v_gt = 1/self.props_gtE['d']
        v_lt = v_l0 # assume: liquid phase incompressible
               
        # @throat t 
        ## properties from equilibrium state E
        s_gtE = self.props_gtE['s']
        s_ltE = self.props_ltE['s']
        Ds_l__Dp__tE = self.props_ltE['Ds/Dp|E']
        #v_gtE = 1 / self.props_gtE['d']
        #v_ltE = v_l0
        s_tE = s_0
        x_tE = (s_tE - s_ltE) / (s_gtE - s_ltE)          
        
        ## correlation for non equilibrium state NE        
        ### gas phase follows the polytropic process at throat
        n = self.n
        Dv_g__Dp__t = - v_gt / (n * p_t) 
        cp_gt = self.props_gtE['cp']
        Ds_g__Dp__t = - cp_gt / p_t * (1 / n - 1 / Kappa)
        
        ### liquid phase incompressible (partial phase change rate N)
        Dv_l__Dp__t = 0              
        Ds_l__Dp__t = self.N_HF(x_tE, self.C_NE) * Ds_l__Dp__tE \
                        * (s_gt - s_lt) / (s_gtE - s_ltE)
        
        Dx__DP__t = - (x_t * Ds_g__Dp__t + (1 - x_t) * Ds_l__Dp__t) / (s_gt - s_lt)
        
        Gc2 = - 1 / (x_t * Dv_g__Dp__t + (1 - x_t) * Dv_l__Dp__t + (v_gt - v_lt) * Dx__DP__t)
        return Gc2
    
    def G2Mom_2PH(self, p_t):
        p_0 = self.p_0
        x_0 = self.x_0
        Kappa = self.Kappa(self.p_0)
        
        v_g0 = 1 / self.props_g0['d']
        
        v_l0 = 1 / self.props_l0['d']
        
        #p_t = self.p_t
        x_t = x_0
        v_gt = self.v_gt_is(v_g0, p_t, p_0, Kappa)
        #v_gt = 1/self.props_gtE['d']
        v_lt = v_l0 # assume: liquid phase incompressible
        v_t = x_t * v_gt + (1 - x_t) * v_lt
        
        G2Mom = 2 * (x_0 * Kappa / (Kappa - 1) * (p_0 * v_g0 - p_t * v_gt) \
                + (1 - x_0) * v_l0 * (p_0 - p_t)) / v_t ** 2  
           
        return G2Mom
    
    def GC2_SC(self, p_t):
        """
        assumption:
        x_t = x_0 = 0
        dx/dp_t != 0
        """
        # @stagnation 0
        p_0 = self.p_0
        v_l0 = 1 / self.props_0['d']
        s_0 = self.props_0['s']
        
        # @throat t
        self.p_t = p_t
        v_gtE = 1 / self.props_gtE['d']
        v_gt = v_gtE # assume: saturated vapor is produced 
        v_lt = v_l0
        
        s_tE = s_0         
        s_gtE = self.props_gtE['s']
        s_ltE = self.props_ltE['s']
        x_tE = (s_tE - s_ltE) / (s_gtE - s_ltE) 
        Ds_l__Dp__tE = self.props_ltE['Ds/Dp|E']
        Dx__DP__t = - self.N_HF(x_tE, self.C_NE) * Ds_l__Dp__tE \
                        / (s_gtE - s_ltE)
        Gc2 = - 1 / ((v_gt - v_lt) * Dx__DP__t)
        return Gc2
    
    def G2Mom_SC(self, p_t):
        p_0 = self.p_0   
        v_l0 = 1 / self.props_0['d']
        G2Mom = 2 * (p_0 - p_t) / v_l0 
        return G2Mom
    
    def G2Mom(self, p_t):
        if self.x_0 <= 0: # stagnation point: subcooled
            return self.G2Mom_SC(p_t)
        else: # stagnation point: two-phase
            return self.G2Mom_2PH(p_t)
        
    def solve_GC2PH(self):  
        
        p_tmin = p_min(self.fluid)
        p_tmax = self.p_0 * 0.99
        k = 0.87
        while True:
            p_t_guess = k * p_tmax  
            if p_tmax < p_tmin:
                #trend
                break
            res = minimize(self.__funcGc_zero_2PH, 
                           p_t_guess,
                           #method = self.solver_par['method'],
                           bounds = ((p_tmin, p_tmax),),
                           **self.solver_par
                           
            )
            self.p_t = res.x[0]
            #print(self.G2Mom_2PH(self.p_t) ** 0.5)
            #print(self.GC2_2PH(self.p_t) ** 0.5)
            self.Gc = self.G2Mom_2PH(self.p_t) ** 0.5
            self.rel_dev = (self.__funcGc_zero_2PH(self.p_t)) ** 0.5 \
                                    / self.Gc
            if self.__valid_results():
                break
            else:
                k = k * 0.95
                if k < 0.6:
                    print('result not convergent')
                    break
     
    def solve_GCSC(self):
          
        p_tmin = p_min(self.fluid)
        p_tmax = self.p_sat_T0
        
        p_t_guess = 0.87 * p_tmax
        """
        if p_t_guess < p_tmin:
            #trend
            return
        """
        res = minimize(self.__funcGc_zero_SC, 
                       p_t_guess,
                       #method = self.solver_par['method'],
                       bounds = ((p_tmin, p_tmax),),
                       **self.solver_par
        )
        self.p_t = res.x[0]
        #print(self.G2Mom_SC(self.p_t) ** 0.5)
        #print(self.GC2_SC(self.p_t) ** 0.5)
        self.Gc = self.G2Mom_SC(self.p_t) ** 0.5
        self.rel_dev = (self.__funcGc_zero_SC(self.p_t)) ** 0.5 \
                                / self.Gc
        if self.__valid_results():
            pass
        else:
            print('result not convergent')
    
    def main(self):          
        #print(self.solver_par)
        try:
            self.check_props()
        except ValueError:
            return
        if self.x_0 <= 0: # stagnation point: subcooled
            try:
                self.solve_GCSC()
            except ValueError:
                return
        else: # stagnation point: two-phase
            try:
                self.solve_GC2PH()
            except ValueError:
                return
    
    def __funcGc_zero_2PH(self, p_t):
        return np.abs(self.GC2_2PH(p_t) - self.G2Mom_2PH(p_t))   
    
    def __funcGc_zero_SC(self, p_t):
        return np.abs(self.GC2_SC(p_t) - self.G2Mom_SC(p_t))
         
    def __valid_results(self):
        if self.rel_dev > self.__max_dev:            
            return False
        else:
            return True    
        
    def __default_solver_par(self):
        par = {
            'method' : 'trust-constr',               
        }
        return par

"""
class Henry_Fauske_modified(crit_flow_model):        
    
    def __init__(self, fluid, p_0, h_0 = None, x_0 = None, 
                 C_NE = None, 
                 **solver_par):
        super().__init__(fluid, p_0, h_0 = h_0, x_0 = x_0)
        if C_NE is None:
            C_NE = 0.14 
            # partial phase change rate, 
            # 0.14 by default, given by Henry&Fauske 1971, 
            # according to exp.res. from Starkman et al. 1964
        self.C_NE = C_NE
        #self.p_t = 3.9e5
        self.solver_par = self.__default_solver_par()
        if solver_par:
            self.solver_par.update(solver_par)
        self.rel_dev = 1
        self.main()
            
    @property
    def props_gtE(self):
        props_gtE = props_E(self.fluid, self.p_t, 1, in1 = 'p', in2 = 'Q',
                            inkl = ['cp', 'cv', 'd', 's'])
        return props_gtE
    @property
    def props_ltE(self):

        props_ltE = props_E(self.fluid, self.p_t, 0, in1 = 'p', in2 = 'Q',
                                inkl = ['cp', 'Ds/Dp|E', 's']) 
        return props_ltE
    @property 
    def Kappa(self):
        cp_g = self.props_gtE['cp']
        cv_g = self.props_gtE['cv']
        return cp_g / cv_g
    
    @property
    def n(self):
        #### thermal equilibrium polytropic exponent (tangren et al. 1949)
        cp_l = self.props_ltE['cp']
        cp_g = self.props_gtE['cp']
        #cv_g = self.props_g0['cv']
        return ((1 - self.x_0) * cp_l / cp_g + 1) / \
                ((1 - self.x_0) * cp_l / cp_g + 1 / self.Kappa ) 
    
    def N_HF(self, x_tE, C_NE):
        if x_tE < C_NE:
            N = x_tE / C_NE
        else:
            N = 1  
        return N
    
    def v_gt_is(self, v_g0, p_t, p_0, Kappa):
        return v_g0 * (p_t / p_0) ** (-1 / Kappa)
    
    def GC2_2PH(self, p_t):
        # @stagnation 0
        p_0 = self.p_0
        x_0 = self.x_0
        
        cp_g = self.props_g0['cp']
        cv_g = self.props_g0['cv']
        v_g0 = 1 / self.props_g0['d']
        s_g0 = self.props_g0['s']
        
        cp_l = self.props_l0['cp']
        v_l0 = 1 / self.props_l0['d']
        s_l0 = self.props_l0['s']
        s_0 = self.props_0['s']

        # 0 -> t
        self.p_t = p_t
        Kappa = self.Kappa # isentropic exponent        
        x_t = x_0 # assume: no mass transfer during expansion  
        s_lt = s_l0
        s_gt = s_g0 # assume: each phase expands isentropically

        v_gtE = 1/self.props_gtE['d'] 
        v_gt = v_gtE # ***modified*** saturated gas
        v_lt = v_l0 # assume: liquid phase incompressible
               
        # @throat t 
        ## properties from equilibrium state E
        s_gtE = self.props_gtE['s']
        s_ltE = self.props_ltE['s']
        Ds_l__Dp__tE = self.props_ltE['Ds/Dp|E']

        #v_ltE = v_l0
        s_tE = s_0
        x_tE = (s_tE - s_ltE) / (s_gtE - s_ltE)          
        
        ## correlation for non equilibrium state NE        
        ### gas phase follows the polytropic process at throat
        n = self.n
        Dv_g__Dp__t = - v_gt / (n * p_t) 
        cp_gt = self.props_gtE['cp']
        Ds_g__Dp__t = - cp_gt / p_t * (1 / n - 1 / Kappa)
        
        ### liquid phase incompressible (partial phase change rate N)
        Dv_l__Dp__t = 0              
        Ds_l__Dp__t = self.N_HF(x_tE, self.C_NE) * Ds_l__Dp__tE \
                        * (s_gt - s_lt) / (s_gtE - s_ltE)
        
        Dx__DP__t = - (x_t * Ds_g__Dp__t + (1 - x_t) * Ds_l__Dp__t) / (s_gt - s_lt)
        
        Gc2 = - 1 / (x_t * Dv_g__Dp__t + (1 - x_t) * Dv_l__Dp__t + (v_gt - v_lt) * Dx__DP__t)
        return Gc2
    
    def G2Mom_2PH(self, p_t):
        p_0 = self.p_0
        x_0 = self.x_0
        Kappa = self.Kappa
        v_g0 = 1 / self.props_g0['d']
        v_l0 = 1 / self.props_l0['d']
        
        self.p_t = p_t
        x_t = x_0
        v_gt = self.v_gt_is( v_g0, p_t, p_0, Kappa) # ***modified*** saturated gas
        v_lt = v_l0 # assume: liquid phase incompressible
        v_t = x_t * v_gt + (1 - x_t) * v_lt
                
        G2Mom = 2 * (x_0 * Kappa / (Kappa - 1) * (p_0 * v_g0 - p_t * v_gt) \
                + (1 - x_0) * v_l0 * (p_0 - p_t)) / v_t ** 2  
           
        return G2Mom
    
    def GC2_SC(self, p_t):

        # @stagnation 0
        v_l0 = 1 / self.props_0['d']
        s_0 = self.props_0['s']
        
        # @throat t
        self.p_t = p_t
        v_gtE = 1 / self.props_gtE['d']
        v_gt = v_gtE # assume: saturated vapor is produced
        v_lt = v_l0
        
        s_tE = s_0
        s_gtE = self.props_gtE['s']
        s_ltE = self.props_ltE['s']
        x_tE = (s_tE - s_ltE) / (s_gtE - s_ltE)
        Ds_l__Dp__tE = self.props_ltE['Ds/Dp|E']
        Dx__DP__t = - self.N_HF(x_tE, self.C_NE) * Ds_l__Dp__tE \
                        / (s_gtE - s_ltE)
        Gc2 = - 1 / ((v_gt - v_lt) * Dx__DP__t)
        return Gc2
    
    def G2Mom_SC(self, p_t):
        p_0 = self.p_0   
        v_l0 = 1 / self.props_0['d']
        G2Mom = 2 * (p_0 - p_t) / v_l0 
        return G2Mom
    
    def solve_GC2PH(self):  
        
        p_tmin = p_min(fluid)
        p_tmax = self.p_0 * 0.99
        k = 0.87
        while True:
            p_t_guess = k * p_tmax  
            if p_tmax < p_tmin:
                #trend
                break
            res = minimize(self.__funcGc_zero_2PH, 
                           p_t_guess,
                           method = self.solver_par['method'],
                           bounds = ((p_tmin, p_tmax),)
            )
            self.p_t = res.x[0]
            #print(self.G2Mom_2PH(self.p_t) ** 0.5)
            #print(self.GC2_2PH(self.p_t) ** 0.5)
            self.Gc = self.G2Mom_2PH(self.p_t) ** 0.5
            self.rel_dev = (self.__funcGc_zero_2PH(self.p_t)) ** 0.5 \
                                    / self.Gc
            if self.__valid_results():
                break
            else:
                k = k * 0.95
                if k < 0.6:
                    print('result not convergent')
                    break
     
    def solve_GCSC(self):
        p_tmin = p_min(fluid)
        p_tmax = self.p_sat_T0
        p_t_guess = 0.9 * p_tmax
        if p_t_guess < p_tmin:
            #trend
            return
        res = minimize(self.__funcGc_zero_SC, 
                       p_t_guess,
                       method = self.solver_par['method'],
                       bounds = ((p_tmin, p_tmax),)
        )
        self.p_t = res.x[0]
        #print(self.G2Mom_SC(self.p_t) ** 0.5)
        #print(self.GC2_SC(self.p_t) ** 0.5)
        self.Gc = self.G2Mom_SC(self.p_t) ** 0.5
        self.rel_dev = (self.__funcGc_zero_SC(self.p_t)) ** 0.5 \
                                / self.Gc
        if self.__valid_results():
            pass
        else:
            print('result not convergent')
    
    def main(self):          
        #print(self.solver_par)
        
        if self.x_0 <= 0: # stagnation point: subcooled
            try:
                self.solve_GCSC()
            except:
                return
        else: # stagnation point: two-phase
            try:
                self.solve_GC2PH()
            except:
                return
    
    def __funcGc_zero_2PH(self, p_t):
        #self.p_t = p_t
        return np.abs(self.GC2_2PH(p_t) - self.G2Mom_2PH(p_t))   
    
    def __funcGc_zero_SC(self, p_t):
        #self.p_t = p_t
        return np.abs(self.GC2_SC(p_t) - self.G2Mom_SC(p_t))
         
    def __valid_results(self):
        if self.rel_dev > self.solver_par['max_dev']:            
            return False
        else:
            return True
        
    def __default_solver_par(self):
        par = {
            'method' : 'trust-constr',
            'max_dev' : 0.001                
        }
        return par   
"""
class Henry_Fauske_modified(Henry_Fauske):        
    """
    modified according to Chung et al. 2010 (Mars Code Manual Volume V)
    """
    def __init__(self, fluid, p_0, h_0 = None, x_0 = None, 
                 C_NE = None, 
                 **solver_par):
        try:
            super().__init__(fluid, p_0, h_0 = h_0, x_0 = x_0, C_NE = C_NE, 
                 run = True,
                 **solver_par)
        except ValueError:
            return
        
    @property
    def props_gtE(self):        
        return super().props_gtE
    
    @property
    def props_ltE(self):
        return super().props_ltE

    @property
    def n(self):
       return super().n
    
    def GC2_2PH(self, p_t):
        # @stagnation 0
        p_0 = self.p_0
        x_0 = self.x_0
        
        cp_g = self.props_g0['cp']
        cv_g = self.props_g0['cv']
        v_g0 = 1 / self.props_g0['d']
        s_g0 = self.props_g0['s']
        
        cp_l = self.props_l0['cp']
        v_l0 = 1 / self.props_l0['d']
        s_l0 = self.props_l0['s']
        s_0 = self.props_0['s']

        # 0 -> t
        self.p_t = p_t
        Kappa = (self.Kappa(self.p_0) + self.Kappa(p_t)) / 2 # isentropic exponent        
        x_t = x_0 # assume: no mass transfer during expansion  
        s_lt = s_l0
        s_gt = s_g0 # assume: each phase expands isentropically

        v_gtE = 1/self.props_gtE['d'] 
        v_gt = v_gtE # ***modified*** saturated gas
        v_lt = v_l0 # assume: liquid phase incompressible
               
        # @throat t 
        ## properties from equilibrium state E
        s_gtE = self.props_gtE['s']
        s_ltE = self.props_ltE['s']
        Ds_l__Dp__tE = self.props_ltE['Ds/Dp|E']

        #v_ltE = v_l0
        s_tE = s_0
        x_tE = (s_tE - s_ltE) / (s_gtE - s_ltE)          
        
        ## correlation for non equilibrium state NE        
        ### gas phase follows the polytropic process at throat
        n = self.n
        Dv_g__Dp__t = - v_gt / (n * p_t) 
        cp_gt = self.props_gtE['cp']
        Ds_g__Dp__t = - cp_gt / p_t * (1 / n - 1 / Kappa)
        
        ### liquid phase incompressible (partial phase change rate N)
        Dv_l__Dp__t = 0              
        Ds_l__Dp__t = self.N_HF(x_tE, self.C_NE) * Ds_l__Dp__tE \
                        * (s_gt - s_lt) / (s_gtE - s_ltE)
        
        Dx__DP__t = - (x_t * Ds_g__Dp__t + (1 - x_t) * Ds_l__Dp__t) / (s_gt - s_lt)
        
        Gc2 = - 1 / (x_t * Dv_g__Dp__t + (1 - x_t) * Dv_l__Dp__t + (v_gt - v_lt) * Dx__DP__t)
        return Gc2
         
class Henry_Fauske_modified2(Henry_Fauske):        

    def __init__(self, fluid, p_0, h_0 = None, x_0 = None, 
                 C_NE = None,
                 **solver_par):        
        try:
            super().__init__(fluid, p_0, h_0 = h_0, x_0 = x_0, C_NE = C_NE, 
                  run = False,
                 **solver_par)
            
        except ValueError:
            return           
        self.main()
        
    @property
    def props_gtE(self):        
        return super().props_gtE
    
    @property
    def props_ltE(self):
        return props_E(self.fluid, self.p_t, 0, in1 = 'p', in2 = 'Q',
                                inkl = ['cp', 'Ds/Dp|E', 'd', 's', 'betas'])   
    
    @property
    def n(self):
       return super().n    

    def Kappa_l0(self):
        if self.x_0 < 0:
            return props_E(self.fluid, self.p_0, self.s_0, in1 = 'p', in2 = 's',
                                inkl = ['ks'])['ks']
        else:
            return props_E(self.fluid, self.p_0, 0, in1 = 'p', in2 = 'Q',
                                inkl = ['ks'])['ks']
        
    """
    def props_lMeta(self, p, s):
        return props_lMeta(self.fluid, p, s, in1 = 'p', in2 = 's', 
                     inkl = ['cp', 'cv' , 'd'])

    def integ_vdp_s_lMeta(self, p1, p2, s, step = 10):
        dp = (p2 - p1) / step    
        p = np.arange(p1, p2, dp)        
        v = np.zeros(step-1)
        gam = np.zeros(step-1)
        vdp = np.zeros(step-1)
        for i in np.arange(step - 1):
            p_m = (p[i] + p[i + 1]) / 2
            gam[i] = self.props_lMeta(p_m, s)['cp'] / \
                self.props_lMeta(p_m, s)['cv']
            try:
                v[i] = 1 / self.props_lMeta(p_m, s)['d']
            except ZeroDivisionError:
                v[i] = v[i-1]
            #pdV[i] = gam[i] / (gam[i] - 1) * (p[i+1] * v[i+1] - p[i] * v[i])
            vdp[i] = v[i] * dp
        res = sum(vdp)
        return res
    """
    def GC2_2PH(self, p_t):
        # @stagnation 0
        p_0 = self.p_0
        x_0 = self.x_0
        
        #v_g0 = 1 / self.props_g0['d']
        s_g0 = self.props_g0['s']
        
        v_l0 = 1 / self.props_l0['d']
        s_l0 = self.props_l0['s']
        s_0 = self.props_0['s']

        # 0 -> t
        self.p_t = p_t
        x_t = x_0 # assume: no mass transfer during expansion  
        s_lt = s_l0
        s_gt = s_g0 # assume: each phase expands isentropically

        v_gtE = 1/self.props_gtE['d'] 
        #v_ltE = 1/self.props_ltE['d']
        v_gt = v_gtE # ***modified*** saturated gas
        Kappa_l = self.Kappa_l0()
        v_lt = self.v_gt_is(v_l0, p_t, p_0, Kappa_l)
        #print(v_l0, v_lt,self.v_gt_is(v_l0, p_t, p_0, self.Kappa_l))               
        # @throat t 
        ## properties from equilibrium state E
        Kappa_t =  self.Kappa(p_t)
        s_gtE = self.props_gtE['s']
        s_ltE = self.props_ltE['s']
        Ds_l__Dp__tE = self.props_ltE['Ds/Dp|E']

        #v_ltE = v_l0
        s_tE = s_0
        x_tE = (s_tE - s_ltE) / (s_gtE - s_ltE)          
        
        ## correlation for non equilibrium state NE        
        ### gas phase follows the polytropic process at throat
        n = self.n
        Dv_g__Dp__t = - v_gt / (n * p_t) 
        cp_gt = self.props_gtE['cp']
        Ds_g__Dp__t = - cp_gt / p_t * (1 / n - 1 / Kappa_t)
        

        Dv_l__Dp__t = - self.props_ltE['betas'] * v_lt  # liquid compressibilty refers to saturated point      
        Ds_l__Dp__t = self.N_HF(x_tE, self.C_NE) * Ds_l__Dp__tE \
                        * (s_gt - s_lt) / (s_gtE - s_ltE)
        
        Dx__DP__t = - (x_t * Ds_g__Dp__t + (1 - x_t) * Ds_l__Dp__t) / (s_gt - s_lt)
        #print(x_t * Dv_g__Dp__t, (1 - x_t) * Dv_l__Dp__t, (v_gt - v_lt) * Dx__DP__t)
        Gc2 = - 1 / (x_t * Dv_g__Dp__t + (1 - x_t) * Dv_l__Dp__t + (v_gt - v_lt) * Dx__DP__t)
        return Gc2
    
    def G2Mom_2PH(self, p_t):
        p_0 = self.p_0
        x_0 = self.x_0
        s_0 = self.s_0
        Kappa = (self.Kappa(self.p_0) + self.Kappa(p_t)) / 2
        Kappa_l = self.Kappa_l0()
        
        v_g0 = 1 / self.props_g0['d']        
        v_l0 = 1 / self.props_l0['d']
        
        #p_t = self.p_t
        x_t = x_0
        v_gt = self.v_gt_is(v_g0, p_t, p_0, Kappa)
        #v_gt = 1/self.props_gtE['d']

        v_lt = self.v_gt_is(v_l0, p_t, p_0, Kappa_l) 
        v_t = x_t * v_gt + (1 - x_t) * v_lt
        integ_vdp_gt0 = Kappa / (Kappa - 1) * (p_0 * v_g0 - p_t * v_gt)
        integ_vdp_lt0 = Kappa_l / (Kappa_l - 1) * (p_0 * v_l0 - p_t * v_lt)
        #print(integ_vdp_gt0,integ_vdp_lt0, v_l0 * (p_0 - p_t))
        G2Mom = 2 * (x_0 * integ_vdp_gt0 + (1 - x_0) * integ_vdp_lt0) / v_t ** 2  
        print(G2Mom ** 0.5  * v_t)
        return G2Mom
    
    def GC2_SC(self, p_t):
        """
        assumption:
        x_t = x_0 = 0
        dx/dp_t != 0
        """
        # @stagnation 0
        p_0 = self.p_0
        s_0 = self.props_0['s']
        v_l0 = 1 / self.props_0['d']
        # @throat t
        self.p_t = p_t
        v_gtE = 1 / self.props_gtE['d']
        v_gt = v_gtE # assume: saturated vapor is produced 
        
        Kappa_l = self.Kappa_l0()
        v_lt = self.v_gt_is(v_l0, p_t, p_0, Kappa_l)     
        s_tE = s_0         
        s_gtE = self.props_gtE['s']
        s_ltE = self.props_ltE['s']
        x_tE = (s_tE - s_ltE) / (s_gtE - s_ltE) 
        Dv_l__Dp__t = - self.props_ltE['betas'] * v_lt  # liquid compressibilty refers to saturated point 
        Ds_l__Dp__tE = self.props_ltE['Ds/Dp|E']
        Dx__DP__t = - self.N_HF(x_tE, self.C_NE) * Ds_l__Dp__tE \
                        / (s_gtE - s_ltE)
        Gc2 = - 1 / (Dv_l__Dp__t + (v_gt - v_lt) * Dx__DP__t)
        return Gc2
    
    def G2Mom_SC(self, p_t):
        p_0 = self.p_0   
        s_0 = self.s_0
        Kappa_l = self.Kappa_l0()        
        v_l0 = 1 / self.props_0['d']
        v_lt = self.v_gt_is(v_l0, p_t, p_0, Kappa_l) 
        integ_vdp_lt0 = Kappa_l / (Kappa_l - 1) * (p_0 * v_l0 - p_t * v_lt)
        G2Mom = 2 * (integ_vdp_lt0) / v_lt ** 2
        #print(G2Mom ** 0.5  * v_lt)
        return G2Mom     
    
class HEM(crit_flow_model):
    """
    homogenous equilibrium model
    
    """
    def __init__(self, fluid, p_0, h_0 = None, x_0 = None,
                 **solver_par):

        super().__init__(fluid, p_0, h_0 = h_0, x_0 = x_0)
        self.solver_par = self.__default_solver_par()
        if solver_par:
            self.solver_par.update(solver_par)
        self.main()

    def rho_tE(self):
        props_tE = props_E(self.fluid, self.p_t, self.s_0, in1 = 'p', in2 = 's',
                                inkl = ['d'])['d']
        return props_tE

    @property
    def props_gtE(self):

        props_gtE = props_E(self.fluid, self.p_t, 1, in1 = 'p', in2 = 'Q',
                                inkl = ['d', 'h', 's'])
        return props_gtE
    @property
    def props_ltE(self):
        props_ltE = props_E(self.fluid, self.p_t, 0, in1 = 'p', in2 = 'Q',
                                inkl = ['d', 'h', 's']) 
        return props_ltE

    @property
    def x_tE(self):
        s_0 = self.props_0['s']
        s_tE = s_0 # assume isentropic expension
        s_gtE = self.props_gtE['s']
        s_ltE = self.props_ltE['s']
        if s_tE > s_ltE:
            x_tE = (s_tE - s_ltE) / (s_gtE - s_ltE)
        else:
            x_tE = 0
        return x_tE
    
    def G2Nrg(self, p_t):
        h_0 = self.h_0
        self.p_t = p_t
        #rho_gtE = self.props_gtE['d']
        #rho_ltE = self.props_ltE['d']        

        x_tE = self.x_tE
        
        h_gtE = self.props_gtE['h']
        h_ltE = self.props_ltE['h']
        #rho_tHE = 1 / (x_tE / rho_gtE + (1 - x_tE) / rho_ltE)
        rho_tHE = self.rho_tE()
        h_tHE = x_tE * h_gtE + (1 - x_tE) * h_ltE
        
        G2 = 2 * (h_0 - h_tHE) * rho_tHE ** 2        
        return G2
        
    def solve_GC(self): 
        if self.x_0 < 0:
            p_tmin = p_min(self.fluid)
            p_tmax = self.p_sat_T0
        else:
            p_tmin = p_min(self.fluid)
            p_tmax = self.p_0 * 0.99

        k = 0.87
        p_t_guess = k * p_tmax  
        if p_tmax < p_tmin:
            #trend
            return
        res = minimize(self.__funcGc_min, 
                       p_t_guess,
                       method = self.solver_par['method'],
                       bounds = ((p_tmin, p_tmax),),
                       tol = 1e-4
        )
        self.p_t = res.x[0]
        self.Gc = self.G2Nrg(self.p_t) ** 0.5
    
    def main(self):
        try:
            self.solve_GC()
        except ValueError:
            return
        
    def __funcGc_min(self, p_t): 
        return -self.G2Nrg(p_t)
    
    def __default_solver_par(self):
        par = {
            'method' : 'trust-constr',
            'max_dev' : 0.001                
        }
        return par
    
    
class Isenthalpic(crit_flow_model):
    """
    ISENTHAPLIC MODEL (non-homogeneous equilibrium) 
    E. ELIASI and G. S. LELLOUCHE2
    Two-phase critical flow
    International Journal of Multiphase Flow
    Volume 20, Supplement 1, 1994, Pages 91-168, ISSN 0301-9322,
    """
    def __init__(self, fluid, p_0, h_0 = None, x_0 = None, 
                 slip_model = None, 
                 **solver_par):
        super().__init__(fluid, p_0, h_0 = h_0, x_0 = x_0)
        if slip_model is None:
            self.slip_model = 'HEM' #
        else:
            self.slip_model = slip_model
        
        self.solver_par = self.__default_solver_par()
        if solver_par:
            self.solver_par.update(solver_par)
        self.main()
        
    @property
    def props_gtE(self):
        return props_E(self.fluid, self.p_t, 1, in1 = 'p', in2 = 'Q',
                            inkl = ['d','h', 's'])
    @property
    def props_ltE(self):
        return props_E(self.fluid, self.p_t, 0, in1 = 'p', in2 = 'Q',
                                inkl = ['d', 'h', 's'])   
    @property
    def props_t(self):
        return props_E(self.fluid, self.p_t, self.x_t, in1 = 'p', in2 = 'Q',
                                inkl = ['d', 'h', 's']) 
    @property
    def slip_ratio(self):
        if self.slip_model == 'HEM':
            return 1
        if self.slip_model == 'Moddy':
            return (self.props_ltE['d']/self.props_gtE['d']) ** (1 / 3)
        
    def slip_density(self, x_t):
        S = self.slip_ratio
        rho_gtE = self.props_gtE['d']
        rho_ltE = self.props_ltE['d']
        rho_tm = 1 / ((x_t / rho_gtE + S * (1 - x_t) / rho_ltE) *\
                       (x_t + (1 - x_t) / S ** 2) ** 0.5)
        return rho_tm
    
    def G2(self, p_t):
        self.p_t = p_t
        x_tmin = 0.00001
        x_tmax = 0.99999       
        x_t_guess = x_tmin
        res = minimize(self.__funcG_zero, 
                       x_t_guess,
                       method = 'TNC',                      
                       bounds = ((x_tmin, x_tmax),)
                       #options={'accuracy': 1e-4,}

        )
        
        self.x_t = res.x[0]
        G2 = self.G2Nrg(self.x_t)
        self.rel_dev = (self.__funcG_zero(self.x_t)/ G2) ** 0.5
        return G2
    
    def G2Nrg(self, x_t):
        h_0 = self.h_0
        h_gtE = self.props_gtE['h']
        h_ltE = self.props_ltE['h']
        h_t = x_t * h_gtE +  h_ltE * (1 - x_t)
        rho_tm = self.slip_density(x_t)
        G2 =  2 * (h_0 - h_t) * rho_tm ** 2
        return G2
        
    def G2Mom(self, x_t):
        p_0 = self.p_0
        p_t = self.p_t
        rho_tm = self.slip_density(x_t)
        G2 = 2 * (p_0 - p_t) * rho_tm  # assume G_0 = 0 (A_0 >> A_t), no friction
        return G2
    
    def solve_GC(self): 
        if self.x_0 < 0:
            p_tmin = p_min(self.fluid)
            p_tmax = self.p_sat_T0
        else:
            p_tmin = p_min(self.fluid)
            p_tmax = self.p_0 * 0.99

        k = 0.87
        p_t_guess = k * p_tmax  
        if p_tmax < p_tmin:
            #trend
            return
        res = minimize(self.__funcGc_min, 
                       p_t_guess,
                       method = self.solver_par['method'],
                       bounds = ((p_tmin, p_tmax),),
                       tol = 1e-4
        )
        #print(res)
        self.p_t = res.x[0]
        self.Gc = self.G2(self.p_t) ** 0.5
    
    def main(self):
        try:
            self.solve_GC()
        except ValueError:
            pass
    
    def __funcG_zero(self, x_t):
        #print(abs(self.G2Mom(x_t) - self.G2Nrg(x_t)))
        return abs(self.G2Mom(x_t) - self.G2Nrg(x_t))
    
    def __funcGc_min(self, p_t): 
        return -self.G2(p_t)
    #def solve_x_t()
    def __default_solver_par(self):
        par = {
            'method' : 'trust-constr',               
        }
        return par

def Test_CFM_general(crit_flow_model):
    try:
        print('stagnation pressure:', p_0 / 1e5, 'bar')
        print('stagnation enthalpy:',  crit_flow_model.h_0/ 1e3, 'kJ/kg')
        print('stagnation temperature:',  crit_flow_model.props_0['T'] - 273.15, '°C')
        print('throat pressure:', crit_flow_model.p_t / 1e5, 'bar')
        print('mass flux:', crit_flow_model.Gc, 'kg/m2s')
        print('dev:', crit_flow_model.rel_dev)
    except:
        pass
                 
def Test_HF(fluid, p_0, x_0):
    print('========= Henry-Fauske model ===========')
    crit_flow_model =  Henry_Fauske(fluid, p_0, x_0 = x_0, 
                                                 #method = 'TNC'
                                                )
    Test_CFM_general(crit_flow_model)
    return crit_flow_model

def Test_HF_mod(fluid, p_0, x_0):
    print('========= Henry-Fauske model (mod.) ===========')
    crit_flow_model =  Henry_Fauske_modified2(fluid, p_0, x_0 = x_0, 
                                                 #method = 'TNC'
                                                )
    Test_CFM_general(crit_flow_model)
    return crit_flow_model

def Test_HEM(fluid, p_0, x_0):
    print('========= HEM =========')
    crit_flow_model =  HEM(fluid, p_0, x_0 = x_0, 
                                             #method = 'TNC'
                                            )
    Test_CFM_general(crit_flow_model)    
    return crit_flow_model

def Test_ihM(fluid, p_0, x_0):
    print('========= isenthalpic model ===========')
    crit_flow_model = Isenthalpic(fluid, p_0, x_0 = x_0, slip_model = 'HEM'
                                                )    
    Test_CFM_general(crit_flow_model)    
    return crit_flow_model


if __name__ == '__main__':
    start = timeit.default_timer()  
    fluid = 'CO2' 
    p_0 = 8e5
    x_0 = -0.03
    """
    #ih =  Test_ihM(fluid, p_0, x_0)
    #ihl = [ih]
    
    HFm = Test_HF_mod(fluid, p_0, x_0, )
    HFml = [HFm]
    
    HF = Test_HF(fluid, p_0, x_0, )
    HFl = [HF]
    
    hem = Test_HEM(fluid, p_0, x_0)
    heml = [hem]
     
    """
    #p = [2e5, 4e5, 8e5, 16e5, 20e5, 60e5]
    start = timeit.default_timer()
    """
    x = [-0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 
         -0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001,
         0, 
         0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,
         0.01, 0.02, 0.03, 0.04,0.05,0.06,0.07,0.08,0.09,
         0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
         0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99]
    #x = [-0.009, -0.008, -0.007, -0.006, -0.005, -0.004, -0.003, -0.002, -0.001]
    #x = [-0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02]
    Gc = np.zeros(len(x))
    pt = np.zeros(len(x))
    HFl = []
    for i,x_0 in enumerate(x):
        HF = Test_HEM(fluid, p_0, x_0)
        HFl.append(HF)
        Gc[i] = HF.Gc
        pt[i] = HF.p_t / 1e5
    d = 0.1e-3
    #mdot = crit_mass_flow_model.Gc * (1/4*np.pi*d ** 2)
    
    5.2, 5.25, 5.3, 5.35, 5.4, 5.45, 5.5, 5.55, 5.6, 5.65, 
                  5.7, 5.8, 5.9, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5,
                  10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                  20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70
    
    """
    
    p = np.array([5.3,
          
         20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70])
    p = p*1e5
    Gc = np.zeros(len(p))
    pt = np.zeros(len(p))
    HFl = []
    for i,p_0 in enumerate(p):
        HF = Test_HF_mod(fluid, p_0, x_0)
        HFl.append(HF)
        Gc[i] = HF.Gc
        pt[i] = HF.p_t / 1e5
    d = 0.1e-3
    print('Time: ', timeit.default_timer() - start)
    """
    start = timeit.default_timer()
    fluid = 'CO2'
    output_code = 't'
    input_code = 'ps' 
    p = 10e5
    h = 200000  
    
    props_t1 = props_E(fluid, p * 1, h, 
                                    in1 = 'p', in2 = 'h',
                                    inkl = ['T', 'd', 'h', 's','Ds/Dp|E']
                                    )
    props_t2 = props_E(fluid, p * 2, h, 
                                    in1 = 'p', in2 = 'h',
                                    inkl = ['T', 'd', 'h', 's','Ds/Dp|E']
                                    )
    props_t3 = props_E(fluid, p* 3, h, 
                                    in1 = 'p', in2 = 'h',
                                    inkl = ['T', 'd', 'h', 's','Ds/Dp|E']
                                    )
    
    #sxl = CP.PropsSI('d(Smass)/d(P)|Sigma', 'P', p, 'Q', 0, 'CO2')
    
    
    stop = timeit.default_timer()
    print('Time: ', stop - start)
    """