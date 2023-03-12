# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:01:40 2019

@author: Yixia
"""

import pandas as pd
import trend
import os

default_path = "C:\\Stoffdaten\\Trend 3.0\\" 
error_code_csv_name = "error_code.csv"
error_code_csv = os.path.join(default_path, error_code_csv_name)

#Prints all available TREND functions to the console.
# print(trend.__doc__)            


class trendprop():
    valid_input_code = [
        'TP', 'TD', 'PH', 'PS', 'TQ', 'PQ',
        'PT', 'DT', 'HP', 'SP', 'QT', 'QP',                            
        'PLIQ', 'PVAP', 'TLIQ', 'TVAP', 
        'TSUBV+', 'TSUBS+', 'PSUBV+', 'PSUBS+', 
        'TMLTL+', 'TMLTS+', 'PMLTS+', 'PMLTL+'
        ]
    
    t = 0
    d = 0
    p = 0
    u = 0 
    h = 0
    s = 0
    g = 0
    a = 0
    cp = 0
    cv = 0
    ws = 0
    b = 0
    c = 0
    cp0 = 0
    q = 0
    hom = 0
    k_s = 0
    def __init__(self, 
                 fluidname, moles = '1', eos_indicator = '1', 
                 input_code = None, prop1 = 0, prop2 = 0,
                 unit_T = 'K', unit_p = 'MPa', unit_m = 'mol', 
                 path = default_path,
                 err_capt = None,
                 ):        
        if err_capt is None:
            err_capt = False
        self.err_capt = err_capt
        self.path = path 
        self.fluidname = fluidname
        self.moles = moles       
        self.eos_indicator = eos_indicator
        if input_code is not None:
            ic1, self.prop1 = self.__prepare_unit(input_code[0], prop1, 
                unit_T = unit_T, unit_p = unit_p, unit_m = unit_m)
            if len(input_code) <= 3:
                ic2, self.prop2 = self.__prepare_unit(input_code[1], prop2, 
                unit_T = unit_T, unit_p = unit_p, unit_m = unit_m)
                self.input_code = ic1 + ic2 + input_code[2:]
            else:
                self.prop2 = prop2
                self.input_code = ic1 + input_code[1:]
                  

        
    def allprop(self):
        self.t,self.d,self.p,self.u,self.h,self.s,self.g,self.a,self.cp,self.cv,self.ws,self.b,self.c,self.cp0,self.q = trend.allprop(
                    self.input_code, 
                    self.prop1,
                    self.prop2,
                    self.fluidname,
                    self.moles,
                    self.eos_indicator,
                    self.path
                )
        self.__check_error(self.t)
        self.__check_error(self.d)
        self.__check_error(self.p)
        self.__check_error(self.u)
        self.__check_error(self.h)
        self.__check_error(self.s)
        self.__check_error(self.g)
        self.__check_error(self.a)
        self.__check_error(self.cp)
        self.__check_error(self.cv)
        self.__check_error(self.ws)
        self.__check_error(self.b)
        self.__check_error(self.c)
        self.__check_error(self.cp0)
        self.__check_error(self.q)
        
        
    def allprop_hom(self):
        """
        only supports inputcode "TD" and "TP"
        """

        phase, self.t,self.d,self.p,self.u,self.h,self.s,self.g,self.a,self.cp,self.cv,self.ws,self.b,self.c,self.cp0 = trend.allprop_hom(
                    self.input_code, 
                    self.prop1,
                    self.prop2,
                    self.fluidname,
                    self.moles,
                    self.eos_indicator,
                    self.path
                ) 
        self.__check_error(self.t)
        self.__check_error(self.d)
        self.__check_error(self.p)
        self.__check_error(self.u)
        self.__check_error(self.h)
        self.__check_error(self.s)
        self.__check_error(self.g)
        self.__check_error(self.a)
        self.__check_error(self.cp)
        self.__check_error(self.cv)
        self.__check_error(self.ws)
        self.__check_error(self.b)
        self.__check_error(self.c)
        self.__check_error(self.cp0)
        self.hom = 1
        
    def molar_mass(self):
        self.M = trend.mw( 
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        )
        self.__check_error(self.M)    
        return self.M
    
    
    def temperature(self):
        self.t = trend.t_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        )
        
        self.__check_error(self.t)        
        return self.t
    
    
    def density(self):
        self.d = trend.d_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        )
        
        self.__check_error(self.d)        
        return self.d
    
    
    def pressure(self):
        self.p = trend.p_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
         )
        
        self.__check_error(self.p)
        return self.p
    
    
    def enthalpy(self):
        self.h = trend.h_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        )
        self.__check_error(self.h)
        return self.h
    
    
    def entropy(self):
        self.s = trend.s_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        )
        self.__check_error(self.s)
        return self.s
    
    
    def quality(self):
        self.q = trend.q_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        )
        self.__check_error(self.q)
        return self.q
    
    def isobaric_heat_capacity(self):
        self.cp = trend.cp_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        ) 
        self.__check_error(self.cp)
        return self.cp
    
    def isochoric_heat_capacity(self):
        self.cv = trend.cv_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        ) 
        self.__check_error(self.cv)
        return self.cv
    
    def isentropic_expansion_coefficient(self):
        self.k_s = trend.expans_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        ) 
        self.__check_error(self.k_s)
        return self.k_s
    
    def isentropic_compressibility(self):
        self.beta_s = trend.comps_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        ) 
        self.__check_error(self.beta_s)
        return self.beta_s
    
    def viscosity(self):
        self.vs = trend.visdyn_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        ) 
        self.__check_error(self.vs)
        return self.vs
    
    def thermal_conductivity(self):
        self.l = trend.tcx_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        ) 
        self.__check_error(self.l)
        return self.l
    
    def Dp_DT_v(self):
        self.Dp_Dt_v = trend.dpdt_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        )
        self.__check_error(self.Dp_Dt_v)
        return self.Dp_Dt_v
    
    def Dp_Dd_T(self):
        self.Dp_Dd_T = trend.dpdd_eos(
            self.input_code, 
            self.prop1,
            self.prop2,
            self.fluidname,
            self.moles,
            self.eos_indicator,
            self.path
        )
        self.__check_error(self.Dp_Dd_T)
        return self.Dp_Dd_T
    
    def triple(self):
        triple = trend.triple_point(self.fluidname, self.path)
        return triple
    
    
    def critical(self):
        t_cri = trend.t_crit(self.fluidname, self.eos_indicator, self.path)
        p_cri = trend.p_crit(self.fluidname, self.eos_indicator, self.path)
        return t_cri, p_cri
    
    def Prop(self, output_code, unit_T = 'K', unit_p = 'MPa', unit_m = 'mol'):
        if output_code.upper() == 'T':
            if unit_T == 'C':
                return self.temperature() - 273.15
            elif unit_T == 'K':
                return self.temperature()
            else:
                print('defalut unit K is used')
                return self.temperature()
            
        if output_code.upper() == 'P':
            if unit_p == 'bar':
                return self.pressure() * 10
            elif unit_p == 'Pa':
                return self.pressure() * 1e6
            elif unit_p == 'MPa':
                return self.pressure()
            else:
                print('defalut unit MPa is used')
                return self.pressure()

        if output_code.upper() == 'D':
            if unit_m == 'kg':
                return self.density() * self.molar_mass()
            elif unit_m == 'mol':
                return self.density()
            else:
                print('defalut unit mol/m3 is used')
                return self.density()
            
        if output_code.upper() == 'V':
            if unit_m == 'kg':
                return 1 / self.density() / self.molar_mass()
            elif unit_m == 'mol':
                return 1 / self.density()
            else:
                print('defalut unit m3/mol is used')
                return 1 / self.density()
            
        if output_code.upper() == 'H':
            if unit_m == 'kg':
                return self.enthalpy() / self.molar_mass()
            elif unit_m == 'mol':
                return self.enthalpy()
            else:
                print('defalut unit J/mol is used')
                return self.enthalpy()
            
        if output_code.upper() == 'S':
            if unit_m == 'kg':
                return self.entropy() / self.molar_mass()
            elif unit_m == 'mol':
                return self.entropy()
            else:
                print('defalut unit J/(mol K) is used')
                return self.entropy()
            
        if output_code.upper() == 'CP':
            if unit_m == 'kg':
                return self.isobaric_heat_capacity() / self.molar_mass()
            elif unit_m == 'mol':
                return self.isobaric_heat_capacity()
            else:
                print('defalut unit J/(mol K) is used')
                return self.isobaric_heat_capacity()
        
        if output_code.upper() == 'CV':
            if unit_m == 'kg':
                return self.isochoric_heat_capacity() / self.molar_mass()
            elif unit_m == 'mol':
                return self.isochoric_heat_capacity()
            else:
                print('defalut unit J/(mol K) is used')
                return self.isochoric_heat_capacity() 
                
        if output_code.upper() == 'BETAS':
            if unit_p == 'bar':
                return self.isentropic_compressibility() / 10
            elif unit_p == 'Pa':
                return self.isentropic_compressibility() / 1e6
            elif unit_p == 'MPa':
                return self.isentropic_compressibility()
            else:
                print('defalut unit MPa is used')
                return self.pressure() 
            
        if output_code.upper() == 'VS':
            if unit_p == 'bar':
                return self.viscosity() / 1e11
            elif unit_p == 'Pa':
                return self.viscosity() / 1e6
            elif unit_p == 'MPa':
                return self.viscosity() / 1e12
            else:
                print('defalut unit uPa s is used')
                return self.viscosity()
        
        if output_code.upper() == 'L':
            return self.thermal_conductivity()
            
        if output_code.upper() == 'Q':
            return self.quality()        
        
        if output_code.upper() == 'KS':
            return self.isentropic_expansion_coefficient()
        
        if output_code == 'Dp/DT|d' or output_code == 'Dp/DT|v':
            if unit_p == 'bar':
                return self.Dp_DT_v() * 10
            elif unit_p == 'Pa':
                return self.Dp_DT_v() * 1e6
            elif unit_p == 'MPa':
                return self.Dp_DT_v()
        
        if output_code == 'Dp/Dd|T':
            
            if unit_m == 'kg':
                res = self.Dp_Dd_T() / self.molar_mass()
            if unit_p == 'bar':
                res = res * 10
            elif unit_p == 'Pa':
                res = res * 1e6
            elif unit_p == 'MPa':
                pass
            return res
        
    def Prop_hom(self, output_code, unit_T = 'K', unit_p = 'MPa', unit_m = 'mol'):
        if self.hom == 0:
            self.allprop_hom()
        if output_code.upper() == 'T':
            if unit_T == 'C':
                return self.t - 273.15
            elif unit_T == 'K':
                return self.t
            else:
                print('defalut unit K is used')
                return self.t
            
        if output_code.upper() == 'P':
            if unit_p == 'bar':
                return self.p * 10
            elif unit_p == 'Pa':
                return self.p * 1e6
            elif unit_p == 'MPa':
                return self.p
            else:
                print('defalut unit MPa is used')
                return self.pressure()

        if output_code.upper() == 'D':
            if unit_m == 'kg':
                return self.d * self.molar_mass()
            elif unit_m == 'mol':
                return self.d
            else:
                print('defalut unit mol/m3 is used')
                return self.d
            
        if output_code.upper() == 'V':
            if unit_m == 'kg':
                return 1 / self.d / self.molar_mass()
            elif unit_m == 'mol':
                return 1 / self.d
            else:
                print('defalut unit m3/mol is used')
                return 1 / self.d
            
        if output_code.upper() == 'H':
            if unit_m == 'kg':
                return self.h / self.molar_mass()
            elif unit_m == 'mol':
                return self.h
            else:
                print('defalut unit J/mol is used')
                return self.h
            
        if output_code.upper() == 'S':
            if unit_m == 'kg':
                return self.s / self.molar_mass()
            elif unit_m == 'mol':
                return self.s
            else:
                print('defalut unit J/(mol K) is used')
                return self.s
            
        if output_code.upper() == 'CP':
            if unit_m == 'kg':
                return self.cp / self.molar_mass()
            elif unit_m == 'mol':
                return self.cp
            else:
                print('defalut unit J/(mol K) is used')
                return self.cp
        
        if output_code.upper() == 'CV':
            if unit_m == 'kg':
                return self.cv / self.molar_mass()
            elif unit_m == 'mol':
                return self.cv
            else:
                print('defalut unit J/(mol K) is used')
                return self.cv
    

    def __prepare_unit(self, 
            input_code, input_prop, 
            unit_T = None, unit_p = None, unit_m = None):
        
        if input_code.upper() == 'T':
            if unit_T == 'C':
                input_prop = input_prop + 273.15
            elif unit_T == 'K':
                input_prop = input_prop
            else:
                print('defalut unit K is used')

            
        if input_code.upper() == 'P':
            if unit_p == 'bar':
                input_prop = input_prop / 10
            elif unit_p == 'Pa':
                input_prop = input_prop / 1e6
            elif unit_p == 'MPa':
                pass
            else:
                print('defalut unit MPa is used')


        if input_code.upper() == 'D':
            if unit_m == 'kg':
                input_prop = input_prop / self.molar_mass()
            elif unit_m == 'mol':
                pass
            else:
                print('defalut unit mol/m3 is used')
            
        if input_code.upper() == 'V':
            if unit_m == 'kg':
                input_prop = 1 / input_prop * self.molar_mass()
            elif unit_m == 'mol':
                input_prop = 1 / input_prop
            else:
                print('defalut unit m3/mol is used')
                input_prop = 1 / input_prop
            input_code = 'D'
            
        if input_code.upper() == 'H':
            if unit_m == 'kg':
                input_prop = input_prop * self.molar_mass()
            elif unit_m == 'mol':
                pass
            else:
                print('defalut unit J/mol is used')
            
        if input_code.upper() == 'S':
            if unit_m == 'kg':
                input_prop = input_prop * self.molar_mass()
            elif unit_m == 'mol':
                pass
            else:
                print('defalut unit J/(mol K) is used')
        
        return input_code, input_prop
    
    def __check_error(self, output):
        if self.err_capt is True:
            if output % 1  == 0:
                output = int(output)
            if output in error_code.index:
                error = error_code.loc[output, 'Comment']
                raise ValueError(error)
        return
            
def error_code():
    
    error_code = pd.read_csv(
            error_code_csv,
            sep = ";",
            header = 0,
            index_col = 0
        )
    """
    error_code = pd.DataFrame(
        data = {
                '-1111' : 'MIX: Wrong input(s) to a routine. These are internal errors, for example negative temperatures, wrong iFlash values etc. are caught ',
                '-1234' : 'PURE: PhaseDet, p(rho_it,T_it) --> press_input ',
                '-1235' : 'PURE: Temperature iteration T(p,h) failed ',
                '-1444' : 'Function Fug_DryIce: FUGCOPURE_CALC: Exponent > 700. Argument to exponential function too big in dry ice fugacity calculation ',
                '-2211' : 'Subroutine PHASEDET_PURE: Flash calculation did not converge ',
                '-2212' : 'Subroutine MAXWELL: Flash calculation with MAXWELL did not converge ',
                '-2222' : 'FLASH: Flash calculation did not converge ',
                '-2223' : 'Function TSUB_EQ: SUBLIMATION TEMPERATURE iteration did not converge (ancillary equation) ',
                '-2224' : 'PHASEDET_SOL: MELTING TEMPERATURE iteration did not converge (ancillary equation) ',
                '-2225' : 'PHASEDET_SOL: SUBLIMATION PRESSURE calculation failed ',
                '-2226' : 'PHASEDET_SOL: MELTING PRESSURE calculation failed ',
                '-2230' : 'PHASEDET_SOL: MELTING LINE: Given pressure not valid ',
                '-2231' : 'PHASEDET_SOL: MELTING LINE: Given temperature not valid ',
                '-2232' : 'PHASEDET_SOL: SUBLIMATION LINE: Given pressure not valid ',
                '-2233' : 'PHASEDET_SOL: SUBLIMATION LINE: Given temperature not valid ',
                '-2908' : 'cubic, LKP or generalized equation: No cp0 parameters set ',
                '-3333' : 'MIX: Step size reduction error in flash algorithm (step size becomes too small) ',
                '-4321' : 'MIX: Flash algorithm converged to unreasonable temperatures ',
                '-4322' : 'Subroutine FLASH_PHASEBOUNDARY_CALC: MIX flash algorithm converged to unreasonable pressures ',
                '-4323' : 'MIX: Flash algorithm converged to two phases with the same compositions ',
                '-4401' : 'Subroutine FLASH_PHASEBOUNDARY: MIX p,x(liq)-flash failed: no bubble point found ',
                '-4402' : 'Subroutine FLASH_PHASEBOUNDARY: MIX p,x(vap)-flash failed: no dew point found ',
                '-4403' : 'Subroutine FLASH_PHASEBOUNDARY: MIX t,x(liq)-flash failed: no bubble point found ',
                '-4404' : 'Subroutine FLASH_PHASEBOUNDARY: MIX t,x(vap)-flash failed: no dew point found ',
                '-4405' : 'MIX: ph-flash internal error ',
                '-4406' : 'Subroutine PHASEDET_PS: MIX ps-flash internal error ',
                '-4407' : 'MIX: td-flash internal error ',
                '-4444' : 'MIX: LUDECOMP failed to invert matrix ',
                '-5211' : 'Surface tension not implemented for mixtures ',
                '-5222' : 'Surface tension not defined in the homogenous region ',
                '-5223' : 'Property for mixtures not implemented ',
                '-5234' : 'No model for this property and fluid available ',
                '-5235' : 'Dielectric constant models cannot be mixed with model DE2 ',
                '-5242' : 'Error in viscosity calculation ',
                '-5243' : 'Existing model for the viscosity not implemented yet ',
                '-5244' : 'Existing model for the thermal conductivity not implemented yet ',
                '-5248' : 'Property only available for mixtures ',
                '-5249' : 'Property only available for binary mixtures ',
                '-5250' : 'Quality not implemented as input parameter for mixtures ',
                '-5300' : 'ECS: invalid reference fluid ',
                '-5301' : 'ECS: internal error ',
                '-5302' : 'ECS: iteration failed ',
                '-5501' : 'Creation of four points at low pressures failed: generation of starting point failed. Stability analysis failed. ',
                '-5503' : 'Phase envelope: generation of start values failed ',
                '-5507' : 'Subroutine PHASENV: MIX: phase envelope calculation failed (step reduction becomes too large) ',
                '-5508' : 'Subroutine PHASENV: MIX: phase envelope calculation failed ',
                '-5520' : 'Subroutine PHASENV: MIX: too many points found -> exceeds the vector storage capacity ',
                '-5530' : 'Subroutine pxdiag: T is larger than tc(1) and tc(2) ',
                '-5555' : 'Subroutine SATPLOT: MIX: Phase envelope calculation failed ',
                '-5566' : 'Subroutine SYSOFFEQS_PT: MIX: Internal error in pt-flash ',
                '-5666' : 'Function DSPIN_EOS: Spinodals: No valid input code, only TLIQ and TVAP possible ',
                '-5667' : 'Function DSPIN_CALC: Spinodals: Change in curvature between minimum and phase boundary (unreasonable result) ',
                '-5668' : 'Function DSPIN_EOS: Spinodales not implemented for mixtures ',
                '-5700' : '"Function PROP_EOS: Valid input combination for PNUMER_EOS is ""TD"" only "',
                '-5777' : 'Costald EOS: Different EOS types cannot be combined here ',
                '-5778' : 'Function PROP_EOS: Costald EOS: No valid input code, only TLIQ or TP possible or wrong property (only density possible) ',
                '-5888' : 'Function VSATTAIT: Costald EOS: Temperature out of valid temperature range ',
                '-5997' : 'Function PROP_EOS: RKM EOS: No valid input code, only TLIQ possible or wrong property (only density possible) ',
                '-5998' : 'Function PROP_EOS: ERKM EOS: No valid input code, only TLIQ or TP possible or wrong property (only density possible) ',
                '-5999' : 'Function v_RKM: Revised Klosek-McKinley EOS: Temperature and or pressure out of valid range ',
                '-6000' : 'Subroutine REF_CALC: No reference state set ',
                '-6001' : 'Acentric factor for generalized negative => calculation not possible ',
                '-6600' : 'No strict equation boundaries available ',
                '-6666' : 'Property in two-phase region undefined ',
                '-6667' : 'Property in homogeneous region undefined ',
                '-7000' : 'UNCERTAINTIES: Estimation of uncertainties only possible for pure fluids ',
                '-7001' : 'UNCERTAINTIES: Uncty File not found ',
                '-7002' : 'UNCERTAINTIES: Estimation of uncertainties oNly possible for Helmholtz EOS (see uncty file) ',
                '-7003' : 'UNCERTAINTIES: estimation in 2phase region not yet implemented ',
                '-7777' : 'MIX: Calculation of vapor fugacities failed ',
                '-7778' : 'MIX: Calculation of liquid fugacities failed ',
                '-7779' : 'MIX: Calculation of solid fugacities failed ',
                '-7877' : 'DLL could not be opened in Excel ',
                '-7878' : 'SETUP: Fluid name or path wrong ',
                '-7879' : 'SETUP Error: Opening fluid file ',
                '-7880' : 'SETUP Error: End of fluid-file is reached during read ',
                '-7881' : 'SETUP Error: Opening .MIX file ',
                '-7882' : 'SETUP Error: End of .MIX-file is reached during read ',
                '-7883' : 'SETUP Error: Opening atcoeff file ',
                '-7884' : 'SETUP Error: CAS list was not found ',
                '-7885' : 'SETUP Error: Opening SRK file ',
                '-7886' : 'SETUP Error: Opening PR file ',
                '-7887' : 'SETUP Error: Opening LKP file ',
                '-7888' : 'SETUP Error: Opening BIN SRK or PR file ',
                '-7889' : 'SETUP Error: Opening BIN LKP file ',
                '-7890' : 'SETUP Error: CAS-ID not found ',
                '-7891' : 'SETUP Error: Opening COSTALD file ',
                '-7892' : 'SETUP Error: Only Helmholtz EOS, PR, SRK, LKP and the corresponding mixing rules are hardcoded at the moment ',
                '-7893' : 'SETUP Error: Opening RKM files ',
                '-7999' : 'Wrong equation format (e.g. ECS) ',
                '-8877' : 'PURE: Density iteration failed ',
                '-8878' : 'Cubic EOS or LKP: Fluid not found ',
                '-8880' : 'Costald fluid not found ',
                '-8881' : 'Costald parameter file not found ',
                '-8888' : 'MIX: Density iteration failed ',
                '-8889' : 'find_crit_tpd: maximum number of iterations reached ',
                '-9000' : 'SETUP Error: Parameter file for generalized EOS not found ',
                '-9876' : 'NAN in FNRDERIVSMIX - fugcoef_pure, SRK, dlog( <0) ',
                '-9901' : 'Calculation not possible for pure fluids ',
                '-9902' : 'Calculation not possible for mixtures ',
                '-9903' : 'Calculation not possible for pure substances or binary mixtures ',
                '-9904' : 'Calculation not possible for solid phase of chosen substance ',
                '-9911' : 'Temperature input <= 0 ',
                '-9912' : 'Temperature <= Tminfluid ',
                '-9913' : 'Temperature >= Tmaxfluid ',
                '-9914' : 'Density >= Rhomaxfluid ',
                '-9915' : 'p >= pmelt --> solid phase ',
                '-9916' : 'Input parameters out of range (viscosity of water) ',
                '-9917' : 'ALLPROP_HOM: phase_ind out of range (phase_ind is an input parameter that specifies the phase which is tried) ',
                '-9918' : 'Temperature <= Tmin of transport equations ',
                '-9919' : 'Temperature >= Tmax of transport equations ',
                '-9920' : 'Density >= rhomax of transport equations ',
                '-9922' : 'Density input <= 0 ',
                '-9932' : 'Pressure input >= Pmaxfluid ',
                '-9933' : 'Pressure input <= 0 ',
                '-9934' : 'Property not available for this equation type ',
                '-9942' : 'Entropy input out of range ',
                '-9941' : 'Enthalpy input out of range ',
                '-9943' : 'T rho flash, density cannot be found in the given range of pmin to pmax ',
                '-9951' : 'One mole fraction exceeds 0 < x < 1 ',
                '-9952' : 'Sum of the mole fractions is not 1 ',
                '-9953' : 'Number of fluids =/ number of mole fractions =/ number of eqtypes ',
                '-9955' : 'Wrong Input (combination) ',
                '-9954' : 'Error in composition entry ',
                '-9957' : 'EOS_Indicator does not exist or wrong mixtype ',
                '-9958' : 'Cubic EOS cannot be mixed with Helmholtz EOS ',
                '-9959' : 'Error in EOS_Indicator entry, wrong format ',
                '-9960' : 'Error in fluid entry ',
                '-9961' : 'Option only valid for Helmholtz, PR, SRK, and LKP ',
                '-9971' : 'Negative radical in calculation of speed of sound ',
                '-9981' : 'Call of Flash_pure with Temp <= ttp or Temp >= tc ',
                '-9982' : 'Call of Flash_pure with press <= ptp or press >= pc ',
                '-12900' : 'Property not provided / implemented for solid equation ',
                '-12901' : 'Input combination not implemented for solids ',
                '-12902' : 'No solid model available for chosen substances ',
                '-14444' : 'Wrong value in solidtype(1). Either internal error, or the solid former does not exist ',
                '-15566' : 'Internal error in PhaseDet_sol ',
                '-15567' : 'Internal error in PhaseDet_ph_ps_sol ',
                '-15570' : 'Internal error in Flash_PhaseBoundary_sol ',
                '-15571' : 'Internal error in ptflash_sol_NC_4P ',
                '-18867' : 'Pressure iteration for solid failed ',
                '-19900' : 'Prediction of solids only possible in combination with Helmholtz EOS ',
                '-19912' : 'Temperature <= tmin of solid equation ',
                '-19914' : 'Density >= rhomax for solid ',
                '-19915' : 'p >= pmelt --> solid phase, no solid equation available or limits of solid equation exceeded ',
                '-19932' : 'Pressure >= pmax for solid equation ',
                '-19941' : 'Enthalpy or entropy input out of range ',
                '-44444' : 'PURE: Density iteration failed (rho_calc) ',
                '-101010' : 'Internal error in Phasedet_sol for mixtures ',
                '-111111' : 'Solid function not activated ',
                '-121212' : 'Error while generating start values for VLE pure, press or Temp <= 0 ',
                '-121213' : 'Invalid iFlash entry in VLE pure ',
                '-898964' : 'Flash_Pure_PhaseBoundary_calc: Density iteration failed ',
                '-898965' : 'Flash_Pure_PhaseBoundary_calc did not converge in specified number of iterations ',
                '-898966' : 'Flash_Pure_PhaseBoundary_calc: safety caution for numerical effects (see oil, Pvap is 1.d-10 --> very small steps) ',
                '-898967' : 'Flash_Pure_PhaseBoundary_calc: wrong solution (wrong slope) '
                },
        columns=['Error code ' , 'Comment']
    )
    """
    return error_code
error_code = error_code()

if __name__ == '__main__':
    import timeit
    start = timeit.default_timer()  
    """
    x = trendprop('co2', input_code = 'QP', prop1 = 1 , prop2 = 20, err_capt = False, unit_p = 'bar')   
    
    M = x.molar_mass()
    L = x.thermal_conductivity()
    T = x.temperature()
    Q = x.quality()
    D = x.Prop('D', unit_m = 'kg')
    H = x.Prop('H', unit_m = 'kg')
    S = x.Prop('S', unit_m = 'kg')
    V = x.Prop('V', unit_m = 'kg')
    Cp = x.Prop('CP', unit_m = 'kg')
    Cv = x.Prop('CV', unit_m = 'kg')
    Triple = x.triple()
    xls = [x]
    
    """
    
    
    # a = trendprop('co2', input_code = 'PQ', prop1 = 1.23e5, prop2 = 1, err_capt = True, unit_p = 'Pa', unit_T = 'C').Prop('H', unit_T = 'C',unit_m = 'kg')
    b = trendprop('co2', input_code = 'PSUBS+', prop1 = 1.01325e5, prop2 = 1, err_capt = True, unit_p = 'Pa', unit_T = 'C').Prop('T', unit_T = 'K',unit_m = 'kg')
    
    stop = timeit.default_timer()
    
    print('Time: ', stop - start)



"""        
#Set input parameters
input = 'TP'
prop1 = 300
prop2 = 0.1
fluidname = 'water'
fluids = fluidname
moles = '
eos_indicator = '1'
#Path to Trend Directory
path = 'C:\\Apps\\Trend 3.0\\'

values = ''
fluids = ['hcl', 'cl2', 'h2s', 'methane', 'hydrogen', 'co', 'argon', 'oxygen', 'nitrogen', 'water']
for i in range(len(fluids)):
    fluidname = fluids[i]
    t_crit = trend.t_crit(fluidname,eos_indicator,path)
    d_crit = trend.d_crit(fluidname,eos_indicator,path)
    values = values + fluidname + ";" + str(round(t_crit,2)) + ";" + str(round(d_crit,3)) + "\n"
print(values)

#Triple point
triple = trend.triple_point(fluidname,path)
#Allprop sub
t,d,p,u,h,s,g,a,cp,cv,ws,b,c,cp0,q = trend.allprop(input,prop1,prop2,fluidname,moles,eos_indicator,path)

#Density function
deos = trend.d_eos(input,prop1,prop2,fluidname,moles,eos_indicator,path)
#Flash Calculation
phasetype,x_phase,prop_phase,prop_overall,lnfug_phase,chempot_phase,phasefrac = trend.flash3(input,prop1,prop2,fluidname,moles,eos_indicator,path)
"""


