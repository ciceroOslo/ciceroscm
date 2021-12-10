import numpy as np
import pandas as pd

import sys

def read_components(filename):
    """
    Read in components to be considered
    """
    df_gas = pd.read_csv(filename, delim_whitespace=True, index_col=0)
    df_gas.rename(
        columns={"TAU1(YEARS)": "TAU1", "NATURAL_EMISSIONS": "NAT_EM"}, inplace=True
    )
    return df_gas

def _rs_function(time):
    """
    Calculate pulse response function for mixed layer
    time is the year index*idtm + i, i.e. the month number
    """
    if time < 2. :
        rs = 0.12935 + 0.21898*np.exp(-time/0.034569) + 0.17003*np.exp(-time/0.26936) + 0.24071*np.exp(-time/0.96083) + 0.24093*np.exp(-time/4.9792)
    else:
        rs = 0.022936 + 0.24278*np.exp(-time/1.2679) + 0.13963*np.exp(-time/5.2528) + 0.089318*np.exp(-time/18.601) + 0.03782*np.exp(-time/68.736) + 0.035549*np.exp(-time/232.3)
    return rs

def _rb_function(time):
    """
    Calculate biotic decay function
    time is the year index*idtm + i, i.e. the month number
    """
    rb = 0.70211*np.exp(-0.35*time) +13.4141E-3*np.exp(-time/20.) - 0.71846*np.exp(-55*time/120.) + 2.9323E-3*np.exp(-time/100.) 
    return rb

def read_natural_emissions(filename, component, startyear=1750, endyear=2500):
    """
    Reading in single column natural emissions file
    """
    df = pd.read_csv(filename, header = None, names=[component], index_col = False)
    df['year'] = np.arange(startyear, endyear+1)
    df = df.set_index('year')
    return df

class ConcentrationsEmissionsHandler:
    """
    Class to handle concentrations
    and emissions input for ciceroscm
    """
    def __init__(self, gaspam_file, conc_file, emis_file, pamset, nat_ch4_file, nat_n2o_file,  conc_run = True, ref_year = 2010, lifetime_mode = 'TAR'):
        self.conc_run = conc_run
        self.df_gas = read_components(gaspam_file)
        print(self.df_gas.head())
        #sys.exit(4)
        self.conc = {}
        self.forc = {}
        if conc_file is not None:
            self.conc_in = self.read_inputfile(conc_file)
        if emis_file is not None:
            self.emis = self.read_inputfile(emis_file)
            for tracer in self.df_gas.index:
                if tracer != "CO2.1":
                    self.conc[tracer] = {}
                    self.forc[tracer] = []
            self.forc["Total_forcing"] =[]                
            self.emis.rename(columns={"CO2":"CO2_FF", "CO2.1":"CO2_AFOLU"}, inplace=True)
        self.nat_emis_ch4 = read_natural_emissions(nat_ch4_file, "CH4")
        self.nat_emis_n2o = read_natural_emissions(nat_n2o_file, "N2O")
        self.ref_yr = ref_year
        self.pamset = pamset
        self.years = []
        self.idtm = 12
        self.yCO2 = 0.0
        self.sCO2 = []
        self.xCO2 = 278.
        self.emCO2_prev = 0.
        self.dfnpp = []
        self.dfnpp.append(0)
        self.ss1 = 0.5*(self.emis["CO2_FF"].values[0]+self.emis["CO2_AFOLU"].values[0]+self.df_gas['NAT_EM']["CO2"])/2.123
        self.lifetime_mode = lifetime_mode
        self.sums = 0.

    def read_inputfile(self, input_file):
        """
        Read input from emissions or concentrations file
        """
        df_input = pd.read_csv(input_file, delim_whitespace=True, index_col=0, skiprows=[1,2,3])
        return df_input
        

    def calculate_strat_quantities(self, yr):
        """
        Calculating sumcl and sumbr in stratosphere
        """
        if yr < self.years[0]:
            return 0.0, 0.0
        chlor_dict = {"CFC-11":3, "CFC-12":2, "CFC-113":3, "CFC-114":2, "CFC-115":1, "CCl4":4,"CH3CCl3":3, "HCFC-22":22, "HCFC-141b":2, "HCFC-123":2,"HCFC-142b":1}
        #More Halons?
        brom_dict = {"H-1211":1, "H-1301":1}
        sumcl = 0
        sumbr = 0
        for comp,mult in chlor_dict.items():
            sumcl = sumcl + (mult*self.conc[comp][yr])**1.7
        for comp,mult in chlor_dict.items():
            sumbr = sumbr + (mult*self.conc[comp][yr])       
        return sumcl, sumbr
    
    def conc2forc(self, yr, rf_luc, rf_sun):
        """
        Calculate forcing from concentrations
        """
        #Setting some constants. What are they? Should they be parametrisable?
        #(Etminan I guess...)
        a1 = -2.4e-7
        b1= 7.2e-4
        c1 = -2.1e-4
        a2 = -8.0e-6
        b2 = 4.2e-6
        c2 = -4.9e-6
        a3 = -1.3e-6
        b3 = -8.2e-6
        tot_forc = 0
        FN = 0
        FS = 0
        yr_ix = yr -self.years[0]
        c0_n2o = self.conc["N2O"][self.years[0]]
        c_n2o = self.conc["N2O"][yr]
        
        
        # Finish per tracer calculations, add per tracer to printable df and sum the total
        for tracer  in self.df_gas.index:
            q = 0
            try:
                value = self.conc[tracer][yr]
                c0 = self.conc[tracer][self.years[0]]
            except:
                #Tracer with no concentration values
                value = 0
                c0 = 0
            if tracer == "CO2":
                q = (a1*(value-c0)**2+b1*np.abs(value-c0)+c1*0.5*(c0_n2o+c_n2o)+5.36)*np.log(value/c0) #Etminan et al 2016 RF
            elif tracer == "CH4":
                fmn = 0.47*np.log(1+2.01E-5*(value*c0_n2o)**(0.75) + 5.31E-15*value*(value*c0_n2o)**(1.52))
                fmn0 = 0.47*np.log(1+2.01E-5*(c0*c0_n2o)**(0.75) + 5.31E-15*c0*(c0*c0_n2o)**(1.52))
        
                #q = 0.036*(SQRT(c)-SQRT(c0))-(fmn-fmn0) ! TAR forcing
                q =  (a3*0.5*(value+c0)+b3*0.5*(c_n2o+c0_n2o)+0.043)*(np.sqrt(value)-np.sqrt(c0)) #Etminan et al 2016 RF
                #Feedback factor: Smith et al 2018
                q = 1.0/1.14*q #  + FORC_PERT(yr_ix,trc_ix))
            elif tracer == "N2O":
                c0_co2 = self.conc["CO2"][self.years[0]]
                c0_ch4 = self.conc["CH4"][self.years[0]]
                fmn = 0.47*np.log(1+2.01E-5*(value*c0_co2)**(0.75) + 5.31E-15*c0_co2*(value*c0_co2)**(1.52))
                fmn0 = 0.47*np.log(1+2.01E-5*(c0*c0_co2)**(0.75) + 5.31E-15*c0_co2*(c0*c0_co2)**(1.52))
                #q = 0.12*(SQRT(c)-SQRT(c0))-(fmn-fmn0) ! TAR FORCING
                q = (a2*0.5*(self.conc["CO2"][yr]+c0_ch4)+b2*0.5*(value+c0)+c2*0.5*(self.conc["CH4"][yr]+c0_co2) + 0.117)*(np.sqrt(value)-np.sqrt(c0)) #Etminan et al 2016 RF
            elif tracer == "SO2":
                #Natural emissions
                #(after IPCC TPII on simple climate models, 1997)
                #enat = 42.0 Not used, why is this here?
                #Emission in reference year
                erefyr = self.emis[tracer][self.ref_yr]
                if erefyr != 0:
                    frac_em = self.emis[tracer][yr]/erefyr
                    q = self.pamset["qdirso2"]*frac_em
            elif tracer == "SO4_IND":
                #Natural emissions
                #(after IPCC TPII on simple climate models, 1997)
                #enat = 42.0 Not used, why is this here?
                #Emission in reference year
                erefyr = self.emis["SO2"][self.ref_yr]
                if erefyr != 0:
                    frac_em = self.emis["SO2"][yr]/erefyr
                    q = self.pamset["qindso2"]*frac_em
            elif tracer == "LANDUSE":
                q = rf_luc 
            elif tracer == "BC":
                erefyr = self.emis[tracer][self.ref_yr]
                if erefyr == 0:
                    q = 0
                else:
                    frac_em = self.emis[tracer][yr]/erefyr
                    q = self.pamset["qbc"]*frac_em
            elif tracer == "OC":
                erefyr = self.emis[tracer][self.ref_yr]
                if erefyr == 0:
                    q = 0
                else:
                    frac_em = self.emis[tracer][yr]/erefyr
                    q = self.pamset["qoc"]*frac_em
                    
            elif tracer in self.df_gas.index and self.df_gas['ALPHA'][tracer] != 0:
                q = (value - c0)*self.df_gas['ALPHA'][tracer] #+forc_pert
            elif tracer == "TROP_O3":
                if self.conc_run or yr_ix == 0:
                    # Uses change in CO2_FF emissions
                    q = (self.emis["CO2_FF"][yr]-self.emis["CO2_FF"][self.years[0]])/(self.emis["CO2_FF"][self.years[0]]-self.emis["CO2_FF"][self.years[0]])*self.pamset["qo3"]
                else:
                   
                    #RBS101115     
                    #IF (yr_ix.LT.yr_2010) THEN ! Proportional to TROP_O3 build-up
                    # Rewritten a bit, place to check for differences...
                    q1 = self.forc[tracer][0]
                    q = q1 +(value-c0)/(30.0-c0)* (self.pamset["qo3"]-q1)
            elif tracer == "STRAT_O3":
                sumcl, sumbr = self.calculate_strat_quantities(yr-3)
                q = -(0.287737*(0.000552*(sumcl)+3.048*sumbr))/1000.
            elif tracer == "STRAT_H2O":
                q = 0.15*1.14*self.forc["CH4"][yr-self.years[0]]# + FORC_PERT(yr_ix,trc_ix)
            elif tracer == "OTHER":
                #Possible with forcing perturbations for other
                #cmponents such as contrails, cirrus etc...
                pass
            self.forc[tracer].append(q)  #+ FORC_PERT(yr_ix,trc_ix)
            #Calculating hemispheric forcing:
            if tracer == "SO2" or tracer == "SO4_IND":
                FN = FN + q*1.6
                FN = FS + q*0.4
            elif tracer == "TROP_O3":
                #1.29+0.74 != 2, update/check
                FN = FN + q*1.29 #1.3
                FN = FS + q*0.74 #0.7
            else:
                FN = FN + q
                FN = FS + q
            tot_forc = tot_forc + q
        # Adding solar forcing
        tot_forc = tot_forc + rf_sun
        self.forc["Total_forcing"].append(tot_forc)
        FN = FN + rf_sun
        FS = FS + rf_sun
               
        return tot_forc, FN, FS

    def emi2conc(self, yr):
        """
        Calculate concentrations from emissions
        """
        #Do per tracer emissions to concentrations, update concentrations df
        # NBNB! Remember to move calculation of  Trop_O3 concentration
        # here from conc2forc.
        """
         # First calculate O3 concentrations (AS IN TAR p269, table 4.11 B)
        CONC(yr_ix,trc_ix) =  30.0 + 6.7 * &
          (ALOG(CONC(yr_ix,trcID("CH4")))-ALOG(1832.0)) &   !ALOG(1700.0))  !Concentration in 2010 &
          + 0.17 * (EMISSIONS(yr_ix,trcID("NOx"))-EM2010(trcID("NOx"))) &
          + 0.0014 * (EMISSIONS(yr_ix,trcID("CO"))-EM2010(trcID("CO"))) &
          + 0.0042 *(EMISSIONS(yr_ix,trcID("NMVOC"))-EM2010(trcID("NMVOC")))
        """
        self.years.append(yr)
        if self.conc_run:
            self.fill_one_row_conc(yr)
            return
        self.add_row_of_zeros_conc(yr)

        for tracer in self.df_gas.index:
            #something something postscenario...?
            if self.df_gas['CONC_UNIT'][tracer] == "-":
                #Forcing calculated from emissions
                continue
            if tracer == "CO2":
                self.co2em2conc(yr)
                continue
            if yr > self.years[0]:
                xc = self.conc[tracer][yr -1]
            else:
                xc = self.conc_in[tracer][yr]

            q = 1.0/self.df_gas["TAU1"][tracer]
            
            if tracer == "CH4":
                self.df_gas['NAT_EM'][tracer] = self.nat_emis_ch4["CH4"][yr]
                q = self.methane_lifetime(q,xc,yr)
            if tracer == "N2O":
                self.df_gas['NAT_EM'][tracer] = self.nat_emis_n2o["N2O"][yr]

            q = q/self.idtm
            for i in range(self.idtm):
                em = self.emis[tracer][yr]
                em = em + self.df_gas['NAT_EM'][tracer] #natural emissions, from gasspamfile
                ach = xc
                em = em/self.idtm
                pc = em/self.df_gas["BETA"][tracer]
                ach = pc/q + (ach-pc/q)*np.exp(-q)
                xc = ach
            self.conc[tracer][yr] = xc

    def methane_lifetime(self,q,xc, yr):
        """
        Calculate methane concentrations from emissions
        """
        ch4_wigley_exp = -0.238
        if self.lifetime_mode == 'TAR':
            #1751 is reference conc in 2000
            dlnOH =  -0.32*(np.log(xc) - np.log(1751.0)) + 0.0042 *(self.emis["NOx"][yr]- self.emis["NOx"][2000]) - 0.000105 * (self.emis["CO"][yr]- self.emis["CO"][2000]) - 0.000315 * (self.emis["NMVOC"][yr]- self.emis["NMVOC"][2000])
            q = q * (dlnOH + 1)
        elif self.lifetime_mode == 'CONSTANT':
            q = 1.0/12.0
        else:
            q = q*(((xc/1700.))**(ch4_wigley_exp))

        q = q + 1./self.df_gas["TAU2"]["CH4"] + 1./self.df_gas["TAU3"]["CH4"]
        
        return q
        

    def co2em2conc(self, yr):
        """
        Calculate co2 concentrations from emissions
        """  
        #Fertilisation factor
        beta_f = 0.287

        #Area of the ocean (m^2)      
        ocean_area = 3.62E+14           

        #Gas exchange coefficient (air <--> ocean) (yr^-1*m^-2)
        coeff = 1.0/(ocean_area*9.06)

        #TIMESTEP (YR)
        dt = 1.0/self.idtm

        #Conversion factor ppm/kg --> umol*/m3
        c = 1.722E+17

        #USING MIXED LAYER DEPTH = 75 metres
        h = 75.0

        cc1 = dt*ocean_area*coeff/(1 + dt*ocean_area*coeff/2.)
        yr_ix = yr - self.years[0]
        #Monthloop:
        for i in range(self.idtm):
            it = yr_ix*self.idtm + i
            sumf = 0.0
            #Net emissions, including biogenic fertilization effects
            if it > 0:
                self.dfnpp.append(60*beta_f*np.log(self.xCO2/278.))
            if it > 1:
                for j in range(1,it):
                    sumf = sumf + self.dfnpp[it]*_rb_function(it)
            ffer = self.dfnpp[it] - dt*sumf
            emCO2 = self.emis["CO2_FF"][yr] + self.emis["CO2_AFOLU"][yr] - ffer
            emCO2 = emCO2 + self.df_gas['NAT_EM']["CO2"]
            emCO2 = emCO2/2.123

            if it == 0:
                ss2 = self.ss1
                self.sums = 0.
            else:
                ss2 = 0.5*emCO2/(ocean_area*coeff) - self.yCO2/(dt*ocean_area*coeff)
                self.sums = self.sums + self.emCO2_prev/(ocean_area*coeff) - self.sCO2[it-1]
            self.sCO2.append(cc1*(self.sums + self.ss1 + ss2))
            self.emCO2_prev = emCO2
            sumz = 0.
            for j in range(it):
                sumz = sumz + self.sCO2[j]*_rs_function(it-j)

            zCO2 = c*coeff*dt/h*(sumz+0.5*self.sCO2[it])
            self.yCO2 = 1.3021*zCO2 + 3.7929E-3*(zCO2**2) + 9.1193E-6*(zCO2**3) + 1.488E-8*(zCO2**4) + 1.2425E-10*(zCO2**5)
            self.xCO2 = self.sCO2[it] + self.yCO2 + 278.
            
        self.conc["CO2"][yr] = self.xCO2
            

    def fill_one_row_conc(self, yr):
        """
        Fill in one row of concentrations in conc_dict
        """
        for tracer in self.conc:
            if tracer in self.conc_in:
                self.conc[tracer][yr] = self.conc_in[tracer][yr]
            else:
                self.conc[tracer][yr] = 0
            
    def add_row_of_zeros_conc(self,yr):
        """
        Fill in one row of concentrations in conc_dict
        """
        for tracer in self.conc:
            self.conc[tracer][yr] = 0
