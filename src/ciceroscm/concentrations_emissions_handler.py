import numpy as np
import pandas as pd


def read_components(filename):
    """
    Read in components to be considered
    """
    df_gas = pd.read_csv(filename, delim_whitespace=True, index_col=0)
    df_gas.rename(
        columns={"TAU1(YEARS)": "TAU1", "NATURAL_EMISSIONS": "NAT_EM"}, inplace=True
    )
    return df_gas


class ConcentrationsEmissionsHandler:
    """
    Class to handle concentrations
    and emissions input for ciceroscm
    """
    def __init__(self, gaspam_file, conc_file, emis_file, pamset, conc_run = True, ref_year = 2010):
        self.conc_run = conc_run
        self.df_gas = read_components(gaspam_file)
        self.conc = {}
        self.forc = {}
        if conc_file is not None:
            self.conc_in = self.read_inputfile(conc_file)
        if emis_file is not None:
            self.emis = self.read_inputfile(emis_file)
            for tracer in self.emis.columns:
                print(tracer)
                if tracer != "CO2.1":
                    self.conc[tracer] = []
                    self.forc[tracer] = []
            self.forc["Total_forcing"] =[]
                
            self.emis.rename(columns={"CO2":"CO2_FF", "CO2.1":"CO2_AFOLU"}, inplace=True)
        self.ref_yr = ref_year
        self.pamset = pamset
        self.years = []
        

    def read_inputfile(self, input_file):
        """
        Read input from emissions or concentrations file
        """
        df_input = pd.read_csv(input_file, delim_whitespace=True, index_col=0, skiprows=[1,2,3])
        print(df_input.describe())
        return df_input

    def calculate_strat_quantities(self, yr_ix):
        """
        Calculating sumcl and sumbr in stratosphere
        """
        if yr_ix < 0:
            return 0.0, 0.0
        chlor_dict = {"CFC-11":3, "CFC-12":2, "CFC-113":3, "CFC-114":2, "CFC-115":1, "CCl4":4,"CH3CCl3":3, "HCFC-22":22, "HCFC-141b":2, "HCFC-123":2,"HCFC-142b":1}
        #More Halons?
        brom_dict = {"H-1211":1, "H-1301":1}
        sumcl = 0
        sumbr = 0
        for comp,mult in chlor_dict.items():
            sumcl = sumcl + (mult*self.conc[comp][yr_ix])**1.7
        for comp,mult in chlor_dict.items():
            sumbr = sumbr + (mult*self.conc[comp][yr_ix])       
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
        c0_n2o = self.conc["N2O"][0]
        c_n2o = self.conc["N2O"][yr_ix]
        
        
        # Finish per tracer calculations, add per tracer to printable df and sum the total
        for tracer  in self.conc:
            q = 0
            try:
                value = self.conc[tracer][yr_ix]
                c0 = self.conc[tracer][0]
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
                c0_co2 = self.conc["CO2"][0]
                c0_ch4 = self.conc["CH4"][0]
                fmn = 0.47*np.log(1+2.01E-5*(value*c0_co2)**(0.75) + 5.31E-15*c0_co2*(value*c0_co2)**(1.52))
                fmn0 = 0.47*np.log(1+2.01E-5*(c0*c0_co2)**(0.75) + 5.31E-15*c0_co2*(c0*c0_co2)**(1.52))
                #q = 0.12*(SQRT(c)-SQRT(c0))-(fmn-fmn0) ! TAR FORCING
                q = (a2*0.5*(self.conc["CO2"][yr_ix]+c0_ch4)+b2*0.5*(value+c0)+c2*0.5*(self.conc["CH4"][yr_ix]+c0_co2) + 0.117)*(np.sqrt(value)-np.sqrt(c0)) #Etminan et al 2016 RF
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
                    q = (self.emis["CO2_FF"][yr]-self.emis["CO2_FF"][self.yr_0])/(self.emis["CO2_FF"][self.ref_yr]-self.emis["CO2_FF"][self.yr_0])*self.pamset["qo3"]
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
                q = 0.15*1.14*self.forc["CH4"][yr]# + FORC_PERT(yr_ix,trc_ix)
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
        yr_ix = yr - self.years[0]
        idtm = 12
            
        #Fertilisation factor
        beta_f = 0.287

        #Area of the ocean (m^2)      
        ocean_area = 3.62E+14           

        #Gas exchange coefficient (air <--> ocean) (yr^-1*m^-2)
        coeff = 1.0/(ocean_area*9.06)

        #TIMESTEP (YR)
        dt = 1.0/idtm

        #Conversion factor ppm/kg --> umol*/m3
        c = 1.722E+17

        #USING MIXED LAYER DEPTH = 75 metres
        h = 75.0

    def fill_one_row_conc(self, yr):
        """
        Fill in one row of concentrations in conc_dict
        """
        for tracer in self.conc:
            if tracer in self.conc_in:
                self.conc[tracer].append(self.conc_in[tracer][yr])
            else:
                self.conc[tracer].append(0)       
