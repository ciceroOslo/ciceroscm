import numpy as np


sek_day = 86400
day_year = 365.

def _coefic(s, t, p):
    """
    Calculate denisty coefficients
    """
    coefic = (
        19652.21
        + 148.4206 * t
        - 2.327105 * t ** 2
        + 1.360477e-2 * t ** 3
        - 5.155288e-5 * t ** 4
        + 3.239908 * p
        + 1.43713e-3 * t * p
        + 1.16092e-4 * t ** 2 * p
        - 5.77905e-7 * t ** 3 * p
        + 8.50935e-5 * p ** 2
        - 6.12293e-6 * t * p ** 2
        + 5.2787e-8 * t ** 2 * p ** 2
        + 54.6746 * s
        - 0.603459 * t * s
        + 1.09987e-2 * t ** 2 * s
        - 6.1670e-5 * t ** 3 * s
        + 7.944e-2 * s ** 1.5
        + 1.6483e-2 * t * s ** 1.5
        - 5.3009e-4 * t ** 2 * s ** 1.5
        + 2.2838e-3 * p * s
        - 1.0981e-5 * t * p * s
        - 1.6078e-6 * t ** 2 * p * s
        + 1.91075e-4 * p * s ** 1.5
        - 9.9348e-7 * p ** 2 * s
        + 2.0816e-8 * t * p ** 2 * s
        + 9.1697e-10 * t ** 2 * p ** 2 * s
    )
    return coefic


def _denso(s, t):
    """
    Calculate density at p0=0
    """
    denso = (
        999.842594
        + 6.793952e-2 * t
        - 9.095290e-3 * t ** 2
        + 1.001685e-4 * t ** 3
        - 1.120083e-6 * t ** 4
        + 6.536332e-9 * t ** 5
        + 8.24493e-1 * s
        - 4.0899e-3 * t * s
        + 7.6438e-5 * t ** 2 * s
        - 8.2467e-7 * t ** 3 * s
        + 5.3875e-9 * t ** 4 * s
        - 5.72466e-3 * s ** 1.5
        + 1.0227e-4 * t * s ** 1.5
        - 1.6546e-6 * t ** 2 * s ** 1.5
        + 4.8314e-4 * s ** 2
    )
    return denso


def _density(p0, t0):
    """
    Calculate water denisity from equation of state
    """
    s = 35.0
    return _denso(s, t0) / (1.0 - p0 / _coefic(s, t0, p0))


class UpwellingDiffusionModel:
    """
    Class to handle energy budget upwelling and downwelling
    """

    def __init__(self, params):
        """
        Intialising
        """
        self.rlamdo = params["rlamdo"]
        self.rakapa = 1.e-4*params["akapa"]
        self.cpi = params["cpi"]
        self.w = params["w"]
        self.beto = params["beto"]
        self.threstemp = params["threstemp"]
        self.mixed = params["mixed"]
        self.rlamda = 1.0/params["LAMBDA"]
        self.foan = 0.61  # ocean fraction in the Nothern Hemisphere make changable?
        self.foas = 0.81  # ocean fraction in the Southern Hemisphere make changable?
        self.ebeta = 0.0 #Make changable?
        self.fnso =0.735
        self.lm = 40
        
        # Setting up dz height difference between ocean layers
        self.dz = np.ones(lm)
        self.dz[0] = self.mixed
        
        self.ldtime = 12.
        self.dt = 1/ldtime*sek_day*year_day
        self.setup_ebud()
        self.FNOLD = 0.0
        self.FSOLD = 0.0
        self.setup_sea_level_rise()

    def _band(self, a, b, c, d):
        """
        Calculate band
        """
        alfa = np.zeros(self.lm - 1)
        ans = np.zeros(self.lm)
        bbeta = np.zeros(self.lm)

        alfa[0] = -b[0] / a[0]
        bbeta[0] = d[0] / a[0]
        print(alfa)
        print(bbeta)
        for i in range(1, lm - 1):
            tem = a[i] * alfa[i - 1] + b[i]
            alfa[i] = -c[i] / tem
            bbeta[i] = (d[i] - a[i] * bbeta[i - 1]) / tem
        tem = a[self.lm-1] * alfa[self.lm - 2] + b[self.lm-1]
        ans[self.lm - 1] = (d[self.lm-1] - a[self.lm-1] * bbeta[self.lm - 1]) / tem

        for i in range(1, self.lm):
            j = self.lm - 1 - i
            ans[j] = alfa[j] * ans[j + 1] + bbeta[j]

        return ans
    
        
    def get_gam_and_fro_factor_ns(self,nh):
        """
        Get correct cam and fro variables depending on
        whether Northern or Southern hemispher is 
        considered
        """
        blm = self.ebeta/self.rlamdo
        if nh:
            gam1 = self.gamn
            gam2 = self.gams
            fro1 = self.foan
        else:
            gam1 = self.gams
            gam2 = self.gamn
            fro1 = self.foas
        factor = self.rlamdo(1.-fro1*gam2/(gam2*gam1-blm*blm)))*self.dt/(self.c1*self.dz[0])
        return factor
        
    def coeff(self, wcfac, gam_fro_fac):
        """
        Calculate a, b c coefficient arrays for hemispher
        """
        a = np.zeros(self.lm)
        b = np.zeros(self.lm)
        c = np.zeros(self.lm)
        rkapafac = 2*self.rkapa*self.dt
        
        b[0] = rkapafac/(self.dz[0]*(0.*self.dz[0] + self.dz[1])) #Can the 0.*dz(0) term be dropped here?
        a[0] = 1.-b[0]+gam_and_fro_factor-wcfac/self.dz[0]
        a[1] = -rkapafac/(self.dz[1]**2)+wcfac/self.dz[1]
        c[1] = -rkapafac/(self.dz[1]*(self.dz[1]+self.dz[2]))-wcfac/self.dz[1]
        b[1] = 1. - a[1] - c[1] 
        for i in range(2,self.lm-1):
            a[i] =-rkapafac/(self.dz[i]*(self.dz[i-1]+self.dz[i]))
            c[i] = -rkapafac/(self.dz[i]*(self.dz[i]+self.dz[i+1]))-wcfac/self.dz[i]
            b[i] = 1. - a[i] - c[i] 

        a[self.lm-1] = -rkapafac/(self.dz[self.lm-1]*(self.dz[self.lm-2]+self.dz[self.lm-1]))
        b[self.lm-1] = 1. - a[self.lm-1] + wfac/dz[self.lm-1] #Her var det brukt i selvom vi var utenfor løkka, litt uklart hva som er ment...
        return a,b,c
        
    def setup_ebud2(self, temp1N, temp1S):
        """
        Setting up coefficients and more for the two hemispheres 
        to be redone every timestep
        """

        #Northern hemisphere:
        wcfac = self.w/(sekday*day_year)*(1-0.3*temp1N/self.threstemp)*self.dt
        self.dtrm1n = 1. - self.cpi*wcfac/self.dz[1] -self.beto*self.dt/(self.c1*self.dz[1])
        self.dtmnl2 = wfac*self.cpi/self.dz[self.lm]
        self.an, self.bn, self.cn = coeff(wcfac, self.get_gam_and_fro_factor_ns(self,True))

        #Southern hemisphere:
        wcfac = self.w/(sekday*day_year)*(1-0.3*temp1S/self.threstemp)*self.dt
        self.dtrm1n = 1. - self.cpi*wcfac/self.dz[1] - self.fnso*self.beto*self.dt/(self.c1*self.dz[1])
        self.dtmnl2 = wfac*self.cpi/self.dz[self.lm]
        self.as, self.bs, self.cs = coeff(wcfac, self.get_gam_and_fro_factor_ns(self,False))        

        return

    def setup_ebud(self):
        """
        Setting up energy budget before run
        """
        rho = 1.03
        htcpty = 0.955
        cnvrt = 0.485
        self.c1 = rho*htcpty*cnvrt*100.*sek_day

        fnsa = 1.0 #Can it be something else
        c1fac = self.dt/(self.c1*self.dz[0]) 
        
        blm = self.ebbeta/self.rlamdo
        self.gamn = self.foan + self.rlamda/self.rlamdo + blm
        self.gams = self.foas + self.rlamda/self.rlamdo + fnsa*blm

        #Northern hemisphere
        self.dtrm2n = (self.beto + self.foas*self.ebbeta/(self.gams*self.gamn -fnsa*blm**2))*c1fac
        self.dtrm3n = self.gams/(self.gams*self.gamn-fnsa*blm**2)*c1fac
        self.dtrm4n = blm/(self.gams*self.gamn-fnsa*blm**2)*c1fac

        #Southern hemisphere
        self.dtrm2s=(self.fnso*self.beto+self.foan*fnsa*self.ebbeta/(self.gams*self.gamn-fnsa*blm**2))*c1fac
        self.dtrm3s=self.gamn/(self.gams*self.gamn-fnsa*blm**2)*c1fac
        self.dtrm4s=fnsa*blm/(self.gams*self.gamn-fnsa*blm**2)*c1fac

        self.dtmnl3=self.dt*self.beto/(self.c1*self.dz[lm-1])
        self.dtmnl1=1.-self.dtmnl3
        self.dtmsl3=fnso*self.dtmnl3
        self.dtmsl1=1.-self.dtmsl3

        #Intialising temperature values
        self.tn = np.zeros(self.lm)
        self.tn = np.zeros(self.lm)

        self.betag = 3.e-4
        self.betaa = 2.e-4

        self.z0 = 0.5
        self.betas = 0.25
        self.ebtau = 20.0
        self.zso = 0.0
        self.zgo = 0.0
        self.zao = 0.0

        return
    
    def setup_sealevel_rise(self):
        """
        Setting up variables to be used in sea level rise calculations
        """
        self.press = np.zeros(self.lm)
        self.tempunp = np.zeros(self.lm)
        self.press[0] = 35.0*1.e4*1.e-5
        self.tempunp[0] = 19.5
        self.dens0 = np.zeros(self.lm)
        self.dens0[0] = _density(self.press[0], self.tempunp[0])
        
        for i in range(1,self.lm):
            self.press[i] =(120.+100.*(i-2))*1.e4 #Units=Pa
            self.press[i] =self.press[i]*1.e-5 #Units=bar
            z = 120. + 100.*(i-1)
            self.tempunp[i] = 125.98*z**(-0.45952)
            self.dens0[i]= _density(self.press[i], self.tempunp[i])
    
    def compute_sea_level_rise(self,templ,dtemp, dtempprev):
        deltsl = np.zeros(2)

        #Sea level rise from temperature change
        for i in range(lm):
            dens1 = _density(self.press[i], (templ[i] +self.tempunp[i]))
            deldens = dens1 - self.dens0[i]
            deltsl[0] = deltsl[0] - deldens*self.dz[i]/dens1
        
        #Sea level rise from melting Ice sheets
        #Maybe outdated
        #Also why not use hemispheric temperature change?
        #Greenland
        self.zgo = self.zgo + 1.5*self.betag*(dtemp + dtempprev)/2.
        #Antarctica
        self.zao = self.zao + self.betaa*(dtemp + dtempprev)/2.
        # Small glaciers:
        aa = self.zso + self.z0*self.betas*dtemp/self.ebtau
        bb = 1. + (1.+self.betas*dtemp)/self.ebtau
        self.zso = aa/bb

        deltsl[1] = self.zgo + self.zao + self.zso
            
    def energy_budget(self, years_since_start):
        """
        Doing energy budget calculation for single year
        """
        temp1n = 0.0
        temp1s = 0.0

        tempn=0.0
        temps=0.0
        tempn_air=0.0
        temps_air=0.0
        tempn_sea=0.0
        temps_sea=0.0

        dfrcn=0.0
        dfrcs=0.0
        dfrcg=0.0
        templ = np.zeros(self.lm)

        dtyear = 1./ldtime
        dn = np.zeros(self.lm)
        ds = np.zeros(self.lm)
        for im in range(self.ldtime):
            
            self.setup_ebud2(temp1n, temp1s)

            dqn = im*self.FN*dtyear + (1-im*dtyear)*self.FNOLD + self.FN_VOLC(im)
            dqs = im*FS*dtyear + (1-im*dtyear)*self.FSOLD + self.FS_VOLC(im)
            dn[0]=self.dtrm1n*tn[0]+self.dtrm2n*ts[0]+self.dtrm3n*dqn+self.dtrm4n*dqs
            ds[0]=self.dtrm1s*ts[0]+self.dtrm2s*tn[0]+self.dtrm3s*dqs + self.dtrm4s*dqn

            for i in range(1,lm-1):
                dn[i]=tn[i] + self.beto*self.dt/(self.c1*self.dz[i]))*(ts[i]-*tn[i])
                ds[i]=ts[i] + self.fnso*self.beto*self.dt/(self.c1*self.dz[i]))*(tn[i]-*ts[i])

            dn[lm-1] = self.dtmnl1*self.tn[lm-1]+self.dtmnl2*tn[0]+dtmnl3*ts[lm-1]
            ds[lm-1]=self.dtmsl1*ts[lm-1]+self.dtmsl2*ts[0]+self.dtmsl3*tn[lm-1]

            tn = self.band(an, bn, cn, dn)
            ts = self.band(as, bs, cs, ds)

            temp1n = tn[0]
            temp1s = ts[0]
            for i in range(lm):
                templ[i] = templ[i] + 0.5*(tn[i] + ts[i])/12. #skulle 12 her vært ldtime?

            fnx = self.rlamda + self.foan*self.rlamdo + self.ebbeta
            fsx =self.rlamda + self.foas*self.rlamdo + self.ebbeta
            tempan=dqn+self.foan*self.rlamdo*temp1n+self.ebbeta*(dqs+self.foas*self.rlamdo*temp1s)/fsx
            tempan = teman/(fnx -self.ebbeta**2/fsx)
            tempas=dqs+self.foas*self.rlamdo*temp1s+self.ebbeta*(dqn+self.foan*self.rlamdo*temp1n)/fnx
            tempas=tempas/(fsx-self.ebbeta**2/fnx)
            tmpn=self.foan*temp1n+(1.-self.foan)*tempan
            tmps=self.foas*temp1s+(1.-self.foas)*tempas

            x1=1638.+float(years_since_start)+float(im-1)/12.

            tempn=tempn+tmpn/12.
            dfrcn=dfrcn+dqn/12.
            temps=temps+tmps/12.
            dfrcs=dfrcs+dqs/12.

            tempn_air=tempn_air+tempan/12.
            temps_air=temps_air+tempas/12.
            
            tempn_sea=tempn_sea+temp1n/12.
            temps_sea=temps_sea+temp1s/12           

        dtemp[year_since_start]=(tempn+temps)/2.  #Global temp chg
        dtempnh[year_since_start] = tempn  # Temp chg NH
        dtempsh[year_since_start] = temps  # Temp chg SH

        dtemp_air[year_since_start]=(tempn_air+temps_air)/2.     #Global temp
        dtempnh_air[year_since_start] = tempn_air             #Temp chg. NH
        dtempsh_air[year_since_start] = temps_air             #Temp chg. SH

        dtemp_sea[year_since_start]=(tempn_sea+temps_sea)/2.  #Global temp
        dtempnh_sea[year_since_start] = tempn_sea             #Temp chg. NH
        dtempsh_sea[year_since_start] = temps_sea
        
        RIBN[[year_since_start] = FN-rlamda*tempn
        RIBS[[year_since_start] = FS-rlamda*temps
        RIB[year_since_start] = (RIBN[[year_since_start]+RIBS[[year_since_start])/2. 

        if year_since_start == 0:
            deltsl[year_since_start] = compute_sea_level_rise(templ, dtemp[year_since_start], 0) 
        else:
            deltsl[year_since_start] = compute_sea_level_rise(templ, dtemp[year_since_start], dtemp[year_since_start]) 
        self.FNOLD = FN
        self.FSOLD = FS           
