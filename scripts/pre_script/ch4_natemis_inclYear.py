import pandas as pd
import numpy as np

fileDir = "/div/amoc/CSCM/SCM_Linux_v2019/Matlab_pre/NatEmis/"

xl = pd.ExcelFile(fileDir+"em_conc_SCM_vAug2017.xlsx")
sheet_names = xl.sheet_names
sheet_data = xl.parse(sheet_names[1])
df_conc = sheet_data.drop(axis=0,index=[0,1,2])

years = df_conc["Component"].values
conc_CH4 = df_conc["CH4"].values
conc_N2O = df_conc["N2O"].values

xl = pd.ExcelFile(fileDir+"emissions_RCMIP.xlsx")
sheet_names = xl.sheet_names
sheet_data = xl.parse(sheet_names[0])
df_em = sheet_data

years_em = sheet_data["Year"].values
em_CH4 = sheet_data["CH4"].values
em_N2O = sheet_data["N2O"].values
em_NOX = sheet_data["NOx"].values
em_CO = sheet_data["CO"].values
em_NMVOC = sheet_data["NMVOC"].values

startyr = 1751
ax = np.where(years==1751)[0][0]

TAU_1 = 9.6
TAU_2 = 120
TAU_3 = 160

BETA_CH4 = 2.78

TAU_1_N2O = 121

em_nat_ch4 = 275
em_nat_n2o = 9.5
em_nat_ch4_org = 275 
em_nat_n2o_org = 9.5 

em_nat_hist = np.full(len(years), em_nat_ch4,np.double)
IDTM = 24

CONC_NEW_last = conc_CH4[0]
for i in range(ax,len(years)):
    not_done = True
    while not_done:
        dlnOH =  -0.32 \
        * (np.log(CONC_NEW_last) \
        - np.log(1751.0) \
        ) + 0.0042 * (em_NOX[i]  
        -  em_NOX[np.where(years==2000)[0][0]] 
        ) - 0.000105 * (em_CO[i] 
        - em_CO[np.where(years==2000)[0][0]] 
        ) - 0.000315 * (em_NMVOC[i] 
        - em_NMVOC[np.where(years==2000)[0][0]])
        
        q = ((1.0/TAU_1)*(dlnOH + 1) + 1/TAU_2 + 1/TAU_3)/IDTM
        pc = (em_CH4[i] + em_nat_ch4)/(IDTM*BETA_CH4)

        CONC_NEW = CONC_NEW_last
        for _ in range(0,IDTM):  
            CONC_NEW = (pc/q+(CONC_NEW-pc/q)*np.exp(-q*1.0))
     
        if (CONC_NEW/conc_CH4[i] >  1.005): 
            em_nat_ch4 = em_nat_ch4*0.995
        elif (CONC_NEW/conc_CH4[i] < 0.995): 
            em_nat_ch4 = em_nat_ch4*1.005
        else:
            em_nat_hist[i] = em_nat_ch4
            not_done = False
    CONC_NEW_last = CONC_NEW


year_out = range(1750,2501)
em_nat_out = np.concatenate((em_nat_hist,  np.full((len(year_out)-len(em_nat_hist)), np.mean(em_nat_hist[-11:]),np.double)), axis=0)

df = pd.DataFrame(columns=["year","em_nat_out"])
df["year"] = year_out
df["em_nat_out"] = em_nat_out
df.to_csv("natemis_ch4_rcmip_inclYear.txt",index=False)
