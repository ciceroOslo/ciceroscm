import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pertfilename = 'pertforc_test.txt'

#Read components, to get the units:
inputdir = '/div/no-backup/users/ragnhibs/ciceroscm/tests/test-data/'
df_comp =pd.read_csv(inputdir + 'gases_v1RCMIP.txt', sep='\t', index_col=0)

print(df_comp)

#Choose component to perturbe:
comp = 'OTHER'
unit_test = df_comp.loc[comp]['EM_UNIT']
#Should add a check if components are not included.
#But anyway it will crash here. We do not need the unit

unit ='Wm-2'

#Make sure to make the perturbation in the correct unit.

#Can read these from file
pert_years = [2000,2001,2002]
#The forcing will be multiplied by -1 to be subtracted from the file
forcing = np.array([4,4,4])


columns = ['component','unit','source','year','forcing']

df = pd.DataFrame(columns=columns)
df['year']=pert_years
df['forcing']=forcing*-1.0
df['source'] = 'src'
df['component'] = comp
df['unit'] = unit
df = df.set_index('component')

print(df)


df.to_csv(pertfilename)


