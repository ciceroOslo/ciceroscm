import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

pertfilename = 'pertem_test.txt'

#Read components, to get the units:
inputdir = '/div/no-backup/users/ragnhibs/ciceroscm/tests/test-data/'
df_comp =pd.read_csv(inputdir + 'gases_v1RCMIP.txt', sep='\t', index_col=0)

#Choose component to perturbe:
comp = 'CO2'
unit = df_comp.loc[comp]['EM_UNIT']
#Should add a check if components are not included.
#But anyway it will crash here

#Make sure to make the perturbation in the correct unit.

#Can read these from file
pert_years = [2000,2001,2002]
#The emissions will be multiplied by -1 to be subtracted from the file
emissions = np.array([4,4,4])


columns = ['component','unit','source','year','emission']

df = pd.DataFrame(columns=columns)
df['year']=pert_years
df['emission']=emissions*-1.0
df['source'] = 'src'
df['component'] = comp
df['unit'] = unit
df = df.set_index('component')

print(df)


df.to_csv(pertfilename)


