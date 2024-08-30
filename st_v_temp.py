import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# read lammps thermo dump without having to specify names
def read_dat(filename):
    df = pd.read_csv(filename, header = 0,  delim_whitespace=True)
    df_modified = df[df.columns[:-1]] # get rid of the leading #
    df_modified.columns = df.columns[1:] # reset the columns
    return df_modified


def select_second_half(df):
    half = int((df.Step.max()) * 0.9)
    df = df[df.Step>half]
    return df

T_list = np.arange(1.0, 3.1, 0.2)

st_list = []

for T in T_list:
    filename = f"T{T:.1f}/NVT/extracted_data.dat"
    df = select_second_half(read_dat(filename))
    st = -T * np.log(df.exp.mean())/df.DeltaA.mean()
    st_list.append(st)

plt.plot(T_list, st_list)
plt.show()

data = np.column_stack((T_list, st_list))
np.savetxt('T_st_data.dat', data, header='Temp st')

