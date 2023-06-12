#This is a small coed to read and plot the output of
#The scan as a function of different variables

import read_data_results3 as rd
#read the .txt file:
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import getopt
import pandas as pd

plt.style.use(["science", "notebook", "grid"])
R_list = [10, 20, 40, 80, 120, 140, 160,180,210,250]
R_list_n = [10, 20, 40, 80, 120, 140, 160,180,210,250, 0]

scan = []

for R in R_list:
        scan_file = './Maps/2D_ref_maps/absorption_ref_file_R%.2f.txt'%R
        scan_val = np.array(rd.read_data3(scan_file))

        # print(scan_val[0])
        scan.append(scan_val[0])

scan.append(np.arange(0, len(scan_val[0]), 1))

print(scan)
# scan_file = './../InitGuessScanScatteringLength_ManyR_withoutText.txt'
# scan = np.array(rd.read_data3(scan_file))
df = pd.DataFrame()
for i in range(len(scan[0])):
        c = [scan[0][i],scan[1][i],scan[2][i],scan[3][i],scan[4][i],scan[5][i],scan[6][i],scan[7][i],scan[8][i], scan[9][i], scan[10][i]]
        row =  pd.Series(data=c, index=R_list_n, dtype=np.float64)
        df = df.append(row, ignore_index=True)

print(df)
#df = df.drop(df[df[10] < 0.001].index)

binTarget = 1
binTarget2 = 6
df1 = df.drop(df[df[0] != binTarget].index)
df2 = df.drop(df[df[0] != binTarget2].index)

for R in R_list:

        if i!=0:
                plt.plot(R, df1[R], 'bx')
                plt.plot(R, df2[R], 'rx')
                i = 0

        if i==0:
                i = 1

plt.plot(R, df1[R], 'bx',  label = 'Bin %s'%binTarget)
plt.plot(R, df2[R], 'rx',  label = 'Bin %s'%binTarget2)
plt.legend()
plt.xlabel('R (cm)')
plt.ylabel('Number of pe recorded per photon')
plt.title('Study of the stability of the R=0 approximation\n in the binned approach')
plt.show()

raise end

#now the plotting
df = df[df['spline_increment'] > 11]

df['true-reco/true'] = (df['true']-df['reco'])/df['true']
df['true-init/true'] = (df['true']-df['initGuessScat'])/df['true']
df['err/true'] = (df['err'])/df['true']

print(df)

plt.errorbar((df['true']-df['initGuessScat'])/df['true'],df['true-reco/true'], yerr = df['err/true'], fmt ='x')
plt.xlabel('True-initialGuess/True scattering length')
plt.ylabel('true-reco/true Scattering Length')
plt.title('Scattering length scan initial guess')
plt.show()

df_buf = df
#[abs(df['true-init/true'])<1]

plt.errorbar(df_buf['true-reco/true'],df_buf['FVAL'], fmt ='x')
plt.xlabel('True-reco/True scattering length')
plt.semilogy()
plt.ylabel('FVAL')
plt.title('Scattering length scan vs FVAL')
plt.show()


df_buf = df[abs(df['FVAL'])<5]


plt.plot([0.4, 50], [0.4, 50], '--', color = 'darkgray')
plt.errorbar(df_buf['true']/100, df_buf['reco']/100, yerr = df_buf['err']/100, fmt ='x',label = 'Scattering only', markersize = 10)
plt.xlabel('True scattering length (m)')
# plt.semilogy()
plt.ylabel('Reco scattering length (m)')
plt.title('Scattering length estimation 2D binned (FVAL < 5) \n and absorption length estimation 1D non-binned', weight = 'bold')


plt.errorbar(true_att, reco_att, yerr = err_att, fmt ='x', label = 'Absorption only', markersize = 10)
plt.errorbar(true_att_binned, reco_att_binned, yerr = err_att_binned, fmt ='x', label = 'Absorption only\nbinned method', markersize = 10)
plt.legend()
plt.xlabel('True attenuation length (m)')
# plt.semilogx()
# plt.semilogy()
plt.ylabel('Estimate attenuation length (m)')
plt.grid(which = 'minor', color = 'darkgray')


plt.show()
