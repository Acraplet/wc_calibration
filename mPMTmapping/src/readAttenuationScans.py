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

scan_file = './../ScanAttenuation_withoutText.txt'
scan = np.array(rd.read_data3(scan_file))
df = pd.DataFrame()
for i in range(len(scan[0])):
        c = [scan[0][i],scan[1][i],scan[2][i],scan[3][i],scan[4][i],scan[5][i],scan[6][i]]
        row =  pd.Series(data=c, index=['true', 'config', 'reco', 'err', 'Q_thresh', 'initGuessScat', 'FVAL' ], dtype=np.float64)
        df = df.append(row, ignore_index=True)

#Here info about the attenuation
true_att = np.array([50, 60, 80, 100, 125, 150, 200]) * 1/100
reco_att = np.array([48.5453,  64.5826, 81.2833, 97.2022, 121.874, 147.323,208.556 ]) * 1/100
err_att = np.array([0.201828, 2.49692, 0.325958,1.3145, 1.97521, 2.8501, 11.4826]) * 1/100

true_att_binned = np.array([40, 50, 60, 80, 100, 125, 150, 200, 220]) * 1/100
reco_att_binned = np.array([43.3924, 52.4663, 62.7514, 86.1419, 101.487, 140.186, 169.491, 228.738,242.46 ]) * 1/100
err_att_binned = np.array([0.379626, 0.520512, 0.621977, 1.14762, 1.83234, 1.75267, 2.82106, 7.4886, 5.00852]) * 1/100



#df = df[df['spline_increment'] > 11]

df['reco-true/true'] = (df['reco']-df['true'])/df['true']

df = df[abs(df['reco-true/true'])<1.0]
df['true-init/true'] = (df['initGuessScat'])-df['true']/df['true']
df['err/true'] = (df['err'])/df['true']

print(df)

plt.errorbar(df['reco-true/true'], (df['initGuessScat']-df['true'])/df['true'], xerr = df['err/true'], fmt ='x')
plt.xlabel('Initial guess - true/True absorption length')
plt.ylabel('reco-true/true Absorption Length')
plt.title('Absorption length scan initial guess')
plt.show()

df_buf = df
#[abs(df['true-init/true'])<1]

plt.errorbar(df_buf['reco-true/true'],df_buf['FVAL'], fmt ='x')
plt.xlabel('Reco-true/True absorption length')
plt.semilogy()
plt.ylabel('FVAL')
plt.title('Absorption length scan vs FVAL')
plt.show()


df_buf = df[abs(df['FVAL'])<5]


plt.plot([0.4, 5], [0.4, 5], '--', color = 'darkgray')
plt.errorbar(df_buf['true']/100, df_buf['reco']/100, yerr = df_buf['err']/100, fmt ='x', color = 'red', markersize = 10)
plt.xlabel('True absorption length (m)')
# plt.semilogy()
plt.ylabel('Reco absorption length (m)')
plt.title('Absorption length binned and non-binned', weight = 'bold')


plt.errorbar(true_att, reco_att, yerr = err_att, fmt ='x', label = 'Absorption only', markersize = 10)
plt.errorbar(true_att_binned, reco_att_binned, yerr = err_att_binned, fmt ='x', color = 'red', label = 'Absorption only\nbinned method', markersize = 10)
plt.legend()
plt.xlabel('True attenuation length (m)')
# plt.semilogx()
# plt.semilogy()
plt.ylabel('Estimate attenuation length (m)')
plt.grid(which = 'minor', color = 'darkgray')


plt.show()
