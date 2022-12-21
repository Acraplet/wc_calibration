#This is a small programm to chekc how the charge evolves with R at given theta, phi positions for the same scattering length
import read_data_results3 as rd
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd


plt.style.use(["science", "notebook", "grid"])
df = pd.DataFrame()
#First read all the data
for name in sys.argv[1:]:
    data = np.array(rd.read_data3("%s"%name)).T
    for i in range(len(data)): #these are the number of lines in the text which correspond to
        #the number of abwff, rayff data points we have
#The store the attenuation, scatterning, position and R
        abwff = data[i][8]
        rayff = data[i][9]
        theta = data[i][3]
        phi = data[i][4]
        R = data[i][5]
        Q = data[i][6]
        alpha = 0
        if rayff == 0.000555:
            alpha = 5
        if rayff == 0.001111:
            alpha = 10
        if rayff == 0.002222:
            alpha = 20
        if rayff == 0.003333:
            alpha = 40
        if rayff == 0.00555:
            alpha = 60
        if rayff == 0.008325:
            alpha = 100
        if rayff == 0.01221:
            alpha = 220
        if rayff == 0.00445:
            alpha = 80
        if rayff == 0.0067:
            alpha = 120
        if rayff == 0.00945:
            alpha = 170
        if rayff == 0.0106:
            alpha = 190
        if abwff == 0.000486:
            alpha = 20
        if abwff == 0.000972:
            alpha = 40
        if abwff == 0.001458:
            alpha = 60
        if abwff == 0.00243:
            alpha = 100
        if abwff == 0.000243:
            alpha = 10
        if abwff == 0.003645:
            alpha = 150
        if abwff == 0.005346:
            alpha = 220



        c = [abwff, rayff, theta, phi, R, Q, alpha]
        row = pd.Series(data=c, index=['abwff', 'rayff', 'theta', 'phi', 'R', 'Q', 'alpha'], dtype=np.float64)
        df = df.append(row, ignore_index=True)

df = df.sort_values(["R", 'abwff', 'rayff'], axis = 0, ascending = True).reset_index()

#Then plot
i = 0
for t in df['theta'].unique():
    if t != 0.:
        for p in df['phi'].unique():
            plt.figure(figsize = (16, 10))
            for a in df['abwff'].unique():
                if a <= 1e4:
                    df_buf = df[df['theta'] == t]
                    df_buf = df_buf[df_buf['phi'] == p]
                    df_buf = df_buf[df_buf['abwff'] == a]
                    #df_buf = df_buf[df_buf['rayff'] == r]
                    plt.subplot(3,1,1)
                    #print(a)
                    plt.plot(df_buf['R'], df_buf['Q'], 'x-', label = 'abwff = %.2e, alpha = %icm'%(a, df_buf['alpha'].unique()))
            for r in df['rayff'].unique():
                if r <=1e4 and i%2 == 0:
                    df_buf = df[df['theta'] == t]
                    df_buf = df_buf[df_buf['phi'] == p]
                    #df_buf = df_buf[df_buf['abwff'] == a]
                    df_buf = df_buf[df_buf['rayff'] == r]
                    plt.subplot(3,1,2)
                    plt.plot(df_buf['R'], df_buf['Q'], 'x-', label = 'rayff = %.2e, alpha = %icm'%(r, df_buf['alpha'].unique()) )


                if r <=1e4 and i%2 == 1:
                    df_buf = df[df['theta'] == t]
                    df_buf = df_buf[df_buf['phi'] == p]
                    #df_buf = df_buf[df_buf['abwff'] == a]
                    df_buf = df_buf[df_buf['rayff'] == r]
                    plt.subplot(3,1,3)
                    plt.plot(df_buf['R'], df_buf['Q'], 'x-', label = 'rayff = %.2e, alpha = %icm'%(r, df_buf['alpha'].unique()) )

                i += 1


            plt.subplot(3,1,1)
            plt.title('Position: theta %.2f phi %.2f'%(t, p))
            plt.legend(fontsize = 12)
            plt.ylabel('total charge collected\nper 1000 photons ')
            plt.subplot(3,1,3)
            plt.xlabel('R mPMT dome-source distance (cm) ')
            plt.ylabel('total charge collected\nper 1000 photons ')
            plt.legend(fontsize = 12)
            plt.subplot(3,1,2)
            plt.legend(fontsize = 12)
            plt.ylabel('total charge collected\nper 1000 photons ')

            plt.savefig("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Maps/OnePosition_charge_versus_R/OnePosition_theta%.2f_phi%.2f.pdf"%(t, p))
        #plt.show()
