#This is a small test script to try to understand how I am going to do my likelihood tets

#For each bin I have a probability of hit given by my efficiency map
#For a given trial that I will call a 'test' going forward I can calculate the overall probability that we have to obtain that specific test if the true parameter value is the one given by map x

import scipy as sp
import read_data_results3 as rd
#read the .txt file:
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import getopt
import pandas as pd
import scipy.stats
import math

filename = "No filename, please input one with the -f option"
comparision_file = 'none'
outputfile_name = 'none'

plt.style.use(["science", "notebook", "grid"])

def binomial(ki,pi, e = 1000):
    #probability to obtain ki hits in bin i which has probability to be hit pi
    return scipy.stats.binom.pmf(ki, e, pi)

def PoissonLL(ni, ki, e_ref, e): #negative values when ki >> ni not working
    ni = ni * e/e_ref
    #print(ni)
    if ki != 0 and ni!=0:

        a = - ni + ki
        k = ni/ki
        b = ki * np.log(k)
        #print("%.1f %.1f %.2f"%(ni, ki, a + b))
        return ni - ki + ki * np.log(ni/ki)

    if ki == 0 and ni!=0:
        #print("%.1f %.1f %.2f"%(ni, ki, ni))
        return ni # this shouldn't be -ni?
    if ni == 0 and ki!=0:
        #print("%.1f %.1f %.2f"%(ni, ki, ki * np.log(ki) + ki))
        return ki * np.log(ki) + ki
    else:
        return 0

def GaussianLL(ni, ki, e):
    return a*np.sqrt(e)/np.sqrt(2*np.pi) * np.exp(-((ki-ni)**2)/(2))# * ni**(ki)

#def chi2(ki, pi, e):
    #if (pi*e) == 0:
        #return 0
    #else:
        #return (pi*e - ki)**2 / (pi*e)

def LLchi2(expected, observed, e_ref, e):
    if expected == observed * (e_ref/e):
        return 0 #don't penalise if we have the correct thing!
    elif expected == 0:
        observed =  observed * (e_ref/e)
        #print(observed * (e_ref/e) * np.log(observed * (e_ref/e)) + observed * (e_ref/e))
        return observed * np.log(observed) + observed
    else:
        observed =  observed * (e_ref/e)
        #print(expected, observed * (e_ref/e))
        #print(np.log((expected - observed * (e_ref/e)) ** 2 / expected))
        return (expected - observed) ** 2 / (expected)

def make_df(file_name):
    df = pd.DataFrame()
    table = rd.read_data3(file_name)
    x, y, z = np.array(table[0]), np.array(table[1]), np.array(table[2])
    R, theta, phi = np.array(table[5]), np.array(table[3]), np.array(table[4])
    Q_tot = np.array(table[6])
    nEvents = np.array(table[7])
    for i in range(len(x)):
        c = [x[i], y[i],  z[i],  theta[i], phi[i],  R[i],  Q_tot[i],  nEvents[i]]
        row =  pd.Series(data=c, index=['x', 'y', 'z', 'theta', 'phi', 'R', 'Q', 'events'], dtype=np.float64)
        df = df.append(row, ignore_index=True)
    return df

#filename = str(sys.argv[1])
phi_max = 90 #for plotting - get the max phi value to show only the right portion of the circle
argv = sys.argv[1:]
factor_needed = False
opts, args = getopt.getopt(argv, "m:o:p:f:c:a:n:")
for opt, arg in opts:
        if opt in ['-m']:
            phi_max = float(arg)
        elif opt in ['-f']:
            ID_ref = str(arg)
        elif opt in ['-c']:
            ID_compa = str(arg)
        elif opt in ['-o']:
            outputfile_name = str(arg)
        elif opt in ['-a']:
            factor_needed = True
            alpha_test = float(arg)

filename = "maps_txtFiles/mPMT_map_ID%s.txt"%ID_ref
comparision_file = "maps_txtFiles/mPMT_map_ID%s.txt"%ID_compa
df = make_df(filename)
df_test = make_df(comparision_file)

df_test = df_test.sort_values(by='theta')

probability = 0.0
probability_f = 0.0
i = 0
s = []
s_f = []
plt.figure(1,figsize = (20,10))
plt.subplot(211)
plt.figure(2,figsize = (20,10))
plt.subplot(211)



for theta in df_test['theta'].unique(): # theta bins
    df_buf = df[df['theta']==theta] #select the right reference bin
    df_test_buf = df_test[df_test['theta']==theta] #select the right test bin

    for phi in df_test_buf['phi'].unique(): #phi bins
        df_buf2 = df_buf[df_buf['phi']==phi]

        df_test_buf2 = df_test_buf[df_test_buf['phi']==phi]
        #print(df_buf2, "\n",phi, theta, "\n", df_test_buf2, '\n')
        pi = float(df_buf2['Q']/df_buf2['events'])
        ni = float(df_buf2['Q'])
        ki = float(df_test_buf2['Q'])
        if (factor_needed):
            factor = float(1/np.exp(-df_buf2['R']/alpha_test))
        else:
            factor = 1
        ki_f = float(df_test_buf2['Q']) * factor
        e = int(df_test_buf2['events'])
        e_ref = int(df_buf2['events'])
        #
        #binom = binomial(ki, pi, e)

        probability += LLchi2(ni, ki, e_ref, e)
        probability_f += LLchi2(ni, ki_f, e_ref, e)
        plt.figure(1)
        plt.subplot(211)
        plt.plot(i, ki/e, 'kx')
        plt.plot(i, ni/e_ref, 'rx')
        plt.subplot(212)
        if (ni/e_ref) != 0 and (ki) != 0 :
            plt.plot(i, (ki/e)/(ni/e_ref), 'kx')
            s.append((ki/e)/(ni/e_ref))
            a = (ki/e)/(ni/e_ref) #to avoid a zero when plotting the legend
        i = i+1   

        plt.figure(2)
        plt.subplot(211)
        plt.plot(i, ki_f/e, 'kx')
        plt.plot(i, ni/e_ref, 'rx')
        plt.subplot(212)
        if (ni/e_ref) != 0 and (ki_f) != 0 :
            plt.plot(i, (ki_f/e)/(ni/e_ref), 'kx')
            s_f.append((ki_f/e)/(ni/e_ref))
            a_f = (ki_f/e)/(ni/e_ref) #to avoid a zero when plotting the legend
        i = i+1   
        #if binom == 0:
            #print(ki, ni, pi, pi*e, e, e_ref, PoissonLL(int(ni/e_ref * e), int(ki), e))
        #print(binomial(ki, pi, e), ni/e_ref, ki/e)
        #print(chi2(ki, pi, e), ki, pi)
        #if both expected and true are non-zero -> ignore
        #if binom == 0:
            #if ki!=0 and ni==0:
                #probability += (-ki*np.log(ki)-ki)
            #elif ki==0 and ni!=0:
                #probability += (-ni*np.log(ni)-ni)
            #else:
                #print('pi:', pi, 'ki/e:', ki/e, 'ni:', ni, 'ki:', ki, GaussianLL(ni/e_ref * e, ki, e))
        #else:
            #probability += np.log(binom)
        #if (ki < 10 or ni < 10):
            #binom = PoissonLL(int(ni/e_ref * e), int(ki), e)
            #print('Poisson', binom, ni, ki)

        #print(binomial(ki, pi, e), ni/e_ref * e, ki,binom, '\n')
        #if pi!=0:
            #probability += np.log(binom)
        #elif (pi==0 and ki==0):
            #probability += factor_0
        #else:
            #probability += factor_non0 * ki

        #print(GaussianLL(ni/e_ref * e, ki, e), ni/e_ref * e, ki, binomial(ki, pi, e))
        #multiplicatplt.subplot(212)ing the porbablility in each bin to have the total one
plt.figure(1)
plt.subplot(211)
plt.plot(i, ki/e, 'kx', label = 'test Q: %s'%ID_compa)
plt.plot(i, ni/e_ref, 'rx', label = 'reference Q: %s'%ID_ref)
plt.xlabel('Arbitrary source position')
plt.ylabel('Mean charge per photon')
plt.title('Comparision ref: FileID%s compa: FileID%s - total Chi2 %.2f'%(ID_ref, ID_compa, probability))
plt.legend()

plt.subplot(212)
plt.plot(i, a, 'kx', label = 'Mean of exp/obs: %.3f std: %.3f' %(np.array(s).mean(), np.array(s).std()))
plt.xlabel('Arbitrary source position')
plt.ylabel('reference/test mean charge per photon')
plt.legend()
plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Comparisions/Individual_source_positions/Ref%s_Compa%s'%(ID_ref, ID_compa))


plt.figure(2)
plt.subplot(211)
plt.plot(i, ki_f/e, 'kx', label = 'test Q: %s with correction factor: %.3f'%(ID_compa, factor))
plt.plot(i, ni/e_ref, 'rx', label = 'reference Q: %s'%ID_ref)
plt.xlabel('Arbitrary source position')
plt.ylabel('Mean charge per photon')
plt.title('Comparision ref: FileID%s *corrected* compa: FileID%s - total Chi2 %.2f'%(ID_ref, ID_compa, probability_f))
plt.legend()

plt.subplot(212)
plt.plot(i, a_f, 'kx', label = 'Mean of exp/obs_corrected: %.3f std: %.3f' %(np.array(s_f).mean(), np.array(s_f).std()))
plt.xlabel('Arbitrary source position')
plt.ylabel('reference/test mean charge per photon')
plt.legend()
plt.savefig('/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/maps_pictures/Comparisions/Individual_source_positions_corrected/Ref%s_Compa%s_Abwffcorrected'%(ID_ref, ID_compa))
#print('%s %s %.2f\n' %(ID_ref, ID_compa, probability))

print('%s %s %.2f %.1f\n' %(ID_ref, ID_compa, probability, probability_f))


#will have to map the positions properly (down the line we will have a continuous map - much better)

#Import the maps



#do the calculation
#return the proba - maybe a log, no idea - do a log?

#That mutiplication is basically the likelihood - work with the log likelihood



