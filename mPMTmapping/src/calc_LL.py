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

#plt.style.use(["science", "notebook", "grid"])

def binomial(ki,pi, e = 1000):
    #probability to obtain ki hits in bin i which has probability to be hit pi
    return scipy.stats.binom.pmf(ki, e, pi)

def PoissonLL(ni, ki, e): #too large - not working
    return np.exp(-ni) * 1/math.factorial(ki) * ni**(ki)

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
        e_ref = 1000
        #print(observed * (e_ref/e) * np.log(observed * (e_ref/e)) + observed * (e_ref/e))
        return observed * (e_ref/e) * np.log(observed * (e_ref/e)) + observed * (e_ref/e)

    else:
        e_ref2= 1000
        #print(expected, observed * (e_ref/e))
        #print(np.log((expected - observed * (e_ref/e)) ** 2 / expected))
        return (expected * (e_ref2/e_ref) - observed * (e_ref2/e)) ** 2 / (expected * (e_ref2/e_ref))

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

opts, args = getopt.getopt(argv, "m:o:p:f:c:a:")
for opt, arg in opts:
        if opt in ['-m']:
            phi_max = float(arg)
        elif opt in ['-f']:
            filename = str(arg)
        elif opt in ['-c']:
            comparision_file = str(arg)
        elif opt in ['-o']:
            outputfile_name = str(arg)


df = make_df(filename)
df_test = make_df(comparision_file)

probability = 0.0

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
        e = int(df_test_buf2['events'])
        e_ref = int(df_buf2['events'])
        #
        #binom = binomial(ki, pi, e)

        probability += LLchi2(ni, ki, e_ref, e)



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
        #multiplicating the porbablility in each bin to have the total one

print('\nRef_file: %s, compa file: %s, L(chi2): %.1f\n' %(filename, comparision_file, probability))


#will have to map the positions properly (down the line we will have a continuous map - much better)

#Import the maps



#do the calculation
#return the proba - maybe a log, no idea - do a log?

#That mutiplication is basically the likelihood - work with the log likelihood



