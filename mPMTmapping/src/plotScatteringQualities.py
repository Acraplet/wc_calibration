#Read a tsxt file
import read_data_results3 as rd
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import getopt
import pandas as pd

test_or_ref = "ref"
spline_or_5nodesRef = "5nodesRef"
folder2 = "Ref_files"
if test_or_ref == "test":
    folder2 = "Test_files"
folder1 = "Bin24_Scattering_FitQuality_5nodes_ref"
if spline_or_5nodesRef == "SimpleSpline":
    folder1 = "Bin24_Scattering_FitQuality_SimpleSpline"

plt.style.use(["science", "notebook", "grid"])

df = pd.DataFrame()

for test_or_ref in ["test", "ref"]:
    for spline_or_5nodesRef in ["SimpleSpline"]: #for now only simple spline

        folder1 = "Bin24_Scattering_FitQuality_SimpleSpline"
        folder2 = "All_files"

        PMT_position_file = "/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Bin24_ScatteringLengthEstimation_%s_%sRuns.txt"%(spline_or_5nodesRef, test_or_ref)
        print(PMT_position_file)

        all_data  = np.array(rd.read_data3(PMT_position_file))

        for i in range(len(all_data[0])):
            t = 0
            if test_or_ref == "test":
                t = 1
            c = [all_data[0][i],all_data[1][i],all_data[2][i],all_data[3][i], t]
            row =  pd.Series(data=c, index=['alpha_true', 'config', 'alpha_reco', 'err', 'test'], dtype=np.float64)
            df = df.append(row, ignore_index=True)



print(df['alpha_true'].unique())

#here we are plotting everything as a comparision with distances
for a in df['alpha_true'].unique():
    buf = df[df['alpha_true'] == a]
    labels = []
    for c in buf['config']:
        if c == 3210:
            labels.append("close only")
        if c == 6543:
            labels.append("middle only")
        if c == 9876:
            labels.append("far only")
        if c == 9876543:
            labels.append("middle & far")
        if c == 6543210:
            labels.append("middle & close")
        if c == 9876543210:
            labels.append("far, close & middle")
    plt.figure(figsize=(15, 20))

    plt.subplot(2,1,1)
    plt.title("Estimation of the Scattering length - %s dataset with %s fit\nTrue Rayleigh scattering length: %.2fcm" %(test_or_ref, spline_or_5nodesRef,  a))
    plt.errorbar(labels, buf['alpha_reco'], yerr = buf["err"], fmt = 'x',label = '%s dataset fitted with %s'%(test_or_ref, spline_or_5nodesRef))
    plt.plot(labels, buf['alpha_true'], 'k-', label = 'true scattering length: %.2f'%(a))
    plt.ylabel("reco scattering length")
    plt.legend()
    plt.subplot(2,1,2)
    plt.errorbar(labels, (-buf['alpha_true']+buf['alpha_reco'])/buf['alpha_true'], yerr = buf['err']/buf['alpha_true'], fmt = 'x')
    plt.plot(labels, buf['alpha_true'] * 0, 'k-')

    plt.xticks(rotation = 0)
    plt.xlabel("Distances used in the fit")
    plt.ylabel("reco-true/true scattering length")
    #plt.savefig("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Pictures/%s/ScatteringLengthEstimations/%s/ScatteringlengthEst_%sFile_%s_TrueAtt%.2f.png"%(folder1, folder2, test_or_ref, spline_or_5nodesRef, a))
    plt.close()
    #plt.show()

#Here we are plotting all the reconstructions but when we have the same set of distances
for c in df['config'].unique():
    buf1 = df[df['config'] == c]
    plt.figure(figsize = (20, 10))
    plt.subplot(1,1,1)
    if c == 3210:
            labels = "close only"
    if c == 6543:
        labels = "middle only"
    if c == 9876:
        labels = "far only"
    if c == 9876543:
        labels = "middle & far"
    if c == 6543210:
        labels = "middle & close"
    if c == 9876543210:
        labels = "far, close & middle"
    test = 0 #to ony have the legend once
    ref = 0
    plt.title("Scattering length estimations with distances used: %s\nprofile fitting method: %s"% (labels, spline_or_5nodesRef))
    plt.plot([0, 6000], [0,0], '--', color = 'darkgray', label = 'Scat. length correctly estimated')
    for a in buf1['alpha_true']:

        buf = buf1[buf1['alpha_true'] == a]
        if a >= 5500:
            buf['err'] = buf['err']/10

        if buf["test"].unique() == 1:
            if test == 1:
                plt.errorbar(buf['alpha_true'], (buf['alpha_reco']-buf['alpha_true'])/buf['alpha_true'], yerr = buf["err"]/buf['alpha_true'], color = "k", fmt = 'x')
            if test != 1:
                plt.errorbar(buf['alpha_true'], (buf['alpha_reco']-buf['alpha_true'])/buf['alpha_true'], yerr = buf["err"]/buf['alpha_true'], color = "k", fmt = 'x',label = 'Test dataset')
                test = 1
        if buf["test"].unique() == 0:
            if ref != 1:
                plt.errorbar(buf['alpha_true'], (buf['alpha_reco']-buf['alpha_true'])/buf['alpha_true'], yerr = buf["err"]/buf['alpha_true'], color = "r", fmt = 'x',label = 'Ref dataset')
                ref = 1
            if ref == 1:
                plt.errorbar(buf['alpha_true'], (buf['alpha_reco']-buf['alpha_true'])/buf['alpha_true'], yerr = buf["err"]/buf['alpha_true'], color = "r", fmt = 'x')
            if a >= 5500:
                plt.errorbar(buf['alpha_true'], (buf['alpha_reco']-buf['alpha_true'])/buf['alpha_true'], yerr = buf["err"]/buf['alpha_true'], color = "b", fmt = 'x', label = 'Ref dataset, error divided by 10')


    plt.legend(loc = 2)
    plt.xlabel('True scattering length (cm)')
    plt.ylabel('estimated - true / true scattering length')
    plt.savefig("/home/ac4317/Laptops/Year1/WCTE/wc_calibration/mPMTmapping/Pictures/%s/ScatteringLengthEstimations/%s/ScatteringlengthEst_Ratio_%sFile_%s_config%i_ErrShrunk.png"%(folder1, folder2, "all", spline_or_5nodesRef, c))
    #plt.show()

#plot and save up, for every given true att length

