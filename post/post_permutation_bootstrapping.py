# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 14:39:38 2015

@author: fbeyer
"""
from scipy.stats.stats import pearsonr
import numpy as np
import nibabel as nb
import os
import matplotlib.pyplot as plt
from numpy import matrix

import pandas as pd
#
import nipype
import nipype.interfaces.fsl as fsl
from nilearn.input_data import NiftiMasker
import nilearn
from scipy.stats.stats import percentileofscore

##
#######ANALYSIS_DIR##############
analysis_dir='/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9/'

#for analysis_dir in np.array(["/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_medication/",\
#                             "/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_wo_crp/",\
#                             "/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_wo_adipo/", \
#                             "/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_wo_IL6/",\
#                             "/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_w_bp/"]):

print analysis_dir
os.chdir(analysis_dir)

#select component
comp=0

#first check significance
#load singular values distribution and plot
os.chdir(analysis_dir)
sing_p=np.load('singvals_p.npy')
print sing_p
print "p-value of salience of comp %i: %.6f" %(comp, sing_p[comp])

sing_vals=np.load('singular_values.npy')
sing_vals_distr=np.load('singular_values_sampled.npy')
fig2=plt.figure()
counts,bin_edges = np.histogram(sing_vals_distr[:,comp],20)

bin_centres = (bin_edges[:-1] + bin_edges[1:])/2.
plt.plot(bin_centres, counts, 'bo')
plt.plot(sing_vals[comp], 50, 'ro')
#plt.plot(sing_vals_distr[:,comp],'bo')
#plt.plot(sing_vals[comp], 'ro')
#plt.hist(sing_vals_distr[:,comp])
#plt.hist(sing_vals[comp])
#plt.xlim(-50)
plt.ylabel('singular value of first component', fontsize=22)
plt.xlabel('# permutations', fontsize=22)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='minor', labelsize=18)
#plt.show()


##inertia/total amount of variance explained
inertia=np.load('inertia.npy')
inertia_sampled=np.load("singular_values_sampled.npy")
#p-value calculated with the script used to be wrong, inertia_p=np.load('inertia_p.npy')
#corrected:
inertia_p_val = (100.0-percentileofscore(inertia_sampled.sum(axis=1), inertia))/100.0
print "inertia calculated based on original and sum of permutations p: %.5d" %(len(inertia_sampled[inertia_sampled>inertia])/np.shape(inertia_sampled)[0])
print "inertia calculated based on original and sum of permutations p: %.5d" %(inertia_p_val)



plt.figure()
plt.hist(inertia_sampled.sum(axis=1))
#plt.show()

#explained variance
fig3=plt.figure()
plt.plot(sing_vals, 'ro')
exVar=sing_vals[comp]/inertia
print "explained variance of comp %i : %.3f" %(comp, exVar)
#plt.show()

print sing_vals/inertia


##transform brain saliencies and bootstraps.

TBM_sal=np.load('brain_saliences.npy')
TBM_bootstrap_ratios=np.load('brain_salience_bootstrap.npy')


#reshape imaging output
#for this recalculate original mask

if (~os.path.isfile(analysis_dir + 'VBM_salience_bootstraps0.nii.gz')):

    nifti_masker = NiftiMasker(standardize=True,
                   memory="nilearn_cache", memory_level=1)
    brain_data = nifti_masker.fit_transform('/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/N749/VBM_input/for_analysis/mean_merged_VBM_raw_N749_mask0.3.nii.gz')

    #VBM mask:'/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/VBM_input/mean_merged_VBM_raw_N749_mask0.3.nii.gz'
    #TBM mask: '/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/mask_for_Jacdet_GM.nii.gz'

    #revert TBM data to original space
    n_components=3
    for i in range(0,n_components): #):
        temp=nifti_masker.inverse_transform(TBM_sal[:,i])
        temp2=nifti_masker.inverse_transform(TBM_bootstrap_ratios[:,i])
        nb.save(temp, analysis_dir +'VBM_salience%s.nii.gz' %i)
        nb.save(temp2, analysis_dir + 'VBM_salience_bootstraps%s.nii.gz' %i)
else: print "nifti-creation complete"
        
########OBESITY DATA##############
#load obesity measures bottstrap and plot
Ob_boo=np.load('Ob_salience_bootstrap.npy')
Ob_boo_abs=abs(Ob_boo)

Ob_sal=np.load('Ob_salience.npy')
Ob_sal_abs=abs(Ob_sal)

print np.shape(Ob_sal)

#Order from spss file:
#RES_BMI
#RES_WHR
#RES_HBA1C
#RES_CHOL
#RES_HDLC
#RES_IL6
#RES_IL6_wo_lowerthreshold
#RES_CRPHS
#RES_ADIPO
#RES_LEP
#RES_sys_bp

if np.shape(Ob_sal)[1]==9:
    Ob_titles=['BMI', 'WHR', 'HBA1C', 'CHOL', 'HDLC', 'IL6','CRPHS',  'Adiponectin','Leptin']
else:
    Ob_titles=['BMI', 'WHR', 'HBA1C', 'CHOL', 'HDLC', 'IL6','CRPHS',  'Adiponectin','Leptin', 'sys']

x=np.arange(0,np.shape(Ob_titles)[0])

#plotting related to obesity measures
#plotting salience ratios
fig1=plt.figure()
plt.xticks(x, Ob_titles, rotation=90)
plt.bar(x,Ob_sal[:,comp])
#
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on'
    ) # labels along the bottom edge are off
#plt.rc('ytick', labelsize=16)
#plt.rc('ytick', labelsize=16)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.ylabel('Ob saliencies/contributions to latent variable', fontsize=20)
plt.ylim((-0.5,0.5))

print "saliencies of component 0:" 
print Ob_sal[:,0]
print "saliencies of component 1:" 
print Ob_sal[:,1]

#plotting bootstrap ratios
fig2=plt.figure()
plt.xticks(x, Ob_titles, rotation=90)
plt.bar(x,Ob_boo[:,comp])
#
plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='on') # labels along the bottom edge are off
#plt.rc('ytick', labelsize=12)
#plt.rc('ytick', labelsize=12)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.ylabel('Z-scored Ob effect from bootstrapping', fontsize=20)
plt.ylim((-2.6,2.6))
sig_limu=np.ones(len(Ob_titles))*-2.3
plt.plot(sig_limu, 'r')
sig_limup=np.ones(len(Ob_titles))*2.3
plt.plot(sig_limup, 'r')

#plt.show()
##
#
print "bootstraps of component 0:" 
print Ob_boo[:,0]
print "bootstraps of component 1:" 
print Ob_boo[:,1]

#####CALCULATE LATENT VARIABLES#########
#
#calculate latent variables:
brain=np.load("Xscaled.npy")
brain_sal=np.load("brain_saliences.npy")

ob_sal=np.load("Ob_salience.npy")
ob=np.load("Yscaled.npy")



lv_brain=matrix(brain)*matrix(brain_sal)
#plt.figure()
plt.plot(np.arange(0,len(lv_brain[:,comp])),lv_brain[:,comp])
#plt.show()

lv_ob=matrix(ob)*matrix(ob_sal)
#plt.figure()
plt.plot(np.arange(0,len(lv_ob[:,comp])),lv_ob[:,comp])
#plt.show()


corr=pearsonr(lv_ob[:,comp],lv_brain[:,comp])
print corr

#plt.figure()
plt.plot(lv_ob[:,comp],lv_brain[:,comp],'x')
#plt.show()
#

df = pd.read_csv('/afs/cbs.mpg.de/share/gr_agingandobesity/life_shared/Frauke/metabolic_VBM/748_allLabvals_MMST_T1_qual_checked_log_transforms_all_regressed.csv')
#
df["lv_ob"]=lv_ob[:,comp]
df["lv_brain"]=lv_brain[:,comp]
#
df.to_csv('/afs/cbs.mpg.de/share/gr_agingandobesity/life_shared/Frauke/metabolic_VBM/748_allLabvals_MMST_T1_qual_checked_log_transforms_all_regressed_with_LV.csv')
