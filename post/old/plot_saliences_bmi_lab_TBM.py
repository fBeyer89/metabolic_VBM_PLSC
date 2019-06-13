# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:32:14 2015

@author: fbeyer
"""
import scipy as sp
import numpy as np
import nibabel as nb
import os
import pandas as pd
import matplotlib.pyplot as plt


analysis_dir='/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/perm_boot_logtr_all_regressed_n3_medication_corrected//'

#
os.chdir(analysis_dir)
#choose component to display
comp=0

#load obesity measures bottstrap and plot
Ob_boo=np.load('Ob_salience_bootstrap.npy')
Ob_boo_abs=abs(Ob_boo)

Ob_sal=np.load('Ob_salience.npy')
Ob_sal_abs=abs(Ob_sal)


Ob_titles=['BMI', 'WHR', 'HBA1C', 'CHOL', 'HDLC', 'CRPHS', 'IL6', 'Adiponectin','Leptin'] #
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
#sig_limu=np.ones(33)*-2.3
#plt.plot(sig_limu, 'r')
#sig_limup=np.ones(33)*2.3
#plt.plot(sig_limup, 'r')


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

#
#
plt.show()
#

#load singular values distribution and plot
os.chdir(analysis_dir)
sing_p=np.load('singvals_p.npy')
print "p-value of salience of comp %i: %.4f" %(comp, sing_p[comp])

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

inertia=np.load('inertia.npy')
inertia_p=np.load('inertia_p.npy')
print "inertia of original SVD %d and its p-value: %d" %(inertia, inertia_p)

#explained variance
fig3=plt.figure()
plt.plot(sing_vals, 'ro')
exVar=sing_vals[comp]/inertia
print "explained variance of comp %i : %.3f" %(comp, exVar)


