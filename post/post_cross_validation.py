# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 11:12:03 2018

@author: fbeyer
"""

import numpy as np
import matplotlib.pyplot as plt

cv_dir="/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/CV_logtr_all_regressed_n9_wo_adipo_secondrun/"
corrs=np.load(cv_dir + "corr_values.npy")
print np.shape(corrs)


#plt.figure()
p_vals=np.zeros(shape=(np.shape(corrs)[0],np.shape(corrs)[1]))
for comp in np.arange(0,3):
    print "##################"
    print "component %i" %comp
    print "##################"
    for i in np.arange(0,20):
        #plt.subplot(5,4,i+1)
        #plt.hist(corrs[i,comp,1:])
        
        print "true correlation: %.3f" %corrs[i,comp,0]
        print sum(corrs[i,comp,1:]>corrs[i,comp,0])
        print len(corrs[i,comp,1:]>corrs[i,comp,0])
        #print corrs[i,1:]>corrs[i,0]
        p_val=sum(corrs[i,comp,1:]>corrs[i,comp,0])/1000.0
        print "p_value: %.3f" %(p_val)
        p_vals[i,comp]=p_val
        #plt.show()


ob_saliencies=np.load(cv_dir + "Ob_saliencies_splits.npy")
print np.shape(ob_saliencies)

Ob_titles=['BMI', 'WHR', 'HBA1C', 'CHOL', 'HDLC', 'CRPHS', 'IL6', 'Adiponectin','Leptin']
x=np.arange(0,np.shape(Ob_titles)[0])

#plotting related to obesity measures
#plotting salience ratios
#for comp in np.arange(0,2):
#    #plt.figure(num=comp)
#    for rep in np.arange(0,10):
#        plt.subplot(5,2,rep+1)
#        plt.xticks(x, Ob_titles, rotation=90)
#        plt.bar(x,ob_saliencies[rep,comp,:])
#        #
#        plt.tick_params(
#        axis='x',          # changes apply to the x-axis
#        which='both',      # both major and minor ticks are affected
#        bottom='off',      # ticks along the bottom edge are off
#        top='off',         # ticks along the top edge are off
#        labelbottom='on'
#        ) # labels along the bottom edge are off
#        #plt.rc('ytick', labelsize=16) 
#        #plt.rc('ytick', labelsize=16) 
#        plt.tick_params(axis='both', which='major', labelsize=18)
#        plt.ylabel('Ob saliencies/contributions to latent variable', fontsize=20)
#        plt.ylim((-0.5,0.5))

#plt.show()

print np.shape(p_vals)
##############################################################################
#mean correlation in 20 splits
mean_corr=np.mean(corrs[:,0,0])
sd_corr=np.std(corrs[:,0,0])

print "mean correlation across 20 splits for first component: %.4f +/- %.4f " %(mean_corr,sd_corr)
#p-values for individual 
print p_vals[:,0]
print len(p_vals[p_vals[:,0]>0.05,0])


##############################################################################
mean_corr=np.mean(corrs[:,1,0])
sd_corr=np.std(corrs[:,1,0])

print "mean correlation across 20 splits for second component: %.4f +/- %.4f " %(mean_corr,sd_corr)
#p-values for individual 
print p_vals[:,1]
print len(p_vals[p_vals[:,1]>0.05,0])