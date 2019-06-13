# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 12:35:43 2017

@author: fbeyer
"""
from nilearn import plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats, integrate
import seaborn as sb

analysis_dir="/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_medication/"

#with only one colour
#display = plotting.plot_glass_brain(analysis_dir + 'salience0_masked.nii.gz',\
#                                   plot_abs=False, cmap='bwr', colorbar=True)
#display.savefig('/home/raid1/fbeyer/Documents/Results/Metabolic_VBM/main_saliencies_brain_conv.eps')
#display.savefig('/home/raid1/fbeyer/Documents/Results/Metabolic_VBM/saliencies_brain_conv_bp.png')

display = plotting.plot_glass_brain(None)
display.add_contours(analysis_dir + 'salience0_masked.nii.gz',cmap='bwr')
#plt.show()
#display.savefig('/home/raid1/fbeyer/Documents/Results/Metabolic_VBM/main_saliencies_brain.eps')
display.savefig('/home/raid1/fbeyer/Documents/Results/Metabolic_VBM/saliencies_brain_medication_comp1.png')