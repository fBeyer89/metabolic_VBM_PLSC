#!/bin/bash

fslmaths /data/pt_life/LIFE_bl/LIFE_2015_VBM_SPM12/LIFE_2015_VBM_SPM12/Templates/Template_6_GM.nii -thr 0.2 -bin /data/pt_life/LIFE_bl/LIFE_2015_VBM_SPM12/LIFE_2015_VBM_SPM12/Templates/Template_6_GM_mask.nii

for analysis in /data/pt_life/LIFE_bl/publications/2017_beyer_metabolic_pls_ohbm_poster/PLS_on_TBM_obesity_markers/PLS_analysis/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9 

#/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_medication /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_w_bp /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_wo_adipo /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_wo_crp /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_wo_IL6
do

cd $analysis

if [ -f *_salience0.nii.gz ];

then 

fslmaths *_bootstraps0.nii.gz -thr 2.3 -bin bootstraps0_masked.nii.gz
fslmaths *_salience0.nii.gz -mas bootstraps0_masked.nii.gz salience0_masked.nii.gz
#additional masking with GM template mask
fslmaths salience0_masked.nii.gz -mas /data/pt_life/LIFE_bl/LIFE_2015_VBM_SPM12/LIFE_2015_VBM_SPM12/Templates/Template_6_GM_mask.nii.gz salience0_masked_GMonly.nii.gz

fslmaths *_bootstraps1.nii.gz -thr 2.3 -bin bootstraps1_masked.nii.gz
fslmaths *_salience1.nii.gz -mas bootstraps1_masked.nii.gz salience1_masked.nii.gz
fslmaths salience1_masked.nii.gz -mas /data/pt_life/LIFE_bl/LIFE_2015_VBM_SPM12/LIFE_2015_VBM_SPM12/Templates/Template_6_GM_mask.nii.gz salience1_masked_GMonly.nii.gz

else

echo "need to do NIFTI conversion for $analysis"

fi

done
