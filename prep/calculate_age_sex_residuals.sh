#!/bin/bash
##log-transform jacdet-files
fslmaths /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/merged_jacdet749.nii.gz -log /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/merged_jacdet749_log.nii.gz

#calculate age,sex residuals
fsl_glm -i /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/merged_jacdet749_log.nii.gz -d /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/age_sex_regressors_N752.txt -o /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/age_sex_results.nii.gz -m /data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/mask_for_Jacdet_onlyGM.nii.gz --out_res=/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/age_sex_residuals.nii.gz
