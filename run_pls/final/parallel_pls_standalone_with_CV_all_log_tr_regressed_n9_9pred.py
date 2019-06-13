# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 13:51:00 2015

@author: fbeyer
"""
#to activate virtuals environment
#cd ..
#source venvPLS/bin/activate
#!/usr/bin/env python2
import itertools
from multiprocessing import Pool, freeze_support, Array, Process, Queue
from functools import partial
import numpy as np
from sklearn.utils.extmath import randomized_svd
from scipy.stats.mstats_basic import zscore
import pyprind
from sklearn.utils import arpack
from sklearn.decomposition.truncated_svd import TruncatedSVD
from sklearn.preprocessing.data import scale
from scipy.stats.stats import percentileofscore
from functools import partial
import numpy as np
import scipy as sc
from scipy.stats.stats import pearsonr
from time import time
import os
import timeit

def fit_pls(X, Y, n_components, scale=True, algorithm="randomized"):
    #scaling

    print "calculating SVD"
    if scale:
        X_scaled = zscore(X, axis=0, ddof=1)
        Y_scaled = zscore(Y, axis=0, ddof=1)
        covariance = np.dot(Y_scaled.T, X_scaled)
    else:
        covariance = np.dot(Y.T, X)

    print np.shape(covariance)
    sum_var=covariance
    svd = TruncatedSVD(n_components, algorithm)
    #computes only the first n_components largest singular values
    #produces a low-rank approximation of covariance matrix
    Y_saliences, singular_values, X_saliences = svd._fit(covariance)
    X_saliences = X_saliences.T
    inertia = singular_values.sum()
    
    if scale:
        return X_saliences, Y_saliences, singular_values, inertia, X_scaled, Y_scaled, sum_var
    else:
        return X_saliences, Y_saliences, singular_values, inertia


def _procrustes_rotation(fixed, moving, moving_singular_values=None):
    N, _, P = np.linalg.svd(np.dot(fixed.T,moving))
    rotation_matrix = np.dot(N, P)
    rotated = np.dot(moving, rotation_matrix)
    
    if moving_singular_values != None:
        rotated_scaled = np.dot(np.dot(moving, np.diag(moving_singular_values)), rotation_matrix)
        rotated_singular_values = np.sqrt(np.square(rotated_scaled).sum(axis=0))
        return rotated, rotation_matrix, rotated_singular_values 
    else:
        return rotated, rotation_matrix 
   
    
def _permute_and_calc_singular_values_process(X, Y, a, b, n_components, algorithm, output, x): #perm_i
    """ basic version for parallel implementation using processes and output queue
    """
    
    #call random seed so not the same random number is used each time
    #pid = current_process()._identity[0]
    #randst = np.random.mtrand.RandomState(pid)
    np.random.seed( int( time() ) + x +50)
    
    #test how permutation works
    c=np.random.permutation(a)
    print a
    print c
    
    
    if len(X) < len(Y):
        #apply permutation to shorter list
        #print "randomization X<Y"
        X_perm = np.random.permutation(X)
        covariance_perm = np.dot(Y.T, X_perm)
    else:
        #print "other permutation"
        Y_perm = np.random.permutation(Y)
        covariance_perm = np.dot(Y_perm.T, X)

 
    
    svd = TruncatedSVD(n_components, algorithm=algorithm)
    
    #print covariance_perm
    Y_saliences_perm, singular_values_perm, X_saliences_perm =  svd._fit(covariance_perm) 

    output.put(singular_values_perm)

def _permute_and_calc_singular_values_pool(X, Y, X_saliences, Y_saliences, n_components,procrustes, algorithm, perm_i): 
    """ basic version for parallel implementation using pool
    """
    #call random seed so not the same random number is used in each process
    np.random.seed( int( time() ) + perm_i)
        
    
    if len(X) < len(Y):
        #apply permutation to shorter list
        #print "randomization X<Y"
        X_perm = np.random.permutation(X)
        covariance_perm = np.dot(Y.T, X_perm)
    else:
        #print "other permutation"
        Y_perm = np.random.permutation(Y)
        covariance_perm = np.dot(Y_perm.T, X)

 
    
    svd = TruncatedSVD(n_components, algorithm=algorithm)
    
    Y_saliences_perm, singular_values_perm, X_saliences_perm =  svd._fit(covariance_perm) 



    if procrustes:
    #It does not matter which side we use to calculate the rotated singular values
    #let's pick the smaller one for optimization
        if len(X_saliences_perm) > len(Y_saliences_perm):
            _, _, singular_values_perm = _procrustes_rotation(Y_saliences, Y_saliences_perm, singular_values_perm)
        else:
            X_saliences_perm = X_saliences_perm.T
            _, _, singular_values_perm = _procrustes_rotation(X_saliences, X_saliences_perm, singular_values_perm)

    
   
    return singular_values_perm
    
    

def permutation_test(X_scaled, Y_scaled, X_saliences, Y_saliences, singular_values,inertia, n_perm, verbose=True, algorithm="randomized"):
    n_components = X_saliences.shape[1]
    print "Starting permutations"
    #if verbose:
    #    my_perc = pyprind.ProgBar(n_perm, stream=1, title='running permutations', monitor=True)
    
    #import warnings
        #warnings.filterwarnings("ignore")
    
    
    
################################################################################    
    #create pool to run permutations in parallel    
    #procrustes=False
    iterable = np.arange(n_perm)
    P = Pool(processes = 20)
    func = partial(_permute_and_calc_singular_values_pool, X_scaled, Y_scaled,X_saliences, Y_saliences, n_components, True, algorithm)
    results=P.map(func, iterable)
    P.close()
    P.join()
    
    #cpu-count
    #multiprocessing.cpu_count()
    #if verbose:
    #     my_perc.update()   
   

    #if verbose:
    #    print my_perc
    #    print "calculating p values"
################################################################################
################################################################################    
    #use a list of processes and output queue
#    output = Queue()
#     

#    # Setup a list of processes that we want to run
#    processes = [Process(target=_permute_and_calc_singular_values, args=(X_scaled, Y_scaled, a, b, n_components, algorithm, output, x)) for x in range(4)]
#     
#    # Run processes
#    for p in processes:
#        p.start()
#     
#    # Exit the completed processes
#    for p in processes:
#        p.join()
#     
#    # Get process results from the output queue
#    results = [output.get() for p in processes]
#     
#    print(results)    
################################################################################
    print "end permutations"
    singular_values_samples=np.array(results).reshape(n_perm,n_components) #reshape results from list to np.array
    singvals_p_vals = np.zeros((n_components))#initialize saliences_p_vals
    for component_i in range(n_components):
        #percentileofscore compares rank to list of ranks (here singular value of component to bootstrapped
        #list of singular values
        singvals_p_vals[component_i] = (100.0-percentileofscore(singular_values_samples[:,component_i], singular_values[component_i]))/100.0
     
    #inertia describes explained variance 
    inertia_p_val = (100.0-percentileofscore(singular_values_samples.sum(axis=0), inertia))/100.0
    
    return  singvals_p_vals, inertia_p_val, singular_values_samples
    
###------------########----------------#############------------############# 
#BOOTSTRAPPING USING PROCESSES (outdated because of not finishing)
###------------########----------------#############------------#############   
def _boostrap_process(X, Y, X_saliences, Y_saliences, n_components, algorithm, output_X, output_Y, x):
    """ basic version for parallel implementation using processes and output queues
    """
    #call random seed so not the same random number is used in each process
    np.random.seed( int( time() ) + x) 
 
    #choose indices to resample randomly with replacement for a sample of same size
    sample_indices = np.random.choice(range(X.shape[0]), size=X.shape[0], replace=True)
    X_boot = X[sample_indices,:]
    Y_boot = Y[sample_indices,:]
    X_boot_scaled = scale(X_boot)
    Y_boot_scaled = scale(Y_boot)
    
    print np.shape(X_boot)
    
    covariance_boot = np.dot(Y_boot_scaled.T, X_boot_scaled)
    svd = TruncatedSVD(n_components, algorithm=algorithm)
    print "finished calculating SVD"
    Y_saliences_boot, _, X_saliences_boot = svd._fit(covariance_boot)
    X_saliences_boot = X_saliences_boot.T
    
    #It does not matter which side we use to calculate the rotated singular values
    #let's pick the smaller one for optimization
    print "rotating values"
    if len(X_saliences_boot) > len(Y_saliences_boot):
        #use procrustes_rotation on smaller dataset
        Y_bootstraps, rotation_matrix = _procrustes_rotation(Y_saliences, Y_saliences_boot)
        X_bootstraps = np.dot(X_saliences_boot, rotation_matrix)
    else:
        X_bootstraps, rotation_matrix = _procrustes_rotation(X_saliences, X_saliences_boot)
        Y_bootstraps = np.dot(Y_saliences_boot, rotation_matrix)  
       
      
    output_X.put(X_bootstraps)
    output_Y.put(Y_bootstraps)
    print "finished rotating" 

################################################################################   
def bootstrap_test(X, Y, X_saliences, Y_saliences, n_components,n_boot, verbose=True):
    print "starting bootstrap"
################################################################################      
    #use a list of processes and output queue
    output_X = Queue()
    output_Y = Queue() 

    #Setup a list of processes that we want to run
    processes = [Process(target=_boostrap_process, args=(X, Y, X_saliences, Y_saliences, n_components, 'randomized', output_X, output_Y, x)) for x in range(n_boot)]
     
    #Run processes
    for p in processes:
        p.daemon = True
        p.start()
        print "adding process"
    
    #p.close()
    #pool.join()     
    #Exit the completed processes
    #for p in processes:
    #    p.join()
    #    print "exiting processes"
     
    #Get process results from the output queue
    X_boo = [output_X.get() for p in processes]
    Y_boo = [output_Y.get() for p in processes]
       

    
    X_saliences_bootstraps=np.array(X_boo).reshape(n_boot,X_saliences.shape[0], X_saliences.shape[1])
    Y_saliences_bootstraps=np.array(Y_boo).reshape(n_boot,Y_saliences.shape[0], Y_saliences.shape[1])
    
    #sum across rows because the bootstraps are appended "under" each other
    X_saliences_bootstrap_ratios = X_saliences_bootstraps.mean(axis=0)/X_saliences_bootstraps.std(axis=0)
    Y_saliences_bootstrap_ratios = Y_saliences_bootstraps.mean(axis=0)/Y_saliences_bootstraps.std(axis=0)
    
    return X_saliences_bootstrap_ratios, Y_saliences_bootstrap_ratios      

###------------########----------------#############------------#############     
#BOOTSTRAPPING USING POOL    
###------------########----------------#############------------#############    
def _bootstrap_pool(X, Y, X_saliences, Y_saliences, n_components,procrustes, algorithm, boot_i): 
    """ basic version for parallel implementation of bootstrapping using pool
    """
    #call random seed so not the same random number is used in each process
    np.random.seed( int( time() ) + boot_i)
    #choose indices to resample randomly with replacement for a sample of same size
    sample_indices = np.random.choice(range(X.shape[0]), size=X.shape[0], replace=True)
    X_boot = X[sample_indices,:]
    Y_boot = Y[sample_indices,:]
    X_boot_scaled = scale(X_boot)
    Y_boot_scaled = scale(Y_boot)

    covariance_boot = np.dot(Y_boot_scaled.T, X_boot_scaled)
    svd = TruncatedSVD(n_components, algorithm=algorithm)
    Y_saliences_boot, _, X_saliences_boot = svd._fit(covariance_boot)
    X_saliences_boot = X_saliences_boot.T
    
    #It does not matter which side we use to calculate the rotated singular values
    #let's pick the smaller one for optimization
    if len(X_saliences_boot) > len(Y_saliences_boot):
        #use procrustes_rotation on smaller dataset
        Y_bootstraps, rotation_matrix = _procrustes_rotation(Y_saliences, Y_saliences_boot)
        X_bootstraps = np.dot(X_saliences_boot, rotation_matrix)
    else:
        X_bootstraps, rotation_matrix = _procrustes_rotation(X_saliences, X_saliences_boot)
        Y_bootstraps = np.dot(Y_saliences_boot, rotation_matrix)  
         
    
    #print np.shape(X_bootstraps)
    #print np.shape(Y_bootstraps)
   
    return X_bootstraps, Y_bootstraps

def bootstrap_pool(X, Y, X_saliences, Y_saliences, n_components,n_boot, verbose, write_dir):
    
    #bootstrap
    X_saliences_bootstraps = np.zeros(X_saliences.shape + (n_boot,))
    Y_saliences_bootstraps = np.zeros(Y_saliences.shape + (n_boot,))

    print "shape of X bootstraps: "
    print np.shape(X_saliences_bootstraps)    
    print "shape of Y bootstraps: "
    print np.shape(Y_saliences_bootstraps)  
    
    print "starting bootstraping"
    #create pool to run permutations in parallel    
    iterable = np.arange(n_boot)
    P = Pool(processes = 6)
    func = partial(_bootstrap_pool, X, Y, X_saliences, Y_saliences, n_components, True, 'randomized')
    res=P.map(func, iterable)
    P.close()
    P.join()


    X_saliences_bootstraps=[row[0] for row in res]
    Y_saliences_bootstraps=[row[1] for row in res]
    
    
    print "shape of X bootstraps: "
    print np.shape(X_saliences_bootstraps)    
    print "shape of Y bootstraps: "
    print np.shape(Y_saliences_bootstraps)  
    print "type of bootstraps:"
    print type(X_saliences_bootstraps)    
  
    #shape of bootstrapped X,Y is now (n_perm, vox_X/Y, n_comp) 

    #saving subresults        
    #saving results   
#    os.chdir(write_dir) 
#    np.save('Xsal_bootstrap_save_temp.npy', X_saliences_bootstraps)
#    np.save('Ysal_bootstrap_save_temp.npy', Y_saliences_bootstraps)
#    
#    
##    
##    #saving subresults     
##    #loading results
#    X_boo=np.load('Xsal_bootstrap_save_temp.npy')
#    Y_boo=np.load('Ysal_bootstrap_save_temp.npy')
    
    
#    print np.shape(X_boo)    
#    print np.shape(Y_boo)  
    
    #sum across first dimension because the bootstraps are along first dim
    X_saliences_bootstrap_asa = np.asarray(X_saliences_bootstraps)
    Y_saliences_bootstrap_asa = np.asarray(Y_saliences_bootstraps)
    X_saliences_bootstrap_ratios=X_saliences_bootstrap_asa.mean(axis=0)/X_saliences_bootstrap_asa.std(axis=0)
    Y_saliences_bootstrap_ratios=Y_saliences_bootstrap_asa.mean(axis=0)/Y_saliences_bootstrap_asa.std(axis=0)
#    
    print "done with summing across first dimensions"
        
    print np.shape(X_saliences_bootstrap_ratios)    
    print np.shape(Y_saliences_bootstrap_ratios)   
    
   
    
    return X_saliences_bootstrap_ratios, Y_saliences_bootstrap_ratios     
    
  
def run_pls(X, Y, n_components, n_perm, n_boot, write_dir):
    
    X_saliences, Y_saliences, singular_values, inertia, X_scaled, Y_scaled, sum_var = fit_pls(X, Y, 
                                                                                     n_components=n_components, 
     
                                                                  algorithm="randomized")
    #testing subresults
    

    os.chdir(write_dir)
    print "saving inputs"
    
    np.save('Xsal.npy', X_saliences)
    np.save('Ysal.npy', Y_saliences)
    np.save('sing.npy', singular_values)
    np.save('inert.npy', inertia)
    np.save('Xscaled.npy', X_scaled)
    np.save('Yscaled.npy', Y_scaled)
    np.save('sumvar.npy', sum_var)
    os.chdir("/home/raid1/fbeyer/Documents/Scripts/PLS/")
       
   
    print "\nstart permutations"
    singvals_p_vals, inertia_p_val, singular_values_samples = permutation_test(X_scaled, 
                                                                                    Y_scaled, 
                                                                                    X_saliences, 
                                                                                    Y_saliences, 
                                                                                    singular_values, 
                                                                                    inertia,
                                                                                    n_perm, 
                                                                                    verbose=True,
                                                                                    algorithm="randomized")

    print "\npermutations done"
    X_saliences_bootstrap_ratios, Y_saliences_bootstrap_ratios = bootstrap_pool(X, Y,X_saliences, Y_saliences, n_components,n_boot,True, write_dir) 
    
    print "\nbootstrapping done"
    
    

    
    return X_saliences, Y_saliences, singular_values, inertia, singvals_p_vals, inertia_p_val, singular_values_samples, X_saliences_bootstrap_ratios, Y_saliences_bootstrap_ratios   
  

if __name__ == "__main__":
    
    
    start_time = timeit.default_timer()
    write_dir="/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/CV_logtr_all_regressed_n9/"
    n_rep=20
    perm_n=1000
    n_comp=9
       
    analysis_dir='/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/N749/VBM_input/for_analysis/'
    brain_data=np.load(analysis_dir+'regression_N748/brain_data_mask0.3_N748.npy')
    Ob = np.load(analysis_dir + 'BMI_lab_measures_N748_logtr_all_reg.npy') 
   
    print np.shape(brain_data)[0]
    print np.floor(0.8*np.shape(brain_data)[0])

    corr_vals=np.zeros(shape=(n_rep,n_comp,perm_n+1))
    X_saliences_splits = np.zeros(shape=(n_rep,n_comp,np.shape(brain_data)[1]))
    Y_saliences_splits = np.zeros(shape=(n_rep,n_comp,np.shape(Ob)[1]))
    
    for i in np.arange(0,n_rep):
        #create 10 80% to 20% data splits (circa 598 vs 150)
        training=np.zeros(shape=np.shape(brain_data)[0],dtype=bool)
        
        training_ind=np.random.choice(range(np.shape(brain_data)[0]), np.floor(0.8*np.shape(brain_data)[0]), replace=False)
        training[training_ind]=True
        test=np.invert(training)
        
        #calculate SVD from training data set
        X_saliences, Y_saliences, singular_values, inertia, X_scaled, Y_scaled, sum_var = fit_pls(brain_data[training,:], Ob[training,:], n_comp, scale=True, algorithm="randomized")
        
        #multiply saliencies with test (leftout) dataset (first scale this one)
        brain_test_scaled=zscore(brain_data[test,:], axis=0, ddof=1)
        Ob_test_scaled = zscore(Ob[test,:], axis=0, ddof=1)
        
        for comp in np.arange(0,n_comp):
             X_saliences_splits[i,comp,:]=X_saliences[:,comp]    
             Y_saliences_splits[i,comp,:]=Y_saliences[:,comp] 
                   
            
             lv_brain=np.matrix(brain_test_scaled)*np.matrix(X_saliences)
             lv_ob=np.matrix(Ob_test_scaled)*np.matrix(Y_saliences)
            
             corr_vals[i,comp,0]=pearsonr(lv_ob[:,comp],lv_brain[:,comp])[0]
            
             #permute the smaller matrix of the training data and reproject
             for j in np.arange(1,perm_n+1):
                Ob_test_scaled_perm = np.random.permutation(Ob_test_scaled)
                lv_ob_perm=np.matrix(Ob_test_scaled_perm)*np.matrix(Y_saliences)
                
                corr_vals[i,comp,j]=pearsonr(lv_ob_perm[:,comp],lv_brain[:,comp])[0]
            
        

    os.chdir(write_dir)
    np.save('corr_values.npy', corr_vals)
    np.save('brain_saliencies_splits.npy', X_saliences_splits)
    np.save('Ob_saliencies_splits.npy', Y_saliences_splits)
    
    print "time elapsed for whole analysis of %i split repetitions and %i permutations in minutes" %(n_rep, perm_n)
    print((timeit.default_timer() - start_time)/60)
    

