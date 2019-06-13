
# coding: utf-8

# In[ ]:




# In[1]:

import scipy as sp
import numpy as np
import nibabel as nb
import os
import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt
#from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#####
def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)

#choose component to display

for comp in np.arange(0,2):
    
    #os.chdir('/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/prep/VBM_input/for_analysis/')
    #data=np.load('BMI_lab_measures_N748_logtr_all_reg.npy')
    Ob_titles_pretty=['BMI', 'WHR', 'HbA1c', 'cholesterol','HDL', 'IL6','CRP', 'adiponectin',  'leptin']
    #Ob_titles_pretty=['BMI', 'WHR', 'HbA1c', 'cholesterol','HDL', 'IL6','CRP', 'adiponectin',  'leptin', 'sys bp']#with blood pressure
    
    x=np.arange(0,len(Ob_titles_pretty))
    
    #load obesity measures bottstrap and plot
    print "loading files"
    analysis_dir='/data/pt_life/LIFE_bl/publications/2017_beyer_metabolic_pls_ohbm_poster/PLS_on_TBM_obesity_markers/PLS_analysis/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9/'    
    #analysis_dir="/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_medication/"
    #analysis_dir='/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_wo_adipo/'
    #analysis_dir='/data/pt_life/data_fbeyer/PLS_on_TBM_obesity_markers/PLS/results/sampleN748/N748_TIV_regression/perm_boot_logtr_all_regressed_n9_wo_crp/'
    os.chdir(analysis_dir)
    Ob_boo=np.load('Ob_salience_bootstrap.npy')
    Ob_boo_abs=abs(Ob_boo)
    
    Ob_sal=np.load('Ob_salience.npy')
    Ob_sal_abs=abs(Ob_sal)
    print(Ob_sal_abs[:,0])
    print(Ob_sal[:,0])
    # In[57]:
    #
    #sb.set(font_scale=2.8)
    #grid_kws = {"height_ratios": (0.8, .05), "hspace": 0.8}
    #fig, (ax, cbar_ax) = plt.subplots(2, figsize=(7, 7), dpi=300, gridspec_kw=grid_kws)
    ##fig,ax = plt.subplots(figsize=(14, 10), dpi=200)
    #
    ###import pandas as pd
    ##position=fig.add_axes([0.93,0.1,0.02,0.35])
    #
    #df=pd.DataFrame(data, columns=Ob_titles_pretty)
    #df.describe()
    #corr_df = df.corr(method='spearman')
    ##corr_df.
    #mask = np.zeros_like(corr_df)
    #mask[np.triu_indices_from(mask)] = True
    #nice_palette=sb.diverging_palette(220, 20, n=7, as_cmap=True)
    #nice_palette_2="RdBu_r"
    #im=sb.heatmap(corr_df, cmap=nice_palette_2, annot=True,fmt='.1f', annot_kws={"size":14}, mask=mask, linewidths=2.5,ax=ax, cbar_ax=cbar_ax, cbar_kws={"orientation":'horizontal', 'label': 'correlation coefficient'})#, "pad":0.3})
    ##b.tick_params(labelsize=5)
    ##vmax=1.0, vmin=-1.0 , 
    #ax.tick_params(labelsize=14)
    #plt.show()
    #
    #
    #
    ## In[15]:
    #
    
    #
    #fig0=plt.figure()
    #plt.xticks(x, Ob_titles_pretty, rotation=90)
    #plt.bar(x,Ob_sal[:,comp])
    ##
    #plt.tick_params(
    #    axis='x',          # changes apply to the x-axis
    #    which='both',      # both major and minor ticks are affected
    #    bottom='off',      # ticks along the bottom edge are off
    #    top='off',         # ticks along the top edge are off
    #    labelbottom='on'
    #    ) # labels along the bottom edge are off
    ##plt.rc('ytick', labelsize=16) 
    ##plt.rc('ytick', labelsize=16) 
    #plt.tick_params(axis='both', which='major', labelsize=18)
    #plt.ylabel('Ob saliencies/contributions to latent variable', fontsize=20)
    #plt.ylim((-0.5,0.5))
    #
    #plt.show()
    #
    #
    #
    #fig2=plt.figure()
    #plt.xticks(x, Ob_titles_pretty, rotation=90)
    #plt.bar(x,Ob_boo[:,comp])
    ##
    #plt.tick_params(
    #    axis='x',          # changes apply to the x-axis
    #    which='both',      # both major and minor ticks are affected
    #    bottom='off',      # ticks along the bottom edge are off
    #    top='off',         # ticks along the top edge are off
    #    labelbottom='on') # labels along the bottom edge are off
    ##plt.rc('ytick', labelsize=12) 
    ##plt.rc('ytick', labelsize=12) 
    #plt.tick_params(axis='both', which='major', labelsize=18)
    #plt.ylabel('Z-scored Ob effect from bootstrapping', fontsize=20)
    #plt.ylim((-2.6,2.6))
    #sig_limu=np.ones(len(Ob_titles_pretty))*-2.3
    #plt.plot(sig_limu, 'r')
    #sig_limup=np.ones(len(Ob_titles_pretty))*2.3
    #plt.plot(sig_limup, 'r')
    #plt.show()
    # In[21]:
    
    #plotting related to obesity measures
    #plotting salience ratios
    
    #shade the salience values with Z<2.3
    shading=abs(Ob_boo[:,comp])<2.3
    sh_num=np.arange(0,len(Ob_boo[:,comp]))[shading]
    
    print "start plotting"
    w=0.5
    offset = .5
    fig1,ax = plt.subplots(figsize=(14, 10), dpi=200)
    fig1.subplots_adjust(bottom=0.2)
    ls=25
    
    plt.xticks(x-0.5, Ob_titles_pretty, rotation=90,fontsize=ls)
    ax.grid(zorder=0)
    ax.bar(x-1,Ob_sal[:,comp],color='cornflowerblue', width=w,zorder=8)
    for i in sh_num:      
        ax.get_children()[i].set_alpha(0.6)
    
    #ax.get_children()[6].set_alpha(0.6)
    #ax.get_children()[1].set_linewidth(2)
    #ax.get_children()[2].set_linewidth(2)
    #ax.get_children()[7].set_linewidth(2)
    #ax.get_children()[8].set_linewidth(2)
    #ax.get_children()[0].set_color('mediumblue')
    #ax.get_children()[1].set_color('mediumblue')
    #ax.get_children()[2].set_color('mediumblue')
    #ax.get_children()[7].set_color('mediumblue')
    #ax.get_children()[8].set_color('mediumblue')
    
    ax.set_ylabel('Weights (a.u.)', fontsize=ls)
    plt.tick_params(axis='y', which='both',direction='out', length=4, width=1, labelsize=ls)
    plt.tick_params(axis='x', which='both', top='off',direction='out', length=4, width=1)
    ax.set_ylim(min(Ob_sal[:,comp])-0.1,max(Ob_sal[:,comp])+0.1)
    
    ax2=ax.twinx()
    ax2.grid(None)
    plt.xticks(x-0.5, Ob_titles_pretty, rotation=90,fontsize=ls)
    ax2.bar(x-w,Ob_boo[:,comp],color='indianred',width=w,zorder=7)
    for i in sh_num:      
        ax2.get_children()[i].set_alpha(0.6)
    #ax2.get_children()[6].set_alpha(0.6)
    #ax2.get_children()[1].set_color('firebrick')
    #ax2.get_children()[2].set_color('firebrick')
    #ax2.get_children()[7].set_color('firebrick')
    #ax2.get_children()[8].set_color('firebrick')
    ax2.set_ylabel('bootstrapped Z-scores', fontsize=ls)
    plt.tick_params(axis='y', which='both',direction='out', length=4, width=1, labelsize=ls)
    ax2.set_ylim(min(Ob_boo[:,comp])-1,max(Ob_boo[:,comp])+0.5)
    align_yaxis(ax, 0, ax2, 0)
    
    #plotting the Z line but based on saliency scale so that line is behind the bars.
    upp_f=2.3*ax.get_ylim()[1]/ax2.get_ylim()[1]
    low_f=-2.3*ax.get_ylim()[0]/ax2.get_ylim()[0]
    
    xx = np.arange(-1,len(Ob_titles_pretty))
    sig_limu=np.ones(len(Ob_titles_pretty)+1)*upp_f
    ax.plot(xx,sig_limu, 'firebrick',zorder=1)
    sig_limu=np.ones(len(Ob_titles_pretty)+1)*low_f
    ax.plot(xx,sig_limu, 'firebrick',zorder=1)
    
    #plt.show()
    
    #plt.savefig('/home/raid1/fbeyer/Documents/Results/Metabolic_VBM/main_saliencies_m9_pred_9.eps', format='eps', dpi=1000)
    filename='/home/raid1/fbeyer/Documents/Results/Metabolic_VBM/main_saliencies_m9_pred_10_bp_LV%i.jpg' %comp 
    plt.savefig(filename, format='jpg')
    
    
    
    




## In[257]:
#
#os.chdir('/afs/cbs.mpg.de/share/gr_agingandobesity/life_shared/Frauke/Genetics/old/FTO_paper_replication/PLS/')
###import pandas as pd
#df = pd.read_csv('320_older_subjects_usable_complete_BMI_WHR_lab_info_with_glucose_il6.csv')
#df = df.convert_objects(convert_numeric=True)
##temp = list(df.columns.values)
##Ob_titles = temp[1:]
#x=np.arange(-1,10)
#x2=np.arange(-0.5,10)
#Ob_titles = df.columns.values[1:]
#
#Ob_titles_pretty=['BMI', 'WHR', 'log(adiponectine)',  'log(leptin)', 'log(ghrelin)', 'HDL', 'LDL', 'cholesterol', 'Hba1c', 'log(CRP)', 'log(IL-6)']
#
##plotting related to obesity measures
##plotting salience ratios
#w=0.5
#offset = .5
#fig1,ax=plt.subplots()
#fig1.set_figheight(10)
#fig1.set_figwidth(10)
#plt.xticks(x2, Ob_titles_pretty, rotation=90,fontsize=20)
#
#rects1=ax.bar(x,Ob_sal[:,comp],color='cornflowerblue', width=w,zorder=5)
#ax.get_children()[0].set_color('mediumblue')
#ax.get_children()[2].set_color('mediumblue')
#ax.get_children()[3].set_color('mediumblue')
#ax.set_ylabel('Saliencies (a.u.)', fontsize=20)
#
#plt.tick_params(axis='y', which='both',direction='out', length=4, width=1, pad=-6, labelsize=16)
#plt.tick_params(axis='x', which='both', top='off',direction='out', length=4, width=1)
#
#ax2=ax.twinx()
#ax2.grid(None)
##plt.xticks(x, Ob_titles_pretty, rotation=90)
#ax2.bar(x+w,Ob_boo[:,comp],color='indianred',width=w,zorder=40)
#ax2.get_children()
#ax2.get_children()[0].set_color('firebrick')
#ax2.get_children()[3].set_color('firebrick')
#ax2.get_children()[2].set_color('firebrick')
#ax2.set_ylabel('Z-scores from bootstrapping', fontsize=20,labelpad=20)
#plt.tick_params(axis='y', which='both',direction='out', length=4, width=1, labelsize=16)
##ax2.xaxis.set_ticks(x2)
#ax2.set_xticks(x2)
##ax2.set_xticks([0])
#ax2.xaxis.grid(True, which='major')
##ax2.xaxis.set_ticklabels(Ob_titles_pretty,rotation=40,fontsize=20)
#
#xx = np.arange(-1,len(Ob_titles))
#sig_limu=np.ones(len(Ob_titles)+1)*-0.301
#ax.plot(xx,sig_limu, 'firebrick',zorder=1)
#sig_limu=np.ones(len(Ob_titles)+1)*0.301
#ax.plot(xx,sig_limu, 'firebrick',zorder=1)
#plt.subplots_adjust(bottom=0.25)
#fig1.savefig('sal_bootstrap_ohbm_abstract.png')#home/raid1/fbeyer/Documents/
#
#
## In[139]:
#
#fig1,ax=plt.subplots()
#
#ax.bar(x-1,Ob_sal[:,comp],color='cornflowerblue', width=w,zorder=5)
#ax2=ax.twinx()
#
#ax2.bar(x-0.5,Ob_boo[:,comp],color='cornflowerblue', width=w,zorder=4)
#xx = np.arange(-1,len(Ob_titles)+1)
#sig_limu=np.ones(len(Ob_titles)+2)*-0.301
#ax.plot(xx,sig_limu, 'firebrick',zorder=1)
#
##plt.xticks(x, labels, rotation='vertical')
#ax2.grid(None)
#
#
## In[65]:
#
#import matplotlib.pyplot as plt
#import matplotlib.patches as patches
#import matplotlib.colors as colors
#import math
#
#
#fig = plt.figure(figsize=(20,10), dpi=200)
#ax = fig.add_subplot(111)
#
#ratio = 1.0 / 3.0
#count = math.ceil(math.sqrt(len(colors.cnames)))
#x_count = count * ratio
#y_count = count / ratio
#x = 0
#y = 0
#w = 1 / x_count
#h = 1 / y_count
#
#for c in colors.cnames:
#    pos = (x / x_count, y / y_count)
#    ax.add_patch(patches.Rectangle(pos, w, h, color=c))
#    ax.annotate(c, xy=pos)
#    if y >= y_count-1:
#        x += 1
#        y = 0
#    else:
#        y += 1
#
#plt.show()
#
#
## In[199]:
#
#x2=np.arange(-0.5,10)
#print x2
#
#
## In[ ]:
#
#x2=np.arange(-0.5,10)
#
#import pandas as pd
#import seaborn as sns
#sns.set(font="monospace")
#
## Load the brain networks example dataset
#df = sns.load_dataset("brain_networks", header=[0, 1, 2], index_col=0)
#
## Select a subset of the networks
#used_networks = [1, 5, 6, 7, 8, 11, 12, 13, 16, 17]
#used_columns = (df.columns.get_level_values("network")
#                          .astype(int)
#                          .isin(used_networks))
#df = df.loc[:, used_columns]
#
## Create a custom palette to identify the networks
#network_pal = sns.cubehelix_palette(len(used_networks),
#                                    light=.9, dark=.1, reverse=True,
#                                    start=1, rot=-2)
#network_lut = dict(zip(map(str, used_networks), network_pal))
#
## Convert the palette to vectors that will be drawn on the side of the matrix
#networks = df.columns.get_level_values("network")
#network_colors = pd.Series(networks, index=df.columns).map(network_lut)
#
## Create a custom colormap for the heatmap values
#cmap = sns.diverging_palette(h_neg=210, h_pos=350, s=90, l=30, as_cmap=True)
#
## Draw the full plot
#sns.clustermap(df.corr(), row_colors=network_colors, linewidths=.5,
#               col_colors=network_colors, figsize=(13, 13), cmap=cmap)
#
