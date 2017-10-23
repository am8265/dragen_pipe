from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import numpy as np 
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

def plot_pca(training_genotypes,population_labels,output_plot_name):
    """ Plot Principal components

    training_genotypes : np.array (nsamples x 6) ; the scaled genotypes    
    population_labels : list ; the population labels
    output_plot_name : str ; name of the output plot file 
    """    
    df = pd.DataFrame(dict(x=training_genotypes[:,0],y=training_genotypes[:,1],z=training_genotypes[:,2],k=training_genotypes[:,3],
                           l=training_genotypes[:,4],m=training_genotypes[:,5],label=population_labels))
    groups = df.groupby('label')
    cols = dict(zip(np.unique(population_labels),['blue','green','red','cyan','yellow','orange','black']))
    print (cols)
    point_types = dict(zip(population_labels,['o','v','^','x','+','>']))
    ethnicities = {'1':'Caucasian','2':'MiddleEastern','3':'Hispanic','4':'East Asian', '5':'South Asian','6':'African','test_sample':'test'}
    ## Plot
    fig = plt.figure(figsize=(15,10))
    ax = fig.add_subplot(111,projection='3d')
    plot_legend = []
    i=0
    
    for name, group in groups:
        plot_legend.append(ethnicities[name])
        ax.scatter(group.x, group.y, group.z, depthshade=False, s=30,c=(0,0,0,0),marker='o',edgecolor=cols[name])
    ax.legend(plot_legend,loc='best')
    plt.title(' Training Genotypes projected on PC1 vs PC2 vs PC3')
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_zlabel('PC3')
    plt.savefig('training_pca.png')

    fig = plt.figure(figsize=(15,10))
    ax1 = fig.add_subplot(321)
    plot_legend = []
    i=0
    
    for name, group in groups:
        plot_legend.append(ethnicities[name])
        ax1.scatter(group.x, group.y, marker='o',color=cols[name])
    ax1.legend(plot_legend,loc='best')
    ax1.set_title(' PC1 vs PC2 ')
    ax1.set_xlabel('PC1')
    ax1.set_ylabel('PC2')

    ax2 = fig.add_subplot(322)
    plot_legend = []
    i=0
    
    for name, group in groups:
        plot_legend.append(ethnicities[name])
        ax2.scatter(group.y, group.z, marker='o',color=cols[name])
    ax2.legend(plot_legend,loc='best')
    ax2.set_title(' PC2 vs PC3 ')
    ax2.set_xlabel('PC2')
    ax2.set_ylabel('PC3')

    ax3 = fig.add_subplot(323)
    plot_legend = []
    i=0
    
    for name, group in groups:
        plot_legend.append(ethnicities[name])
        ax3.scatter(group.x, group.z, marker='o',color=cols[name])
    ax3.legend(plot_legend,loc='best')
    ax3.set_title(' PC1 vs PC3 ')
    ax3.set_xlabel('PC1')
    ax3.set_ylabel('PC3')

    ax4 = fig.add_subplot(324)
    plot_legend = []
    i=0
    
    for name, group in groups:
        plot_legend.append(ethnicities[name])
        ax4.scatter(group.z, group.k, marker='o',color=cols[name])
    ax4.legend(plot_legend,loc='best')
    ax4.set_title(' PC3 vs PC4 ')
    ax4.set_xlabel('PC3')
    ax4.set_ylabel('PC4')

    ax5 = fig.add_subplot(325)
    plot_legend = []
    i=0
    
    for name, group in groups:
        plot_legend.append(ethnicities[name])
        ax5.scatter(group.k, group.l, marker='o',color=cols[name])
    ax5.legend(plot_legend,loc='best')
    ax5.set_title(' PC4 vs PC5 ')
    ax5.set_xlabel('PC4')
    ax5.set_ylabel('PC5')

    ax6 = fig.add_subplot(326)
    plot_legend = []
    i=0
    
    for name, group in groups:
        plot_legend.append(ethnicities[name])
        ax6.scatter(group.l, group.m, marker='o',color=cols[name])
    ax6.legend(plot_legend,loc='best')
    ax6.set_title(' PC5 vs PC6 ')
    ax6.set_xlabel('PC5')
    ax6.set_ylabel('PC6')

    plt.tight_layout()
    plt.savefig(output_plot_name,dpi=600)    
    
      
    
