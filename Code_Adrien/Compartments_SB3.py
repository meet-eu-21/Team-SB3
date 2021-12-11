import sys
import os
import numpy as np
import scipy
from matplotlib.colors import LogNorm, Normalize
from scipy.spatial import ConvexHull
from scipy import sparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.manifold import MDS #If you want the scikit learn mds
import HiCtoolbox
from pathlib import Path




    
#========================================= Simulation function ======================================

def pipeline(R,HiCfile,EpiGfile) :

    R=25000
    NbmaxEpi=15 #Epi states go from 0 to 15
    alpha=0.227
    selectedmark=1 #index of the selected mark
    HiCfilename= HiCfile # matrix containing the contact Hi-C
    EpiGfilename= EpiGfile #EpiGenetic file associating to each bin an epigenetic type
    
    #Build matrix
    A=np.loadtxt(HiCfilename) #Load the matrix
    A=np.int_(A) #To int
    print('Input data shape : ',np.shape(A))
    
    ## Bin the matrix
    
    A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)#build array at pb resolution
    A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
    binned_map=HiCtoolbox.bin2d(A,R,R) #!become csr sparse array
    LENTEST=np.shape(A)[0]
    print('Input at the good resolution : ',np.shape(binned_map))
        
    ##Filter the matrix
    
    filtered_mat = HiCtoolbox.filteramat(binned_map)[0]
    print('Input at the good resolution : ',np.shape(filtered_mat))
    
    
    
    
    ## Apply the SCN to the matrix
    
    binned_map_scn = HiCtoolbox.SCN(filtered_mat)
    
    print('Input at the good resolution : ',np.shape(binned_map_scn))
    
    
    ## Let us see a heatmap corresponding to the SCN binned_map :
    
    
    hm_scn = sns.heatmap(
        binned_map_scn, 
        vmin=0, vmax=0.1,
        cmap= "YlOrRd",
        square=True,
        norm = LogNorm()
    )
    
    ### Path of my own folder in cluster ##
    
    save_path = "/shared/home/apauron/"
    
    
    
    
    heatmap = HiCfile.replace(".RAWobserved","_heatmap.png")
    heatmap = HiCfile.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    Path(heatmap).mkdir(parents=True, exist_ok=True)
    fig = hm_scn.get_figure()
    fig.savefig(heatmap,dpi = 400)
    
    plt.close()
    
    ## Transform the matrix into the OE  matrix.
    
    oe_mat = binned_map_scn.copy()
    n = np.shape(oe_mat)[0]
    list_diagonals_mean = [np.mean(oe_mat.diagonal(i)) for i in range(-n+1,n)]
    
    for i in range(n):
        for j in range(n) :
            oe_mat[i,j] /= list_diagonals_mean[j-i+n-1]
    
    print(oe_mat)
        
    hm_oe = sns.heatmap(
        oe_mat, 
        vmin=0.01, vmax=100,
        cmap= "coolwarm",
        square=True,
        norm = LogNorm()
    )  
    
    oe_heatmap = HiCfile.replace(".RAWobserved","_oe_heatmap.png")
    oe_heatmap = HiCfile.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    Path(oe_heatmap).mkdir(parents=True, exist_ok=True)
    fig = hm_oe.get_figure()
    fig.savefig(oe_heatmap,dpi = 400)
    
    plt.close()
    
    ## Let us compute the correlation matrix
    
    corr = np.corrcoef(oe_mat)
    hm_corr = sns.heatmap(
        corr, 
        vmin=-1, vmax=1,
        cmap= "coolwarm",
        square=True,
    )  
    
    fig = hm_corr.get_figure()
    
    corr_heatmap = HiCfile.replace(".RAWobserved","_corr_heatmap.png")
    corr_heatmap = HiCfile.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    Path(corr_heatmap).mkdir(parents=True, exist_ok=True)
    fig.savefig(corr_heatmap,dpi = 400)
    
    plt.close()
    
    
    ## Get the svd decomposition
    
    eigenvalues, eigenvectors = np.linalg.eig(corr)
    
    eigenvectors = np.transpose(eigenvectors)
    
    nameplot = HiCfile.replace(".RAWobserved","_ev_plot.png")
    nameplot = HiCfile.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    Path(nameplot).mkdir(parents=True, exist_ok=True)
    plt.plot(np.arange(n),eigenvectors[0])
    plt.savefig(nameplot)
    


    