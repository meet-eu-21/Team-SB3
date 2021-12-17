import os
import numpy as np
import scipy
import h5py
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

def pipeline(R,HiCfile,gene_density_file) :

    alpha=0.227
    HiCfilename= HiCfile # matrix containing the contact Hi-C
    
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
    
    filtered_mat, binsaved = HiCtoolbox.filteramat(binned_map)
    print('Input at the good resolution : ',np.shape(filtered_mat))
    
    
    
    
    ## Apply the SCN to the matrix
    
    binned_map_scn = HiCtoolbox.SCN(filtered_mat)
    
    print('Input at the good resolution : ',np.shape(binned_map_scn))
    
    
    ## Let us see a heatmap corresponding to the SCN binned_map :
    
    
    hm_scn = sns.heatmap(
        binned_map_scn, 
        cmap= "YlOrRd",
        square=True,
        norm = LogNorm()
    )
    
    ### Path of my own folder in cluster ##
    
    save_path = "/shared/home/apauron/"
    
    
    
    
    heatmap = HiCfile.replace(".RAWobserved","_heatmap.png")
    heatmap = heatmap.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    filename =os.path.basename(heatmap)
    directory = heatmap.replace(filename,"")
    Path(directory).mkdir(parents=True, exist_ok=True)
    fig = hm_scn.get_figure()
    fig.savefig(heatmap,dpi = 400)
    
    plt.close()
    
    ## Transform the matrix into the OE  matrix.
    
    oe_mat = binned_map_scn.copy()
    n = np.shape(oe_mat)[0]
    list_diagonals_mean = [np.mean(oe_mat.diagonal(i)) for i in range(-n+1,n)]
    
    for i in range(n):
        for j in range(n) :
            if list_diagonals_mean[j-i+n-1] == 0 :
                oe_mat[i,j] = 0
            else :       
                oe_mat[i,j] /= list_diagonals_mean[j-i+n-1]
    
        
    hm_oe = sns.heatmap(
        oe_mat, 
        cmap= "coolwarm",
        square=True,
        norm = LogNorm()
    )  
    
    oe_heatmap = HiCfile.replace(".RAWobserved","_oe_heatmap.png")
    oe_heatmap = oe_heatmap.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
  
    fig = hm_oe.get_figure()
    fig.savefig(oe_heatmap,dpi = 400)
    
    plt.close()
    
    ## Let us compute the correlation matrix
    
    corr = np.corrcoef(oe_mat)
    hm_corr = sns.heatmap(
        corr,
        vmin = -0.1, vmax = 0.1,
        cmap= "coolwarm",
        square=True,
    )  
    
    fig = hm_corr.get_figure()
    
    corr_heatmap = HiCfile.replace(".RAWobserved","_corr_heatmap.png")
    corr_heatmap = corr_heatmap.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
 
    fig.savefig(corr_heatmap,dpi = 400)
    
    plt.close()
    
    
    ## Get the svd decomposition
    
    eigenvalues, eigenvectors = np.linalg.eig(corr)
    
    eigenvectors = np.transpose(eigenvectors)
    
    s_vector = eigenvectors[0]
        
    
    
    f = h5py.File(gene_density_file, 'r')
    data_name = list(f.keys())[0]
    dset = f[data_name]
    data = dset[...]

    positive = True
    if np.sum( data[np.where(s_vector>0)]) < np.sum( data[np.where(s_vector<0)]) :
        positive = False

    
    z2 = np.zeros(len(s_vector))
    
    plt.plot(np.arange(len(s_vector)), s_vector)
    plt.xlabel("Index of the vector")
    plt.ylabel("First eigenvector value")
    plt.title("First eigenvector value")
    plt.tight_layout()
    
    if positive : 
        plt.fill_between(np.arange(len(s_vector)), s_vector, 0.0,
                         where=(s_vector >= z2),
                         color='red', label = "Compartment A")
        plt.fill_between(np.arange(len(s_vector)), s_vector, 0.0,
                         where=(s_vector < z2),
                         color='blue', label = "Compartment B")

        
        plt.legend()
    else :
        
        plt.fill_between(np.arange(len(s_vector)), s_vector, 0.0,
                         where=(s_vector < z2),
                         color='red',label = "Compartment A")
        plt.fill_between(np.arange(len(s_vector)), s_vector, 0.0,
                         where=(s_vector >= z2),
                         color='blue',label = "Compartment B")
        plt.legend()

    nameplot = HiCfile.replace(".RAWobserved","_compartments.png")
    nameplot = nameplot.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    
    plt.savefig(nameplot)
    plt.close()
    
    #3D
    print('3D')#Here : sparse int64
    contact_map=HiCtoolbox.SCN(filtered_mat.copy())
    contact_map[contact_map==0] = 0.000000001
    print("début_fast_floyd")
    contact_map=np.asarray(contact_map)**alpha #now we are not sparse at all
    dist_matrix = HiCtoolbox.fastFloyd(1/contact_map) #shortest path on the matrix
    dist_matrix=dist_matrix-np.diag(np.diag(dist_matrix))#remove the diagonal
    dist_matrix=(dist_matrix+np.transpose(dist_matrix))/2; #just to be sure that the matrix is symetric, not really usefull in theory
    print("fin_fast_floyd")
"""    
  
    XYZ,E=HiCtoolbox.sammon(dist_matrix, 3)#with the one from tom j pollard
    
    
    if positive :
        list_compartments = np.where(s_vector > 0,"A","B")
    else :
        list_compartments = np.where(s_vector > 0, "B","A")
    
    

            
    ### Convert to pdb file
    
    #point rescale
    hull=ConvexHull(XYZ)
    scale=100/hull.area**(1/3)
    XYZ=XYZ*scale
    
    pdbfilename  = HiCfile.replace(".RAWobserved","_3D.pdb")
    pdbfilename = pdbfilename.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    HiCtoolbox.writePDB(pdbfilename,XYZ,list_compartments)
        
                
"""

