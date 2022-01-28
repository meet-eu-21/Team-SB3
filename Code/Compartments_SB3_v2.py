import os
import numpy as np
import scipy
import h5py
from matplotlib.colors import LogNorm, Normalize
from scipy.spatial import ConvexHull
from scipy import sparse
import matplotlib.pyplot as plt
from hmmlearn import hmm
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import MDS #If you want the scikit learn mds
import HiCtoolbox
from pathlib import Path





    
#========================================= Simulation function ======================================

## Relu function is needed to plot eigenvector
def relu(x) :
    """
    Apply the relu function, i.e 0 if x is negative and x else.
    
    Keyword arguments :
    x -- the number or array to apply the function Relu
    
    Returns :
    
    the number or array modified
    """
    
    return np.maximum(x,np.zeros(len(x)))

def relumin(x) :
    
    """
    Apply the negative version of relu function, i.e 0 if x is positive and x else.
    
    Keyword arguments :
    x -- the number or array to apply the function relumin
    
    Returns :
    
    the number or array modified
    """
    
    return np.minimum(x,np.zeros(len(x)))

def SCN_Val(D, max_iter = 10):
    mat = np.float64(D.toarray())
    
    """
    A different version of SCN found in HiCtoolbox. Take a 2D array and apply the SCN transform on it.
    We do not make the matrix symetric at the end because it is rectangular
    
    
    Keyword arguments :
    D -- the 2D array to be transformed
    max_iter : the number of times you want to apply the SCN transform on your matrix
    
    Returns :
    
    mat -- transformed matrix
    """
    # Iteration over max_iter
    for i in range(max_iter):
        mat /= np.maximum(1, mat.sum(axis = 0))
        mat /= np.maximum(1, mat.sum(axis = 1)[:,None])
    return mat

def filteramat_Val(Hicmat,Filterextremum=True,factor=1.5):
    """
    A different version of filteramat function in HiCtoolbox, adapted for interchromosomal contact (rectangle matrix)
    Take a 2D array and apply a filter on it.
    It filters lines and columns with a 0 sum, and line and columns with a sum over mean + std*factor
    
    
    Keyword arguments :
    Hicmat -- the 2D array to be filtered
    Filterextremum -- boolean, by defaut on True to check if you want to filter extremums.
    factor -- double , parameter which defines bounds of the filter
    
    Returns :
    
    Hicmatreduce -- filtered matrix
    segmenter0 -- list of saved bin for first chromosome (in x axis)
    segmenter1 -- list of saved bins for second chromosome (y axis)
    
    """
        
    Hicmatreduce=Hicmat
    
    ##Filter lines and columns whose sum is 0
    
    A_0 = []
    for i in range(np.shape(Hicmat)[0]):
        if np.sum(Hicmat[i,:]) > 0:
            A_0.append(i)
    
    # Filter the lines contained in A_0
    Hicmatreduce=Hicmatreduce[A_0,:] 
    
    A_1 = []
    for j in range(np.shape(Hicmat)[1]):
        if np.sum(Hicmat[:,j]) > 0:
            A_1.append(j)
    
    # Filter the columns contained in A_0

    Hicmatreduce=Hicmatreduce[:,A_1]

    ##Filter lines and columns whose sum is 0 again, to be sure.

    sumHicmat0= Hicmatreduce.sum(axis = 0)
    segmenter0=sumHicmat0 > 0
    A_0 =np.where(segmenter0)
    
    # Filter the lines contained in A_0

    Hicmatreduce=Hicmatreduce[:,A_0[1]]

    sumHicmat1 = Hicmatreduce.sum(axis = 1)
    segmenter1 = sumHicmat1 > 0
    A_1 = np.where(segmenter1)
    
    # Filter the columns contained in A_0

    Hicmatreduce=Hicmatreduce[A_1[0],:]

    # If you want to filter :
    
    if Filterextremum:
        #second step : filter lower bin for first chromosome
        sumHicmat_0=np.sum(Hicmatreduce,axis=0)
        msum0=np.mean(sumHicmat_0)
        mstd0=np.std(sumHicmat_0)
        #get mean and standard deviation of sums

        mini0 = msum0-mstd0*factor
        maxi0 = msum0+mstd0*factor
        #Make the bolean condition
        newcond0=mini0 < sumHicmat_0
        newcond20=sumHicmat_0 < maxi0
        newcond0=np.logical_and(newcond0,newcond20)
        B0=np.where(newcond0)

        #third step : filter lower bin for second chromosome
        sumHicmat_1=np.sum(Hicmatreduce,axis = 1)
        #get mean and standard deviation of sums
        msum1=np.mean(sumHicmat_1)
        mstd1=np.std(sumHicmat_1)
        mini1 = msum1-mstd1*factor
        maxi1 = msum1+mstd1*factor
        #Make the bolean condition
        newcond1=mini1 < sumHicmat_1
        newcond21=sumHicmat_1 < maxi1
        newcond1=np.logical_and(newcond1,newcond21)
        B1=np.where(newcond1)
        
        #Filter
        Hicmatreduce=Hicmatreduce[:,B0[1]]
        Hicmatreduce=Hicmatreduce[B1[0],:]
        
        segmenter0=A_0[1][B0[1]] #Create the binsaved index for first chromosome
        segmenter1=A_1[0][B1[0]] #Create the binsaved index for second chromosome
    return Hicmatreduce,segmenter0,segmenter1

###  Intrachromosomal Pipeline  --------------------------------------------------------------------
def pipeline_intra(R,HiCfile,gene_density_file) :

    """
    A pipeline for generating intrachromosomal compartments. Mainly inspired by Leopold Carron    
    
    Keyword arguments :
    R -- integer, the resolution of the file
    HiCfile -- path, the relative path of the HiCfile you want to get the compartments from
    density_file -- path , density file corresponding to our HiCfile chromosome
    
    Returns :
    
    None
    
    Save figures in Results_local folder created at the working directory
    
    """
    
    alpha=0.227 #Useful for pdb generation
    HiCfilename= HiCfile # matrix containing the contact Hi-C
    
    #Build matrix
    A=np.loadtxt(HiCfilename) #Load the matrix
    A=np.int_(A) #To int
    print('Input data shape : ',np.shape(A))
    
    ## Bin the matrix
    
    A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)#build array at pb resolution
    A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
    binned_map=HiCtoolbox.bin2d(A,R,R) #!become csr sparse array
    
    print('Resolution after binning : ',np.shape(binned_map))
        
    ##Filter the matrix
    
    filtered_mat, binsaved = HiCtoolbox.filteramat(binned_map)
    print('Resolution after filtering : ',np.shape(filtered_mat))
    
    
    
    
    ## Apply the SCN to the matrix
    
    binned_map_scn = HiCtoolbox.SCN(filtered_mat)
    
    print('Resolution after SCN: ',np.shape(binned_map_scn))
    
    
    ## Let us see a heatmap corresponding to the SCN binned_map :
    
    
    hm_scn = sns.heatmap(
        binned_map_scn, 
        cmap= "YlOrRd",
        square=True,
        norm = LogNorm()
    )
    
    ### Path of my own results folder ##
    
    save_path = os.path.join(os.getcwd(),"Results_local")
    
    
    
    ## Create the name of the heatmap file 
    heatmap = HiCfile.replace(".RAWobserved","_heatmap.png")
    heatmap = os.path.join(save_path,heatmap)
    
    ## Save the file
    filename =os.path.basename(heatmap)
    directory = heatmap.replace(filename,"")
    Path(directory).mkdir(parents=True, exist_ok=True)
    fig = hm_scn.get_figure()
    fig.savefig(heatmap,dpi = 400)
    
    plt.close()
    
    ## Transform the matrix into the OE  matrix.
    
    oe_mat = binned_map_scn.copy() #A copy of SCN generated file
    n = np.shape(oe_mat)[0] #Symetrical matrix
    list_diagonals_mean = [np.mean(oe_mat.diagonal(i)) for i in range(-n+1,n)]
    
    for i in range(n):
        for j in range(n) :
            if list_diagonals_mean[j-i+n-1] == 0 :
                oe_mat[i,j] = 0
            else :       
                oe_mat[i,j] /= list_diagonals_mean[j-i+n-1]
    
    ## Generate a heatmap of the oe matrix
    
    hm_oe = sns.heatmap(
        oe_mat, 
        cmap= "coolwarm",
        square=True,
        norm = LogNorm()
    )  
    
    ## Get the filename of the file
    oe_heatmap = HiCfile.replace(".RAWobserved","_oe_heatmap.png")
    oe_heatmap = os.path.join(save_path,oe_heatmap)
  
    fig = hm_oe.get_figure()
    fig.savefig(oe_heatmap,dpi = 400)
    
    plt.close()
    
    ## Let us compute the correlation matrix
    
    corr = np.corrcoef(oe_mat) #Correlation matrix in numpy format
    
    ## Correlation matrix heatmap style
    hm_corr = sns.heatmap(
        corr,
        vmin = -0.1, vmax = 0.1,
        cmap= "coolwarm",
        square=True,
    )  
        
    fig = hm_corr.get_figure()
    
    ## Get the filename and save it
    
    corr_heatmap = HiCfile.replace(".RAWobserved","_corr_heatmap.png")
    corr_heatmap = os.path.join(save_path,corr_heatmap)
 
    fig.savefig(corr_heatmap,dpi = 400)
    
    plt.close()
    
    
    ## Do the PCA analysis
    
        
    pca = PCA(n_components=2)
    pca.fit(corr)
    

    ## Get the first eigenvector
    
    s_vector = pca.components_[0]
    
    ## Turn it into a list
    list_vector = list(s_vector)
    
    ## Get the name of eigenvector files
    vfile = HiCfile.replace(".RAWobserved","_vp.txt")
    vfile = os.path.join(save_path,vfile)

    ## Write the corresponding eigenvector in a file 

    list_filtered = [i for i in range(1,np.shape(binned_map)[0]) if i not in binsaved ]
    
    for x in list_filtered :
        list_vector.insert(x-1,-1.0)

    
    with open(vfile,'w') as f :
        for x in list_vector :
            f.write(str(x) + "\n")
    
    ## Open the gene density file
         
    f = h5py.File(gene_density_file, 'r')
    data_name = list(f.keys())[0]
    dset = f[data_name]
    data = dset[...]

    ## Check which compartments are euchromatin or heterochromatin based on their density.
    
    positive = True # A = dense and B = not dense
    if np.sum( data[np.where(s_vector>0)]) < np.sum( data[np.where(s_vector<0)]) :
        positive = False # A = not dense and B = dense

    compfile = HiCfile.replace(".RAWobserved","_comp.txt")
    compfile = os.path.join(save_path,compfile)

    ## Write the comps
    ## Do the HMM analysis
    
    remodel = hmm.GaussianHMM(n_components=2, covariance_type="diag", n_iter=1000)

    remodel.fit(corr)

    Z2 = remodel.predict(corr)
    list_compartments = []
    
    ## If positive scenario, write 1.0 for A compartment and 0.0 for B
    
    if positive :
        for i in range(len(Z2)) :
            if Z2[i] == 0 :
                list_compartments.append(0.0)
            if Z2[i] == 1 :
                list_compartments.append(1.0)
    
    ## Do the opposite for negative scenario, write 1.0 for B compartment and 0.0 for A

    else :
        for i in range(len(Z2)) :
            if Z2[i] == 0 :
                list_compartments.append(1.0)
            if Z2[i] == 1 :
                list_compartments.append(0.0)
  
    ## Insert filtered bins as -1.0
    
    for x in list_filtered :
        list_compartments.insert(x-1,-1.0)
    
    ## Write the file

    with open(compfile,'w') as f :
        for x in list_compartments :
            f.write(str(x) + "\n")
    
    ## Plot the first eigenvector and save it
    
    plt.plot(np.arange(len(s_vector)), s_vector)
    plt.xlabel("Index of the vector")
    plt.ylabel("First eigenvector value")
    plt.title("First eigenvector value")
    plt.tight_layout()
    
    if positive :
        plt.plot(np.arange(len(s_vector)), relu(s_vector),'r')

        plt.fill_between(np.arange(len(s_vector)), relu(s_vector), 0,
                         color='red', label = "Compartment A")
        plt.plot(np.arange(len(s_vector)), relumin(s_vector),'b')

        plt.fill_between(np.arange(len(s_vector)), relumin(s_vector), 0,
                         color='blue', label = "Compartment B")

        
        plt.legend()
    else :
        plt.plot(np.arange(len(s_vector)), relumin(s_vector),'r')
        plt.fill_between(np.arange(len(s_vector)), relumin(s_vector), 0,
                         color='red',label = "Compartment A")
        plt.plot(np.arange(len(s_vector)), relu(s_vector),'b')
        plt.fill_between(np.arange(len(s_vector)), relu(s_vector), 0,
                         color='blue',label = "Compartment B")
        plt.legend()

    nameplot = HiCfile.replace(".RAWobserved","_compartments.png")
    nameplot = os.path.join(save_path,nameplot)
    
    plt.savefig(nameplot)
    plt.close()
    

    #Creation of the 3D structure of the chromosome.
    
    print('3D')#Here : sparse int64
    contact_map=HiCtoolbox.SCN(filtered_mat.copy()) ## Re-apply SCN
    contact_map[contact_map==0] = 0.000000001 ## To avoid infinite error during fast floyd
    contact_map=np.asarray(contact_map)**alpha #now we are not sparse at all
    
    print("begin fast floyd")
    dist_matrix = scipy.sparse.csgraph.floyd_warshall(1/contact_map,return_predecessors = False) #shortest path on the matrix
    print("fast floyd ended")
    dist_matrix=dist_matrix-np.diag(np.diag(dist_matrix))#remove the diagonal
    dist_matrix=(dist_matrix+np.transpose(dist_matrix))/2; #just to be sure that the matrix is symetric, not really usefull in theory
    

    
    #MDS to generate 3D coordinates from the distance matrix
    embedding = MDS(n_components=3)#LOAD the MDS #With scikit-learn mds
    XYZ = embedding.fit_transform(dist_matrix) #Make the transform
    print("MDS transform done")

    ## For the pdb file, put A and B as alternate locations for each Carbon    
    
    if positive :
        list_compartments = np.where(s_vector > 0,"A","B")
    else :
        list_compartments = np.where(s_vector > 0, "B","A")
    

            
    ### Convert to pdb file
    
    #point rescale
    hull=ConvexHull(XYZ)
    scale=100/hull.area**(1/3)
    print("ConvexHull done")
    XYZ=XYZ*scale
    
    ## Save the pdb file
    
    pdbfilename  = HiCfile.replace(".RAWobserved","_3D.pdb")
    pdbfilename = os.path.join(save_path,pdbfilename)
    HiCtoolbox.writePDB(pdbfilename,XYZ,list_compartments)
    print("PDB written")


###  Interchromosomal pipeline --------------------------------------------------------------------
    
def pipeline_inter(R,HiCfile,gene_density_file_1,gene_density_file_2) :
    
    """
    A pipeline for generating interchromosomal compartments. Mainly inspired by intra version    
    
    Keyword arguments :
    R -- integer, the resolution of the file
    HiCfile -- path, the relative path of the HiCfile you want to get the compartments from
    density_file_1 -- path , density file corresponding to our HiCfile first chromosome (in order of apparition)
    density_file_2 -- path , density file corresponding to our HiCfile second chromosome (in order of apparition)

    
    Returns :
    
    None
    
    Save figures in Results_local folder created at the working directory

    
    """

    HiCfilename= HiCfile # matrix containing the contact Hi-C
    GDfilename_1 = gene_density_file_1 ## A short for the first gene density file
    GDfilename_2 = gene_density_file_2 ## A short for the second gene density file
    
    ## Get the ID of each chromosome :
    
    num_chr_1 = GDfilename_1.split('chr')[1].split('.')[0]
    num_chr_2 = GDfilename_2.split('chr')[1].split('.')[0]
    
    #Build matrix : 
    
    A=np.loadtxt(HiCfilename) #Load the matrix
    A=np.int_(A) #To int
    print('Input data shape A : ',np.shape(A)) # Nombre de lignes * nombre de colonnes du fichier HiCfilename

    ## Bin the matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION from VAL : No concatenation of A and its transposition because it forces a square matrix :
   
    
    # --------------------------------------------------------------------------
    A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
    print('Input data shape A post coo_matrix : ',np.shape(A))
    binned_map=HiCtoolbox.bin2d(A,R,R) #!become csr sparse array

    
    ## Check the resolution of the binned map
    
    print('Input at the good resolution (binned map) : ',np.shape(binned_map))
        
    ##Filter the matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : A new filtration method
    
    
    filtered_mat, binsaved_Y, binsaved_X = filteramat_Val(binned_map)

    # --------------------------------------------------------------------------
    
    print('Input at the good resolution (filtered map) : ',np.shape(filtered_mat))

    
    ## Apply the SCN to the matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : SCN handling rectangle matrices
    
    binned_map_scn = SCN_Val(filtered_mat)
    
    # --------------------------------------------------------------------------
        
    print('Input at the good resolution after SCN : ',np.shape(binned_map_scn))
    
    ## Let us see a heatmap corresponding to the SCN binned_map :
    
    
    hm_scn = sns.heatmap(
        binned_map_scn, 
        cmap= "YlOrRd",
        square=False,
        norm = LogNorm()
    )
    
    ### Path of our results folder ##
    
    save_path = os.path.joint(os.getcwd(),"Results_local")
        
    ## Prepare file name of SCN matrix and save it
    heatmap = HiCfile.replace(".RAWobserved","_heatmap.png")
    heatmap = os.path.join(save_path,heatmap)
    
    ## Name without full path
    filename =os.path.basename(heatmap)
    directory = heatmap.replace(filename,"")
    
    ## Create directories if they do not exist
    Path(directory).mkdir(parents=True, exist_ok=True)
    fig = hm_scn.get_figure()
    fig.savefig(heatmap,dpi = 400)
    
    plt.close()
    
    ## Transform the matrix into the OE  matrix.
    
    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : Creation of rectangle OE matrix
    
    oe_mat = binned_map_scn.copy()
    
    #Rectangle matrix
    n,m = np.shape(oe_mat)
    
    list_diagonals_mean = [np.mean(oe_mat.diagonal(i)) for i in range(-n+1,m)]
    
    ## Iteration on lines and columns
    for i in range(n):
        for j in range(m) :
            if list_diagonals_mean[j-i+n-1] == 0 :
                oe_mat[i,j] = 0
            else :       
                oe_mat[i,j] /= list_diagonals_mean[j-i+n-1]
    # --------------------------------------------------------------------------

    ## Let us see a heatmap corresponding to the OE binned_map :

    hm_oe = sns.heatmap(
        oe_mat, 
        cmap= "coolwarm",
        square=False,
        norm = LogNorm()
    )  
    
    ## Prepare the file and save it 
    
    oe_heatmap = HiCfile.replace(".RAWobserved","_oe_heatmap.png")
    oe_heatmap = os.path.join(save_path,oe_heatmap)
  
    fig = hm_oe.get_figure()
    fig.savefig(oe_heatmap,dpi = 400)
    
    plt.close()
    
    ## Let us compute the correlation matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : 
    # Two correlation matrices are computed : one for chromosome 1 and the other for chromosome 2
    

    
    corr_1 = np.corrcoef(oe_mat)
    ## Apply correlation on transposed matrix
    corr_2 = np.corrcoef(oe_mat,rowvar = False)
    
    ## Save first correlation matrix as a figure
    hm_corr_1 = sns.heatmap(corr_1, vmin = -0.1, vmax = 0.1, cmap= "coolwarm",square=False)  
    fig = hm_corr_1.get_figure()
    
    corr_heatmap_1 = HiCfile.replace(".RAWobserved","_corr_heatmap_"+str(num_chr_1)+".csv")
    corr_heatmap_1 = os.path.join(save_path,corr_heatmap_1)
    
    fig.savefig(corr_heatmap_1,dpi = 400)
    plt.close()
    
    ## Save second one as a figure 
    hm_corr_2 = sns.heatmap(corr_2, vmin = -0.1, vmax = 0.1, cmap= "coolwarm",square=False)  
    fig = hm_corr_2.get_figure()
    
    corr_heatmap_2 = HiCfile.replace(".RAWobserved","_corr_heatmap_"+str(num_chr_2)+".csv")
    corr_heatmap_1 = os.path.join(save_path,corr_heatmap_2)
    
    fig.savefig(corr_heatmap_2,dpi = 400)
    plt.close()
    
    # --------------------------------------------------------------------------

    print("Computation of Correlation Matrices :",np.shape(corr_1),'/',np.shape(corr_2))
    ## Nan error can happen, meaning a correlation of 0.
    if np.isnan(corr_1).any() :
        corr_1[np.isnan(corr_1)] = 0

    if np.isnan(corr_2).any() :
        corr_2[np.isnan(corr_2)] = 0
        
    
    ## Get the svd decomposition

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : Instead of one, searching for first eigenvector of the two correlation matrices
    
    # First matrix 
    
    eigenvalues_1, eigenvectors_1 = np.linalg.eig(corr_1)
    eigenvectors_1 = np.transpose(eigenvectors_1)
    s_vector_1 = eigenvectors_1[0]

    # Second matrix 
    
    eigenvalues_2, eigenvectors_2 = np.linalg.eig(corr_2)
    eigenvectors_2 = np.transpose(eigenvectors_2)
    s_vector_2 = eigenvectors_2[0]
    
    # --------------------------------------------------------------------------
    
    print("Computation of eigenvectors :",len(s_vector_1),'/',len(s_vector_2))

    ## Gene Density data treatment
    
    # For first chromosome
    
    f1 = h5py.File(GDfilename_1, 'r')
    data_name_1 = list(f1.keys())[0]
    dset_1 = f1[data_name_1]
    data_1 = dset_1[...]
    data_1 = data_1[binsaved_X]

    plt.plot(np.arange(len(data_1)),data_1)
    plt.savefig(save_path + "/gene_density_"+str(num_chr_1)+".png")
    plt.close()
    # --------------

    # For second chromosome

    f2 = h5py.File(GDfilename_2, 'r')
    data_name_2 = list(f2.keys())[0]
    dset_2 = f2[data_name_2]
    data_2 = dset_2[...]
    data_2 = data_2[binsaved_Y]

    plt.plot(np.arange(len(data_2)),data_2)
    plt.savefig(save_path + "/gene_density_"+str(num_chr_2)+".png")
    plt.close()

    print("Study the gene densities")

    ## Find the compartments

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : Study of inversion of compartments A and B 
    
    ## Positive boolean is based on density of chromosome and affect A and B compartments to euchromatin and heterochromatin
    positive_1 = True
    if np.sum( data_1[np.where(s_vector_1>0)]) < np.sum( data_1[np.where(s_vector_1<0)]) :
        positive_1 = False

    print("First Gene Density positivity :",positive_1)
    
    #Plot first chromosome first eigenvector
    
    plt.plot(np.arange(len(s_vector_1)), s_vector_1)
    plt.xlabel("Index of the vector")
    plt.ylabel("First eigenvector value")
    plt.title("First eigenvector value for chromosome "+str(num_chr_1))
    plt.tight_layout()
    
    if positive_1 :
        plt.plot(np.arange(len(s_vector_1)), relu(s_vector_1),'r')
        
        # Fill under the curve to get compartments


        plt.fill_between(np.arange(len(s_vector_1)), relu(s_vector_1), 0,
                         color='red', label = "Compartment A")
        plt.plot(np.arange(len(s_vector_1)), relumin(s_vector_1),'b')
        # Fill under the curve to get compartments

        plt.fill_between(np.arange(len(s_vector_1)), relumin(s_vector_1), 0,
                         color='blue', label = "Compartment B")

        
        plt.legend()
        
    else :
        plt.plot(np.arange(len(s_vector_1)), relumin(s_vector_1),'r')
        ## Fill under the curve to get compartments
        plt.fill_between(np.arange(len(s_vector_1)), relumin(s_vector_1), 0,
                         color='red',label = "Compartment A")
        plt.plot(np.arange(len(s_vector_1)), relu(s_vector_1),'b')
        
        ## Fill under the curve to get compartments

        plt.fill_between(np.arange(len(s_vector_1)), relu(s_vector_1), 0,
                         color='blue',label = "Compartment B")
        plt.legend()

    ## Save the plot
    nameplot = HiCfile.replace(".RAWobserved","_compartments_"+str(num_chr_1)+".png")
    nameplot = os.path.join(save_path,nameplot)
    
    plt.savefig(nameplot)
    plt.close()

    # --------------

    ## Do the same for chromosome 2 
    
    positive_2 = True
    if np.sum( data_2[np.where(s_vector_2>0)]) < np.sum( data_2[np.where(s_vector_2<0)]) :
        positive_2 = False

    print("Second Gene Density positivity :",positive_2)
        
    plt.plot(np.arange(len(s_vector_2)), s_vector_2)
    plt.xlabel("Index of the vector")
    plt.ylabel("First eigenvector value")
    plt.title("First eigenvector value for chromosome "+str(num_chr_2))
    plt.tight_layout()
    
    if positive_2 :
        plt.plot(np.arange(len(s_vector_2)), relu(s_vector_2),'r')

        plt.fill_between(np.arange(len(s_vector_2)), relu(s_vector_2), 0,
                         color='red', label = "Compartment A")
        plt.plot(np.arange(len(s_vector_2)), relumin(s_vector_2),'b')

        plt.fill_between(np.arange(len(s_vector_2)), relumin(s_vector_2), 0,
                         color='blue', label = "Compartment B")

        
        plt.legend()
    else :
        plt.plot(np.arange(len(s_vector_2)), relumin(s_vector_2),'r')
        plt.fill_between(np.arange(len(s_vector_2)), relumin(s_vector_2), 0,
                         color='red',label = "Compartment A")
        plt.plot(np.arange(len(s_vector_2)), relu(s_vector_2),'b')
        plt.fill_between(np.arange(len(s_vector_2)), relu(s_vector_2), 0,
                         color='blue',label = "Compartment B")
        plt.legend()

    nameplot = HiCfile.replace(".RAWobserved","_compartments_"+str(num_chr_2)+".png")
    nameplot = os.path.join(save_path,nameplot)
    
    plt.savefig(nameplot)
    plt.close()
    # --------------------------------------------------------------------------

    print("Find the compartments with HMM")

    ## Get the filtered bins for chr 1 and 2 
    
    list_filtered1 = [i for i in range(1,np.shape(binned_map)[0]) if i not in binsaved_X ]
    list_filtered2 = [i for i in range(1,np.shape(binned_map)[1]) if i not in binsaved_Y ]
    
    
    ## Do the HMM analysis for chr 1
    
    ## Get the filename
    compfile1 = HiCfile.replace(".RAWobserved","_comp_" + str(num_chr_1)+".txt")
    compfile1 = os.path.join(save_path,compfile1)
    
    ## Apply Gaussian HMM
    remodel = hmm.GaussianHMM(n_components=2, covariance_type="diag", n_iter=1000)

    remodel.fit(corr_1)

    Z21 = remodel.predict(corr_1)
    list_compartments1 = []
    
    ## Create list of compartments depending of positive boolean
    
    if positive_1 :
        for i in range(len(Z21)) :
            if Z21[i] == 0 :
                list_compartments1.append(0.0)
            if Z21[i] == 1 :
                list_compartments1.append(1.0)
                
    
    
    else :
        for i in range(len(Z21)) :
            if Z21[i] == 0 :
                list_compartments1.append(1.0)
            if Z21[i] == 1 :
                list_compartments1.append(0.0)
                
    ## Write in a txt file the compartments
    
    for x in list_filtered1 :
        list_compartments1.insert(x-1,-1.0)
    
    with open(compfile1,'w') as f :
        for x in list_compartments1 :
            f.write(str(x) + "\n")
            
    
    ## Do the HMM analysis for chr 2 : same thing
    
    compfile2 = HiCfile.replace(".RAWobserved","_comp_" + str(num_chr_2)+".txt")
    compfile2 = os.path.join(save_path,compfile2)
    
    remodel = hmm.GaussianHMM(n_components=2, covariance_type="diag", n_iter=1000)

    remodel.fit(corr_2)

    Z22 = remodel.predict(corr_2)
    list_compartments2 = []
    
    if positive_2 :
        for i in range(len(Z22)) :
            if Z22[i] == 0 :
                list_compartments2.append(0.0)
            if Z22[i] == 1 :
                list_compartments2.append(1.0)
    
    else :
        for i in range(len(Z22)) :
            if Z22[i] == 0 :
                list_compartments2.append(1.0)
            if Z22[i] == 1 :
                list_compartments2.append(0.0)

  
    
    for x in list_filtered2 :
        list_compartments2.insert(x-1,-1.0)

    with open(compfile2,'w') as f :
        for x in list_compartments2 :
            f.write(str(x) + "\n")
    
    
    


    


### APPLICATION ====================================================================================

mode = "inter"

if mode == "intra" :
    resolution = 100000
    HiC_fic = "chr16_100kb.RAWobserved"
    geneDen_fic = "chr16.hdf5"

    pipeline_intra(resolution,HiC_fic,geneDen_fic)

if mode == "inter" :

    resolution = 100000
    HiC_fic = "chr15_16_100kb.RAWobserved"
    geneDen_fic_1 = "chr15.hdf5"
    geneDen_fic_2 = "chr16.hdf5"

    pipeline_inter(resolution,HiC_fic,geneDen_fic_1,geneDen_fic_2)

