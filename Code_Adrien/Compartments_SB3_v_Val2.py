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

def SCN_Val(D, max_iter = 10):
    mat = np.float64(D.toarray())
    """
    Out  : SCN(D)
    Code version from Vincent Matthys
    """
    # Iteration over max_iter
    for i in range(max_iter):
        mat /= np.maximum(1, mat.sum(axis = 0))
        mat /= np.maximum(1, mat.sum(axis = 1)[:,None])
    # To make matrix symetric again
    return mat

###  Pipeline Intrachromosomique --------------------------------------------------------------------
def pipeline_intra(R,HiCfile,gene_density_file) :

    alpha=0.227
    HiCfilename= HiCfile # matrix containing the contact Hi-C
    
    #Build matrix
    A=np.loadtxt(HiCfilename) #Load the matrix
    A=np.int_(A) #To int
    print('Input data shape A : ',np.shape(A))
    
    ## Bin the matrix

    # 1 : Transpose la matrice A => ligne deviennent colonnes et vice-versa (matrice virtuelle A')
    # 2 : Échange les lignes 0 et 1 de la matrice A'
    # 3 : Transpose la matrice A' en A" <=> A" est une copie de la matrice A qui a échangé ses colonnes 0 et 1
    # 4 : Concaténation des matrices A et A" <=> A" se retrouve collée à A dans le sens des lignes
    # 5 : La nouvelle matrice A possède deux fois plus de lignes qu'avant
    A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)#build array at pb resolution
    A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
    binned_map=HiCtoolbox.bin2d(A,R,R) #!become csr sparse array
    LENTEST=np.shape(A)[0]
    print('Input at the good resolution (binned map) : ',np.shape(binned_map))
        
    ##Filter the matrix
    
    filtered_mat, binsaved = HiCtoolbox.filteramat(binned_map)
    print('Input at the good resolution (filtered mat) : ',np.shape(filtered_mat))
    
    
    
    
    ## Apply the SCN to the matrix
    
    binned_map_scn = HiCtoolbox.SCN(filtered_mat)
    
    print('Input at the good resolution after SCN : ',np.shape(binned_map_scn))
    
    
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
    data = data[binsaved]
    
    plt.plot(np.arange(len(data)),data)
    plt.savefig("gene_density.png")
    
    plt.close()
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
    
    """
        #Build color annotation at desired resolution
    color=pd.read_csv(EpiGfilename,delimiter='\t',header=None,names=[1,2,3,4]) ### 4 columns 
    color=color[color[1]=='chr16']#take only chr of interest
    number=color[4].max() #number of color in the file
    color_vec=np.zeros((LENTEST,number+1)) #build array at pb resolution LENchr * number of color
    i=0
    while i<np.shape(color)[0]:
    	color_vec[color[2].iloc[i]:color[3].iloc[i],color[4].iloc[i]]=1 ## iloc = gets the row
    	i+=1
        
    # At the end you get for each bp the epimark associated with a 1 at the column corresponding.
    # For ex, if the epi mark is 15 at bp 1000, the 1000th row has a 1 on the 16th column and 0 elsewhere.
      
    
    color_bins=HiCtoolbox.bin2d(color_vec,R,1) ## keep 16 epimarks at col but is binned with given resolution
    color_bins=color_bins/np.amax(color_bins) ##To normalize
    
    #The score corresponding to the "density" of each epiGmark"
    print('Bp cover by this mark, has to be >0 :',np.sum(color_bins[:,selectedmark]) )
    
    """
    
    # FILTER like the HiC tool box
    print("before filtering : ",np.shape(binned_map))

    
    print("after filtering : ",np.shape(filtered_mat))#,np.shape(color_vecseg))
    
    """
    ## Change our previous color_bins to get only the bins corresponding to filtered map and the column associated with our epiGmark selected
    color2=color_bins[binsaved[1]] #filter the epi by removed bin in HiC
    color2=color2[:,selectedmark] #now color2 is 1D
    color2=np.float64(color2.todense()) #type issue
    """
    #3D
    print('3D')#Here : sparse int64
    contact_map=HiCtoolbox.SCN(filtered_mat.copy()) 
    contact_map=np.asarray(contact_map)**alpha #now we are not sparse at all
    dist_matrix = HiCtoolbox.fastFloyd(1/contact_map) #shortest path on the matrix
    dist_matrix=dist_matrix-np.diag(np.diag(dist_matrix))#remove the diagonal
    dist_matrix=(dist_matrix+np.transpose(dist_matrix))/2; #just to be sure that the matrix is symetric, not really usefull in theory
    
    
    #MDS
    #embedding = MDS(n_components=3)#LOAD the MDS #With scikit-learn mds
    #XYZ = embedding.fit_transform(dist_matrix) #Make the transform
    #XYZ=np.float64(XYZ)
    XYZ,E=HiCtoolbox.sammon(dist_matrix, 3)#with the one from tom j pollard
    
    
    
    """
    print("Output shape : ",np.shape(XYZ),np.shape(color2))
    """
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


###  Pipeline Interchromosomique --------------------------------------------------------------------
    
def pipeline_inter(R,HiCfile,gene_density_file_1,gene_density_file_2) :

    alpha=0.227
    HiCfilename= HiCfile # matrix containing the contact Hi-C
    GDfilename_1 = gene_density_file_1
    GDfilename_2 = gene_density_file_2
    num_chr_1 = GDfilename_1.split('chr')[1].split('.')[0]
    num_chr_2 = GDfilename_2.split('chr')[1].split('.')[0]
    
    #Build matrix
    A=np.loadtxt(HiCfilename) #Load the matrix
    A=np.int_(A) #To int
    print('Input data shape A : ',np.shape(A)) # Nombre de lignes * nombre de colonnes du fichier HiCfilename
    ## Bin the matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : On ne concatène pas la matrice A et sa transposée car :
    ###     1) Ça force une nature carrée de la matrice alors que les chromosomes ne sont pas forcément de même longueur
    ###     2) Il y a ambigüité des données lorsque deux couples complémentaires de coordonnées n'ont pas la même valeur
    ###         (ex: si les coords (0,1)=5 et (1,0)=10, la concaténation des transposées crééra
    ###                 des coords (1,0)=5 et (0,1)=10 pour la matrice finale)
    
    #A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)#build array at pb resolution
    #print('Input data shape A post concatene : ',np.shape(A))
    #print(A)
    # --------------------------------------------------------------------------
    A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
    print('Input data shape A post coo_matrix : ',np.shape(A))
    binned_map=HiCtoolbox.bin2d(A,R,R) #!become csr sparse array
    LENTEST=np.shape(A)[0]
    print('Input at the good resolution (binned map) : ',np.shape(binned_map))
#    print(binned_map)
#    print(type(binned_map))
        
    ##Filter the matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : On ne filtre pas la matrice binned_map car ça force une nature carrée

    filtered_mat, binsaved = HiCtoolbox.filteramat(binned_map)
    print(type(filtered_mat))
    print('Input at the good resolution (filtered map) : ',np.shape(filtered_mat))
    # --------------------------------------------------------------------------
    
    
    
    ## Apply the SCN to the matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : SCN gérant les matrices non carrées et appliquée à la matrice non filtrée

    #binned_map_scn = HiCtoolbox.SCN(filtered_mat)

    binned_map_scn = SCN_Val(binned_map)
    # --------------------------------------------------------------------------
        
    print('Input at the good resolution after SCN : ',np.shape(binned_map_scn))
    print("NaN in matrix :",np.isnan(binned_map_scn).any())
    print("Inf in matrix :",np.isinf(binned_map_scn).any())
#    print(binned_map_scn)
#    print(np.isnan(binned_map_scn).any())
#    print(np.isinf(binned_map_scn).any())
    
    ## Let us see a heatmap corresponding to the SCN binned_map :
    
    
    hm_scn = sns.heatmap(
        binned_map_scn, 
        cmap= "YlOrRd",
        square=False,
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
    
    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : Création d'une OE matrix non carrée
    oe_mat = binned_map_scn.copy()
    #n = np.shape(oe_mat)[0]
    n,m = np.shape(oe_mat)
    #list_diagonals_mean = [np.mean(oe_mat.diagonal(i)) for i in range(-n+1,n)]
    list_diagonals_mean = [np.mean(oe_mat.diagonal(i)) for i in range(-n+1,m)]
    
    for i in range(n):
        #for j in range(n) :
        for j in range(m) :
            if list_diagonals_mean[j-i+n-1] == 0 :
                oe_mat[i,j] = 0
            else :       
                oe_mat[i,j] /= list_diagonals_mean[j-i+n-1]
    # --------------------------------------------------------------------------

    print("Transformation into OE matrix")
    print("NaN in matrix :",np.isnan(oe_mat).any())
    print("Inf in matrix :",np.isinf(oe_mat).any())
    
    hm_oe = sns.heatmap(
        oe_mat, 
        cmap= "coolwarm",
        square=False,
        norm = LogNorm()
    )  
    
    oe_heatmap = HiCfile.replace(".RAWobserved","_oe_heatmap.png")
    oe_heatmap = oe_heatmap.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
  
    fig = hm_oe.get_figure()
    fig.savefig(oe_heatmap,dpi = 400)
    
    plt.close()
    
    ## Let us compute the correlation matrix
    
    corr = np.corrcoef(oe_mat)
    
    print("Computation of Correlation Matrix")
    print("NaN in matrix :",np.isnan(corr).any())
    print("Inf in matrix :",np.isinf(corr).any())

    
    hm_corr = sns.heatmap(
        corr,
        vmin = -0.1, vmax = 0.1,
        cmap= "coolwarm",
        square=False,
    )  
    
    fig = hm_corr.get_figure()
    
    corr_heatmap = HiCfile.replace(".RAWobserved","_corr_heatmap.png")
    corr_heatmap = corr_heatmap.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
 
    fig.savefig(corr_heatmap,dpi = 400)
    
    plt.close()
    
    """
    ## Get the svd decomposition
    
    eigenvalues, eigenvectors = np.linalg.eig(corr)
    
    eigenvectors = np.transpose(eigenvectors)
    
    s_vector = eigenvectors[0]
    
    # ----------------------------------------------------------------------------
    ### MODIFICATION VAL : Traitement de deux fichiers Gene Density
    # ----------------------------------------------------------------------------
    f1 = h5py.File(GDfilename_1, 'r')
    data_name_1 = list(f1.keys())[0]
    dset_1 = f1[data_name_1]
    data_1 = dset_1[...]
    print("D1 :",len(data_1))
    data_1 = data_1[binsaved]

    plt.plot(np.arange(len(data_1)),data_1)
    plt.savefig("gene_density_"+str(num_chr_1)+".png")
    plt.close()
    
    f2 = h5py.File(GDfilename_2, 'r')
    data_name_2 = list(f2.keys())[0]
    dset_2 = f2[data_name_2]
    data_2 = dset_2[...]
    print("D2 :",len(data_2))
    data_2 = data_2[binsaved]


    plt.plot(np.arange(len(data_2)),data_2)
    plt.savefig("gene_density_"+str(num_chr_2)+".png")
    plt.close()
    
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
    
    
    # FILTER like the HiC tool box
    print("before filtering : ",np.shape(binned_map))

    
    print("after filtering : ",np.shape(filtered_mat))#,np.shape(color_vecseg))
    
    #3D
    print('3D')#Here : sparse int64
    contact_map=HiCtoolbox.SCN(filtered_mat.copy()) 
    contact_map=np.asarray(contact_map)**alpha #now we are not sparse at all
    dist_matrix = HiCtoolbox.fastFloyd(1/contact_map) #shortest path on the matrix
    dist_matrix=dist_matrix-np.diag(np.diag(dist_matrix))#remove the diagonal
    dist_matrix=(dist_matrix+np.transpose(dist_matrix))/2; #just to be sure that the matrix is symetric, not really usefull in theory
    
    
    #MDS
    #embedding = MDS(n_components=3)#LOAD the MDS #With scikit-learn mds
    #XYZ = embedding.fit_transform(dist_matrix) #Make the transform
    #XYZ=np.float64(XYZ)
    XYZ,E=HiCtoolbox.sammon(dist_matrix, 3)#with the one from tom j pollard
    
    
    

#    print("Output shape : ",np.shape(XYZ),np.shape(color2))

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

### APPLICATION ====================================================================================

mode = "inter"

if mode == "intra" :
    resolution = 100000
    HiC_fic = "chr16_100kb.RAWobserved"
    geneDen_fic = "chr16.hdf5"

    pipeline_intra(resolution,HiC_fic,geneDen_fic)

if mode == "inter" :

    resolution = 100000
    HiC_fic = "chr12_16_100kb.RAWobserved"
    geneDen_fic_1 = "chr12.hdf5"
    geneDen_fic_2 = "chr16.hdf5"

    pipeline_inter(resolution,HiC_fic,geneDen_fic_1,geneDen_fic_2)

