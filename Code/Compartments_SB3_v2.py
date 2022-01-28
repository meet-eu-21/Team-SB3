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
    return np.maximum(x,np.zeros(len(x)))

def relumin(x) :
    return np.minimum(x,np.zeros(len(x)))

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

def filteramat_Val(Hicmat,Filterextremum=True,factor=1.5):
    """
    in : a HiCmat without any transformation, factor of reduction
    out : the HiCmatreduce,thevector of his transformation
    THE filter part from the main in one function
    """

        
    Hicmatreduce=Hicmat
    
    A_0 = []
    for i in range(np.shape(Hicmat)[0]):
        if np.sum(Hicmat[i,:]) > 0:
            A_0.append(i)
    
    Hicmatreduce=Hicmatreduce[A_0,:]
    
    A_1 = []
    for j in range(np.shape(Hicmat)[1]):
        if np.sum(Hicmat[:,j]) > 0:
            A_1.append(j)
    
    Hicmatreduce=Hicmatreduce[:,A_1]


    #first step : filter empty bin
    sumHicmat0= Hicmatreduce.sum(axis = 0)
    segmenter0=sumHicmat0 > 0
    A_0 =np.where(segmenter0)
    Hicmatreduce=Hicmatreduce[:,A_0[1]]

    sumHicmat1 = Hicmatreduce.sum(axis = 1)
    segmenter1 = sumHicmat1 > 0
    A_1 = np.where(segmenter1)
    
    Hicmatreduce=Hicmatreduce[A_1[0],:]

    if Filterextremum:
        #second step : filter lower bin for first chromosome
        sumHicmat_0=np.sum(Hicmatreduce,axis=0)
        msum0=np.mean(sumHicmat_0)
        mstd0=np.std(sumHicmat_0)
        mini0 = msum0-mstd0*factor
        maxi0 = msum0+mstd0*factor
        #Make the bolean condition
        newcond0=mini0 < sumHicmat_0
        newcond20=sumHicmat_0 < maxi0
        newcond0=np.logical_and(newcond0,newcond20)
        B0=np.where(newcond0)

        #third step : filter lower bin for second chromosome
        sumHicmat_1=np.sum(Hicmatreduce,axis = 1)
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

###  Pipeline Intrachromosomique --------------------------------------------------------------------
def pipeline_intra(R,HiCfile,gene_density_file) :

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
    
    ### Path of my own folder in cluster ##
    
    save_path = "/shared/projects/form_2021_21/SB3/"
    
    
    
    
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
    
    
    ## Do the PCA analysis
    
        
    pca = PCA(n_components=2)
    pca.fit(corr)
    

    
    s_vector = pca.components_[0]
    list_vector = list(s_vector)
    
    vfile = HiCfile.replace(".RAWobserved","_vp.txt")
    vfile = vfile.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)

    ## Write the corresponding eigenvector in a file 

    list_filtered = [i for i in range(1,np.shape(binned_map)[0]) if i not in binsaved ]
    
    for x in list_filtered :
        list_vector.insert(x-1,-1.0)

    
    with open(vfile,'w') as f :
        for x in list_vector :
            f.write(str(x) + "\n")
    
    
   

    

        
    f = h5py.File(gene_density_file, 'r')
    data_name = list(f.keys())[0]
    dset = f[data_name]
    data = dset[...]

    positive = True
    if np.sum( data[np.where(s_vector>0)]) < np.sum( data[np.where(s_vector<0)]) :
        positive = False

    compfile = HiCfile.replace(".RAWobserved","_comp.txt")
    compfile = compfile.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)

    ## Write the comps
    ## Do the HMM analysis
    
    remodel = hmm.GaussianHMM(n_components=2, covariance_type="diag", n_iter=1000)

    remodel.fit(corr)

    Z2 = remodel.predict(corr)
    list_compartments = []
    
    if positive :
        for i in range(len(Z2)) :
            if Z2[i] == 0 :
                list_compartments.append(0.0)
            if Z2[i] == 1 :
                list_compartments.append(1.0)
    
    else :
        for i in range(len(Z2)) :
            if Z2[i] == 0 :
                list_compartments.append(1.0)
            if Z2[i] == 1 :
                list_compartments.append(0.0)
  
    
    for x in list_filtered :
        list_compartments.insert(x-1,-1.0)
    

    with open(compfile,'w') as f :
        for x in list_compartments :
            f.write(str(x) + "\n")
    
    z2 = np.zeros(len(s_vector))
    
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
    

    
    """
    ## Change our previous color_bins to get only the bins corresponding to filtered map and the column associated with our epiGmark selected
    color2=color_bins[binsaved[1]] #filter the epi by removed bin in HiC
    color2=color2[:,selectedmark] #now color2 is 1D
    color2=np.float64(color2.todense()) #type issue
    """
    #3D
    print('3D')#Here : sparse int64
    contact_map=HiCtoolbox.SCN(filtered_mat.copy())
    contact_map[contact_map==0] = 0.000000001
    contact_map=np.asarray(contact_map)**alpha #now we are not sparse at all
    
    print("begin fast floyd")
    dist_matrix = scipy.sparse.csgraph.floyd_warshall(1/contact_map,return_predecessors = False) #shortest path on the matrix
    print("fast floyd ended")
    dist_matrix=dist_matrix-np.diag(np.diag(dist_matrix))#remove the diagonal
    dist_matrix=(dist_matrix+np.transpose(dist_matrix))/2; #just to be sure that the matrix is symetric, not really usefull in theory
    

    
    #MDS
    embedding = MDS(n_components=3)#LOAD the MDS #With scikit-learn mds
    XYZ = embedding.fit_transform(dist_matrix) #Make the transform
    print("MDS transform done")
    #XYZ=np.float64(XYZ)
    #XYZ,E=HiCtoolbox.sammon(dist_matrix, 3)#with the one from tom j pollard
    
    #print("Sammon done")
    
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
    print("ConvexHull done")
    XYZ=XYZ*scale
    
    pdbfilename  = HiCfile.replace(".RAWobserved","_3D.pdb")
    pdbfilename = pdbfilename.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    HiCtoolbox.writePDB(pdbfilename,XYZ,list_compartments)
    print("PDB written")


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
    ### MODIFICATION VAL : Nouveau filtrage gérant les matrices non carrées

    #filtered_mat, binsaved = HiCtoolbox.filteramat(binned_map)
    
    
    filtered_mat, binsaved_Y, binsaved_X = filteramat_Val(binned_map)

    # --------------------------------------------------------------------------
    
    print('Input at the good resolution (filtered map) : ',np.shape(filtered_mat))

    
    ## Apply the SCN to the matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : SCN gérant les matrices non carrées

    #binned_map_scn = HiCtoolbox.SCN(filtered_mat)
    
    binned_map_scn = SCN_Val(filtered_mat)
    
    # --------------------------------------------------------------------------
        
    print('Input at the good resolution after SCN : ',np.shape(binned_map_scn))
#    print("NaN in matrix :",np.isnan(binned_map_scn).any())
#    print("Inf in matrix :",np.isinf(binned_map_scn).any())
#    print(binned_map_scn)
    
    ## Let us see a heatmap corresponding to the SCN binned_map :
    """
    
    hm_scn = sns.heatmap(
        binned_map_scn, 
        cmap= "YlOrRd",
        square=False,
        norm = LogNorm()
    )
    """
    ### Path of my own folder in cluster ##
    
    save_path = "/shared/home/apauron/"
    """
    heatmap = HiCfile.replace(".RAWobserved","_heatmap.png")
    heatmap = heatmap.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    filename =os.path.basename(heatmap)
    directory = heatmap.replace(filename,"")
    Path(directory).mkdir(parents=True, exist_ok=True)
    fig = hm_scn.get_figure()
    fig.savefig(heatmap,dpi = 400)
    
    plt.close()
    """
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

    print("Transformation into OE matrix :",np.shape(oe_mat))
    print("NaN in matrix :",np.isnan(oe_mat).any())
    print("Inf in matrix :",np.isinf(oe_mat).any())
    """
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
    """
    ## Let us compute the correlation matrix

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : Calcul de deux matrices de corrélation (forcément carrées),
    ###                     Une pour la OE_map et une pour sa transposée
    
    #corr = np.corrcoef(oe_mat)
    
    #hm_corr = sns.heatmap(corr, vmin = -0.1, vmax = 0.1, cmap= "coolwarm",square=False)  
    #fig = hm_corr.get_figure()
    #corr_heatmap = HiCfile.replace(".RAWobserved","_corr_heatmap.png")
    #corr_heatmap = corr_heatmap.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    #fig.savefig(corr_heatmap,dpi = 400)
    #plt.close()

    
    corr_1 = np.corrcoef(oe_mat)
    corr_2 = np.corrcoef(oe_mat,rowvar = False)
    """
    hm_corr_1 = sns.heatmap(corr_1, vmin = -0.1, vmax = 0.1, cmap= "coolwarm",square=False)  
    fig = hm_corr_1.get_figure()
    """
    corr_heatmap_1 = HiCfile.replace(".RAWobserved","_corr_heatmap_"+str(num_chr_1)+".csv")
    corr_heatmap_1 = corr_heatmap_1.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    """
    fig.savefig(corr_heatmap_1,dpi = 400)
    plt.close()
    """
    """
    hm_corr_2 = sns.heatmap(corr_2, vmin = -0.1, vmax = 0.1, cmap= "coolwarm",square=False)  
    fig = hm_corr_2.get_figure()
    """
    corr_heatmap_2 = HiCfile.replace(".RAWobserved","_corr_heatmap_"+str(num_chr_2)+".csv")
    corr_heatmap_2 = corr_heatmap_2.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    """
    fig.savefig(corr_heatmap_2,dpi = 400)
    plt.close()
    """
    # --------------------------------------------------------------------------

    print("Computation of Correlation Matrices :",np.shape(corr_1),'/',np.shape(corr_2))
    print("NaN in matrix 1 :",np.isnan(corr_1).any())
    if np.isnan(corr_1).any() :
        corr_1[np.isnan(corr_1)] = 0
    print("Inf in matrix 1 :",np.isinf(corr_1).any())
    print("NaN in matrix 2 :",np.isnan(corr_2).any())
    if np.isnan(corr_2).any() :
        corr_2[np.isnan(corr_2)] = 0
    print("Inf in matrix 2 :",np.isinf(corr_2).any())
    np.savetxt(corr_heatmap_1,corr_1)
    np.savetxt(corr_heatmap_2,corr_2)
    """
    ## Get the svd decomposition

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : Recherche des vecteurs et valeures propres des deux matrices de corrélation.
    
    #eigenvalues, eigenvectors = np.linalg.eig(corr)
    #eigenvectors = np.transpose(eigenvectors)
    #s_vector = eigenvectors[0]

    eigenvalues_1, eigenvectors_1 = np.linalg.eig(corr_1)
    eigenvectors_1 = np.transpose(eigenvectors_1)
    s_vector_1 = eigenvectors_1[0]

    eigenvalues_2, eigenvectors_2 = np.linalg.eig(corr_2)
    eigenvectors_2 = np.transpose(eigenvectors_2)
    s_vector_2 = eigenvectors_2[0]
    # --------------------------------------------------------------------------
    
    print("Computation of eigenvectors :",len(s_vector_1),'/',len(s_vector_2))

    ## Gene Density data treatment
    
    f1 = h5py.File(GDfilename_1, 'r')
    data_name_1 = list(f1.keys())[0]
    dset_1 = f1[data_name_1]
    data_1 = dset_1[...]
    data_1 = data_1[binsaved_X]

    plt.plot(np.arange(len(data_1)),data_1)
    plt.savefig("gene_density_"+str(num_chr_1)+".png")
    plt.close()
    # --------------

    f2 = h5py.File(GDfilename_2, 'r')
    data_name_2 = list(f2.keys())[0]
    dset_2 = f2[data_name_2]
    data_2 = dset_2[...]
    data_2 = data_2[binsaved_Y]

    plt.plot(np.arange(len(data_2)),data_2)
    plt.savefig("gene_density_"+str(num_chr_2)+".png")
    plt.close()

    print("Study the gene densities")

    ## Find the compartments

    # --------------------------------------------------------------------------
    ### MODIFICATION VAL : Étude de la positivité des Densité de Gènes de chaque chromosome.
    positive_1 = True
    if np.sum( data_1[np.where(s_vector_1>0)]) < np.sum( data_1[np.where(s_vector_1<0)]) :
        positive_1 = False

    print("First Gene Density positivity :",positive_1)
    
    z2_1 = np.zeros(len(s_vector_1))
    
    plt.plot(np.arange(len(s_vector_1)), s_vector_1)
    plt.xlabel("Index of the vector")
    plt.ylabel("First eigenvector value")
    plt.title("First eigenvector value for chromosome "+str(num_chr_1))
    plt.tight_layout()
    
    if positive_1 :
        plt.plot(np.arange(len(s_vector_1)), relu(s_vector_1),'r')

        plt.fill_between(np.arange(len(s_vector_1)), relu(s_vector_1), 0,
                         color='red', label = "Compartment A")
        plt.plot(np.arange(len(s_vector_1)), relumin(s_vector_1),'b')

        plt.fill_between(np.arange(len(s_vector_1)), relumin(s_vector_1), 0,
                         color='blue', label = "Compartment B")

        
        plt.legend()
    else :
        plt.plot(np.arange(len(s_vector_1)), relumin(s_vector_1),'r')
        plt.fill_between(np.arange(len(s_vector_1)), relumin(s_vector_1), 0,
                         color='red',label = "Compartment A")
        plt.plot(np.arange(len(s_vector_1)), relu(s_vector_1),'b')
        plt.fill_between(np.arange(len(s_vector_1)), relu(s_vector_1), 0,
                         color='blue',label = "Compartment B")
        plt.legend()

    nameplot = HiCfile.replace(".RAWobserved","_compartments_"+str(num_chr_1)+".png")
    nameplot = nameplot.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    
    plt.savefig(nameplot)
    plt.close()

    # --------------

    positive_2 = True
    if np.sum( data_2[np.where(s_vector_2>0)]) < np.sum( data_2[np.where(s_vector_2<0)]) :
        positive_2 = False

    print("Second Gene Density positivity :",positive_2)
    
    z2_2 = np.zeros(len(s_vector_2))
    
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
    nameplot = nameplot.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    
    plt.savefig(nameplot)
    plt.close()
    # --------------------------------------------------------------------------

    print("Find the compartments with HMM")

    list_filtered1 = [i for i in range(1,np.shape(binned_map)[0]) if i not in binsaved_X ]
    list_filtered2 = [i for i in range(1,np.shape(binned_map)[1]) if i not in binsaved_Y ]
    
    
    ## Do the HMM analysis for chr 1
    
    compfile1 = HiCfile.replace(".RAWobserved","_comp_" + str(num_chr_1)+".txt")
    compfile1 = compfile1.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    
    remodel = hmm.GaussianHMM(n_components=2, covariance_type="diag", n_iter=1000)

    remodel.fit(corr_1)

    Z21 = remodel.predict(corr_1)
    list_compartments1 = []
    
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
                
  
    for x in list_filtered1 :
        list_compartments1.insert(x-1,-1.0)
    
    with open(compfile1,'w') as f :
        for x in list_compartments1 :
            f.write(str(x) + "\n")
            
    
    ## Do the HMM analysis for chr 2
    
    compfile2 = HiCfile.replace(".RAWobserved","_comp_" + str(num_chr_2)+".txt")
    compfile2 = compfile2.replace("/shared/projects/form_2021_21/trainers/dataforstudent/HiC/",save_path)
    
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
    HiC_fic = "chr21_22_100kb.RAWobserved"
    geneDen_fic_1 = "chr21.hdf5"
    geneDen_fic_2 = "chr22.hdf5"

    pipeline_inter(resolution,HiC_fic,geneDen_fic_1,geneDen_fic_2)

