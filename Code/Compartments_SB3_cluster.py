import os
import numpy as np
import h5py
import scipy
from matplotlib.colors import LogNorm, Normalize
from scipy.spatial import ConvexHull
from scipy import sparse
import matplotlib.pyplot as plt
import seaborn as sns
import HiCtoolbox
from hmmlearn import hmm
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.manifold import MDS



    
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
    
    ## Generate two text files : one for compartments and the other for the first eigenvector
    
    
                


