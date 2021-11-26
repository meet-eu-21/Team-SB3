import sys
import os
import numpy as np
import scipy
from scipy.spatial import ConvexHull
from scipy import sparse
import pandas as pd
from sklearn.manifold import MDS #If you want the scikit learn mds
import HiCtoolbox

os.chdir("D:/Utilisateurs/Adrien/Documents/GitHub/Team-SB3")
R=100000
NbmaxEpi=15 #Epi states go from 0 to 15
alpha=0.227
selectedmark=1 #index of the selected mark
HiCfilename='chr16_100kb.RAWobserved' # matrix containing the contact Hi-C
EpiGfilename='E116_15_coreMarks_dense.txt' #EpiGenetic file associating to each bin an epigenetic type

#Build matrix
A=np.loadtxt(HiCfilename) #Load the matrix
A=np.int_(A) #To int
print('Input data shape : ',np.shape(A))
A=np.concatenate((A,np.transpose(np.array([A[:,1],A[:,0],A[:,2]]))), axis=0)#build array at pb resolution
A = sparse.coo_matrix( (A[:,2], (A[:,0],A[:,1])))
binned_map=HiCtoolbox.bin2d(A,R,R) #!become csr sparse array
LENTEST=np.shape(A)[0]
print('Input at the good resolution : ',np.shape(binned_map))
