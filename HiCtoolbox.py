#python 3
#2019-2020
#CC-By-SA
#Carron Leopold
#specific utils tools for hic need

import h5py
import numpy as np
from scipy import sparse
import scipy


#####
#ADDITIONAL LOADER

def EpiGbyres(EpiGfilename,res,achr,sizeatres,NbEpi):
	"""
	Generate a matrix of repartition of epiG in a chr at the given resolution
	Do it the hard way : parse the file to the desired resolution directly
	"""
	Epimat=np.zeros((NbEpi,sizeatres))
	o=open(EpiGfilename,'r')
	l=o.readline()
	while l:
		ls=l.split()
		if ls[0]==achr:
			#force to scale bp cover in bin cover here
			EpiV=int(ls[3])-1
			b=np.float(ls[1])
			e=np.float(ls[2])
			begin=int(np.floor(b/res)) #force to cast for indexing error
			end=int(np.floor(e/res)) #force to cast for indexing error
			eC=end*res
			bC=begin*res
			val=1
			if begin==end: #1 bin
				Epimat[EpiV,begin]+=val*(e-b)
			else:
				if (end-begin)==1: #2 bin
					Epimat[EpiV,begin]+=val*(eC-b)
					Epimat[EpiV,end]+=val*(e-eC)
				else: #more than 2 bin
					Epimat[EpiV,begin]+=val*(res-(b-bC))
					while begin<end:
						Epimat[EpiV,begin]+=val*res
						begin+=1
					if (e-eC)>0:
						Epimat[EpiV,begin]+=val*(e-eC) #here begin=end
		l=o.readline()
	o.close()
	return Epimat

#####
#Binning tools
def bin2d(Data,p,q):    
    """   
    Data = input matrix
    p,q rescaling factors
    Written for sparse 
    """
    n,m=np.shape(Data);
    s=(int(np.ceil(n/p)),int(np.ceil(m/q)))
    i,j,d = sparse.find(Data);
    i=np.int_(np.ceil(i/p))
    j=np.int_(np.ceil(j/q))
    M=sparse.csr_matrix((d,(i,j)))
    return M



def bin1D(anumpyarray,resolutionfrom,resolutionto):
	"""
	in : A numpy array , number of bin in raw and in col
	out : the matrix binned
	"""
	print(resolutionto,resolutionfrom)
	if resolutionto>resolutionfrom:
		convertionfactor=np.ceil(resolutionto/resolutionfrom)
		s=anumpyarray.shape
		print("dimension du vecteur:",s)
		#has to be identical as result in other function like chrsizedict)
		newsizei=np.ceil(s[0]*resolutionfrom/resolutionto)
		newarray=np.zeros(int(newsizei))
		print("taille du vecteur appres rescale :",newarray.shape)
		i=0
		while i<newsizei:
			ifrom=int(i*convertionfactor)
			ito=int((i+1)*convertionfactor)
			if i==newsizei-1:
				asum=np.sum(anumpyarray[ifrom:])
			else:
				asum=np.sum(anumpyarray[ifrom:ito])
			newarray[i]=asum
			i+=1
		return newarray
	elif resolutionto==resolutionfrom:
		print("no binning")
		return anumpyarray
	else:
		print("wrong resolution parameter in bin1D")

def bin2dfullmat(anumpyarray,resolutionfrom,resolutionto):
	"""
	in : A numpy array , number of bin in raw and in col
	out : the matrix binned
	Written for full
	"""
	print('change of resolution from ',resolutionfrom,' to ',resolutionto)
	if resolutionto>resolutionfrom:
		convertionfactor=np.ceil(resolutionto/resolutionfrom)
		s=anumpyarray.shape
		print("Initial HiC size before binning:",s)
		#has to be identical as result in other function like chrsizedict)
		newsizei=np.ceil(s[0]*resolutionfrom/resolutionto)
		newsizej=np.ceil(s[1]*resolutionfrom/resolutionto)
		newarray=np.zeros((int(newsizei),int(newsizej)))
		print("HiC size after binning :",newarray.shape)
		i=0
		j=0
		while i<newsizei:
			while j<newsizej:
				ifrom=int(i*convertionfactor)
				ito=int((i+1)*convertionfactor)
				jfrom=int(j*convertionfactor)
				jto=int((j+1)*convertionfactor)
				if i==newsizei-1:
					asum=np.sum(anumpyarray[ifrom:,jfrom:jto])
				elif j==newsizej-1:
					asum=np.sum(anumpyarray[ifrom:ito,jfrom:])
				elif i==newsizei-1 and j==newsizej-1:
					asum=np.sum(anumpyarray[ifrom:,jfrom:])
				else:
					asum=np.sum(anumpyarray[ifrom:ito,jfrom:jto])
				newarray[i,j]=asum
				newarray[j,i]=asum
				j+=1
			i+=1
			j=0
		return newarray
	elif resolutionto==resolutionfrom:
		print("No binning")
		return anumpyarray
	else:
		print("Wrong resolution parameter")


#####
#operation tool

def SCN(D, max_iter = 10):
	"""
	Out  : SCN(D)
	Code version from Vincent Matthys
	"""    
	# Iteration over max_iter    
	for i in range(max_iter):        
		D /= np.maximum(1, D.sum(axis = 0))       
		D /= np.maximum(1, D.sum(axis = 1))    
		# To make matrix symetric again   
	return (D + D.T)/2 

def fastFloyd(contact):
	"""
	out : FF(contact)
	Code version from Vincent Matthys
	"""      
	n = contact.shape[0]    
	shortest = contact    
	for k in range(n):        
		i2k = np.tile(shortest[k,:], (n, 1))        
		k2j = np.tile(shortest[:, k], (n, 1)).T        
		shortest = np.minimum(shortest, i2k + k2j)    
	return shortest


def filteramat(Hicmat,Filterextremum=True,factor=1.5):
	"""
	in : a HiCmat without any transformation, factor of reduction
	out : the HiCmatreduce,thevector of his transformation
	THE filter part from the main in one function
	"""
	Hicmatreduce=Hicmat
	#first step : filter empty bin
	sumHicmat=Hicmat.sum(axis = 0)
	segmenter1=sumHicmat>0
	A=np.where(segmenter1)
	Hicmatreduce=Hicmatreduce[A[1],:]
	Hicmatreduce=Hicmatreduce[:,A[1]]
	if Filterextremum:
		#second step : filter lower bin
		sumHicmat=np.sum(Hicmatreduce,0)
		msum=np.mean(sumHicmat)
		mstd=np.std(sumHicmat)
		mini = msum-mstd*factor
		maxi = msum+mstd*factor
		#Make the bolean condition
		newcond=mini < sumHicmat
		newcond2=sumHicmat < maxi
		newcond=np.logical_and(newcond,newcond2)
		B=np.where(newcond)
		#Filter
		Hicmatreduce=Hicmatreduce[B[1],:]
		Hicmatreduce=Hicmatreduce[:,B[1]]
		segmenter1=A[1][B[1]] #Create the binsaved index
	return Hicmatreduce,segmenter1


def sammon(x, n, display = 2, inputdist = 'raw', maxhalves = 20, maxiter = 500, tolfun = 1e-9, init = 'default'):


    """Perform Sammon mapping on dataset x

    y = sammon(x) applies the Sammon nonlinear mapping procedure on
    multivariate data x, where each row represents a pattern and each column
    represents a feature.  On completion, y contains the corresponding
    co-ordinates of each point on the map.  By default, a two-dimensional
    map is created.  Note if x contains any duplicated rows, SAMMON will
    fail (ungracefully). 

    [y,E] = sammon(x) also returns the value of the cost function in E (i.e.
    the stress of the mapping).

    An N-dimensional output map is generated by y = sammon(x,n) .

    A set of optimisation options can be specified using optional
    arguments, y = sammon(x,n,[OPTS]):

       maxiter        - maximum number of iterations
       tolfun         - relative tolerance on objective function
       maxhalves      - maximum number of step halvings
       input          - {'raw','distance'} if set to 'distance', X is 
                        interpreted as a matrix of pairwise distances.
       display        - 0 to 2. 0 least verbose, 2 max verbose.
       init           - {'pca', 'cmdscale', random', 'default'}
                        default is 'pca' if input is 'raw', 
                        'msdcale' if input is 'distance'

    The default options are retrieved by calling sammon(x) with no
    parameters.

    File        : sammon.py
    Date        : 18 April 2014
    Authors     : Tom J. Pollard (tom.pollard.11@ucl.ac.uk)
                : Ported from MATLAB implementation by 
                  Gavin C. Cawley and Nicola L. C. Talbot

    Description : Simple python implementation of Sammon's non-linear
                  mapping algorithm [1].

    References  : [1] Sammon, John W. Jr., "A Nonlinear Mapping for Data
                  Structure Analysis", IEEE Transactions on Computers,
                  vol. C-18, no. 5, pp 401-409, May 1969.

    Copyright   : (c) Dr Gavin C. Cawley, November 2007.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

    """

    # Create distance matrix unless given by parameters
    if inputdist == 'distance':
        D = x
        if init == 'default':
            init = 'cmdscale'
    else:
        D = scipy.spatial.distance.cdist(x, x)
        if init == 'default':
            init = 'pca'

    if inputdist == 'distance' and init == 'pca':
        raise ValueError("Cannot use init == 'pca' when inputdist == 'distance'")

    if np.count_nonzero(np.diagonal(D)) > 0:
        raise ValueError("The diagonal of the dissimilarity matrix must be zero")

    # Remaining initialisation
    N = x.shape[0]
    scale = 0.5 / D.sum()
    D = D + np.eye(N)     

    if np.count_nonzero(D<=0) > 0:
        raise ValueError("Off-diagonal dissimilarities must be strictly positive")   

    Dinv = 1 / D
    if init == 'pca':
        [UU,DD,_] = np.linalg.svd(x)
        y = UU[:,:n]*DD[:n] 
    elif init == 'cmdscale':
        from cmdscale import cmdscale
        y,e = cmdscale(D)
        y = y[:,:n]
    else:
        y = np.random.normal(0.0,1.0,[N,n])
    one = np.ones([N,n])
    d = scipy.spatial.distance.cdist(y,y) + np.eye(N)
    dinv = 1. / d
    delta = D-d 
    E = ((delta**2)*Dinv).sum() 

    # Get on with it
    for i in range(maxiter):

        # Compute gradient, Hessian and search direction (note it is actually
        # 1/4 of the gradient and Hessian, but the step size is just the ratio
        # of the gradient and the diagonal of the Hessian so it doesn't
        # matter).
        delta = dinv - Dinv
        deltaone = np.dot(delta,one)
        g = np.dot(delta,y) - (y * deltaone)
        dinv3 = dinv ** 3
        y2 = y ** 2
        H = np.dot(dinv3,y2) - deltaone - np.dot(2,y) * np.dot(dinv3,y) + y2 * np.dot(dinv3,one)
        s = -g.flatten(order='F') / np.abs(H.flatten(order='F'))
        y_old    = y

        # Use step-halving procedure to ensure progress is made
        for j in range(maxhalves):
            s_reshape = np.reshape(s, (-1,n),order='F')
            y = y_old + s_reshape
            d = scipy.spatial.distance.cdist(y, y) + np.eye(N)
            dinv = 1 / d
            delta = D - d
            E_new = ((delta**2)*Dinv).sum()
            if E_new < E:
                break
            else:
                s = 0.5*s

        # Bomb out if too many halving steps are required
        if j == maxhalves-1:
            print('Warning: maxhalves exceeded. Sammon mapping may not converge...')

        # Evaluate termination criterion
        if abs((E - E_new) / E) < tolfun:
            if display:
                print('TolFun exceeded: Optimisation terminated')
            break

        # Report progress
        E = E_new
        if display > 1:
            print('epoch = %d : E = %12.10f'% (i+1, E * scale))

    if i == maxiter-1:
        print('Warning: maxiter exceeded. Sammon mapping may not have converged...')

    # Fiddle stress to match the original Sammon paper
    E = E * scale
    
    return [y,E]


#####
#Output WRITER


def writePDB(fnameout,value,EpiValue):
	"""
	write a PDB from value contain in value
	just bored to read 200 pages of biopython for a simple writing script
	I use tips from pierre poulain blog
	"""
	Sh=np.shape(value)
	#out
	fout=open(fnameout,'w')
	i=1
	while i<=Sh[0]:
		S="{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}" #PDB format in python
		#print("here",value[i-1,0],value[i-1,1],EpiValue[i-1])
		S=S.format('ATOM',i,'CA','','SOL','A',i,'',float(value[i-1,0]),float(value[i-1,1]),float(value[i-1,2]),1,float(EpiValue[i-1]),'O','')
		#print(S)
		fout.write(S+"\n")
		i+=1
	fout.close()

