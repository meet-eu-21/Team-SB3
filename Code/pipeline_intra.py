#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 09:58:58 2021

@author: apauron
"""
import get_files_cluster ## To get the filename
import Compartments_SB3_cluster


"""
    A pipeline for generating intrachromosomal compartments in the cluster. 
    Keyword arguments :
    None

    
    Returns :
    
    None
    
 """

path_data_cluster = "/shared/projects/form_2021_21/trainers/dataforstudent/HiC/" ## Folder containing the data in cluster
a = get_files_cluster.getfiles(path_data_cluster,".RAWobserved") ## Get all raw data

list_files_intra = [] ## Will get all the filenames
list_resolutions = [] ## Will contain all resolutions
list_genes_density_files = [] ## Will contain all gene density files
for file in a :
	if "intrachromosomal" in file :
		list_files_intra.append(file) ## Get the file
        
        ## According to resolution append different resolutions
		if "25kb_resolution" in file :
			list_resolutions.append(25000)
			density_path = "/shared/projects/form_2021_21/trainers/dataforstudent/genedensity/25kb/"
			chr_name = file.split("/")[-1]
			chr_name = chr_name.split("_")[0]
			gdname = density_path +  chr_name + ".hdf5" ## Get density file name
			list_genes_density_files.append(gdname)
		
		
		if "100kb_resolution" in file :
			list_resolutions.append(100000)
			density_path = "/shared/projects/form_2021_21/trainers/dataforstudent/genedensity/100kb/"
			chr_name = file.split("/")[-1]
			chr_name = chr_name.split("_")[0]
			gdname = density_path +  chr_name + ".hdf5" ## Get density file name
			list_genes_density_files.append(gdname) 
		
	


for (filetocomp,resolution,gdfile) in zip(list_files_intra,list_resolutions,list_genes_density_files) :
    print(filetocomp)
    print(resolution)
    print(gdfile)
    Compartments_SB3_cluster.pipeline_intra(resolution,filetocomp,gdfile) ## Call the main pipeline
    

