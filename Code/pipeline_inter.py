#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 09:58:58 2021

@author: apauron
"""
import get_files_cluster ## To get the filename
import Compartments_SB3_cluster

"""
    A pipeline for generating interchromosomal compartments in the cluster. Mainly inspired by intra version    
    
    Keyword arguments :
    None

    
    Returns :
    
    None
    
 """


path_data_cluster = "/shared/projects/form_2021_21/trainers/dataforstudent/HiC/" ## data containing the HiCfiles
a = get_files_cluster.getfiles(path_data_cluster,".RAWobserved") ##Get all raw observed data

list_files_inter = [] ## Will contain all the interchromosomal files
list_genes_density_files = [] ## Will contain all the gene density files
for file in a :
	if "interchromosomal" in file : ## Get only interchromosomal files
		list_files_inter.append(file)
		density_path = "/shared/projects/form_2021_21/trainers/dataforstudent/genedensity/100kb/" ## Get the density file
		chr_name = file.split("/")[-1] ## The two chr files
        chr1_name = chr_name.split("_")[0]
        chr2_name = chr_name.split("_")[1]
		gdname1 = density_path +  chr1_name + ".hdf5" 
        gdname2 = density_path + chr2_name + ".hdf5"
		list_genes_density_files.append([gdname1,gdname2]) ## Add the two density files
				
	


for (filetocomp,gdfile) in zip(list_files_inter,list_genes_density_files) :
    print(filetocomp) # To get sure
    print(gdfile) # Idem
    Compartments_SB3_cluster.pipeline_inter(100000,filetocomp,gdfile[0],gdfile[1]) ## Call another function
    

