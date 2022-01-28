#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 09:58:58 2021

@author: apauron
"""
import get_files_cluster ## To get the filename
import Compartments_SB3_cluster


path_data_cluster = "/shared/projects/form_2021_21/trainers/dataforstudent/HiC/"
a = get_files_cluster.getfiles(path_data_cluster,".RAWobserved")

list_files_intra = []
list_resolutions = []
list_genes_density_files = []
for file in a :
	if "intrachromosomal" in file :
		list_files_intra.append(file)
		if "25kb_resolution" in file :
			list_resolutions.append(25000)
			density_path = "/shared/projects/form_2021_21/trainers/dataforstudent/genedensity/25kb/"
			chr_name = file.split("/")[-1]
			chr_name = chr_name.split("_")[0]
			gdname = density_path +  chr_name + ".hdf5"
			list_genes_density_files.append(gdname)
		
		
		if "100kb_resolution" in file :
			list_resolutions.append(100000)
			density_path = "/shared/projects/form_2021_21/trainers/dataforstudent/genedensity/100kb/"
			chr_name = file.split("/")[-1]
			chr_name = chr_name.split("_")[0]
			gdname = density_path +  chr_name + ".hdf5"
			list_genes_density_files.append(gdname)
		
	


for (filetocomp,resolution,gdfile) in zip(list_files_intra,list_resolutions,list_genes_density_files) :
    print(filetocomp)
    print(resolution)
    print(gdfile)
    Compartments_SB3_cluster.pipeline_intra(resolution,filetocomp,gdfile)
    

