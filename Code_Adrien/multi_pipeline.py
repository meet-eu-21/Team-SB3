#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 09:58:58 2021

@author: apauron
"""
import get_files ## To get the filename
import Compartments_SB3

url = 'http://www.lcqb.upmc.fr/meetu/dataforstudent/HiC/'
a = get_files.approfondissement_url(url)

path_data_cluster = "/shared/projects/form_2021_21/trainers/dataforstudent/HiC/"
list_files_intra = []
list_resolutions = []

for file in a :
    if "intrachromosomal" in file :
        list_files_intra.append(file.replace("http://www.lcqb.upmc.fr/meetu/dataforstudent/HiC/",path_data_cluster))
    if "25kb_resolution" in file :
        list_resolutions.append(25000)
    if "100kb_resolution" in file :
        list_resolutions.append(100000)

for (filetocomp,resolution) in zip(list_files_intra,list_resolutions) :
    Compartments_SB3.pipeline(resolution,filetocomp,"E116_15_coreMarks_dense.txt")
    

        
    


 ## Get for all the data
    
