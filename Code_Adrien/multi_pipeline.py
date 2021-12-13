#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 09:58:58 2021

@author: apauron
"""
import get_files_cluster ## To get the filename
import Compartments_SB3_cluster


path_data_cluster = "/shared/projects/form_2021_21/trainers/dataforstudent/HiC/"
a = get_files_cluster.getfiles(path_data_cluster)

list_files_intra = []
list_resolutions = []

for file in a :
    if "intrachromosomal" in file :
        list_files_intra.append(file)
    if "25kb_resolution" in file :
        list_resolutions.append(25000)
    if "100kb_resolution" in file :
        list_resolutions.append(100000)

for (filetocomp,resolution) in zip(list_files_intra,list_resolutions) :
    print(filetocomp)
    print(resolution)
    Compartments_SB3_cluster.pipeline(resolution,filetocomp,"E116_15_coreMarks_dense.txt")
    


