#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 14:17:28 2022

@author: apauron
"""

import os
import get_files_cluster
import pandas as pd
from numpy import genfromtxt
import numpy as np

### Create a csv file containing for each chromosome and cellular type the similarity with gold standard. Do it for 100kb and 25kb intra.

### Get the parent folder of the working directory
path_parent = os.path.dirname(os.getcwd())

## Name of the folder containing the data
folder_gold_standards = "Gold_standards"
folder_results = "Results_intra"

files_gold_standard = get_files_cluster.getfiles(os.path.join(path_parent,folder_gold_standards),"comp")
files_results = get_files_cluster.getfiles(os.path.join(path_parent,folder_results),"comp")
list_chr = os.listdir(os.path.join(path_parent,folder_results,"HUVEC","25kb_resolution_intrachromosomal"))


for resolution in ["25kb","100kb"] :
    res_df = pd.DataFrame()
    
    for cell_type in os.listdir(os.path.join(path_parent,folder_gold_standards)):
        list_similarity = []

        for chr in list_chr :
            
            
            results = ""
            gold = ""
            
            for file_gold,file_results in zip(files_gold_standard,files_results) :
                
                if chr in file_results and cell_type in file_results and resolution in file_results :
                    
                    results = file_results
                
                if chr in file_gold and cell_type in file_gold and resolution in file_gold :
                    
                    gold = file_gold
            
            lresults = genfromtxt(results, delimiter='\n')
            lgold = genfromtxt(gold, delimiter='\n')
            
            if np.shape(lresults)[0] > np.shape(lgold)[0] :
                lresults = np.copy(lresults[:np.shape(lgold)[0]])
            
            else :
                lgold = np.copy(lgold[:np.shape(lresults)[0]])
        
            compare_mat = np.equal(lresults,lgold)
            
            count = np.sum(compare_mat)
            list_similarity.append(count/min(np.shape(lresults)[0],np.shape(lgold)[0]))
        
        res_df[cell_type] = np.round(list_similarity,2)
        res_df.index = list_chr
    
    resfilename = resolution + "_similarity.csv"
    res_df.to_csv(resfilename)
            

            
        
        