#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 14:17:28 2022

@author: apauron
"""

import os
import get_files_cluster
import pandas as pd
import re
from numpy import genfromtxt
import numpy as np

### Create a csv file containing for each chromosome and cellular type the similarity with gold standard. Do it for 100kb and 25kb intra.

### Get the parent folder of the working directory
path_parent = os.path.dirname(os.getcwd())

## Name of the folder containing the data
folder_gold_standards = "Gold_standards"
folder_results = "Results_Intra"
SB1_results = "SB1_converted_SB3"

def get_results_intra(result_folder):

    #Get all the files in the results folder and the gold standard folder

    files_gold_standard = get_files_cluster.getfiles(os.path.join(path_parent,folder_gold_standards),"comp")
    files_results = get_files_cluster.getfiles(os.path.join(path_parent,result_folder),"comp")
    #Get the list of the chromosomes thanks to one of the folder
    
    list_chr = os.listdir(os.path.join(path_parent,folder_results,"HUVEC","25kb_resolution_intrachromosomal"))
    
    #For each resolution, we create a csv file giving for each cell type and chromosome the percentage of same compartments.
    
    for resolution in ["25kb","100kb"] :
        res_df = pd.DataFrame()
        
        #For each cell type, chromosome and resolution
        for cell_type in os.listdir(os.path.join(path_parent,result_folder)):
            list_similarity = []
            print(cell_type)
            if os.path.isdir(os.path.join(path_parent,result_folder,cell_type)):
                for chr in list_chr :
                    
                    
                    results = ""
                    gold = ""
                    
                    #Get the same results files and gold standards files
                    
                    for file_gold,file_results in zip(files_gold_standard,files_results) :
                        
                        if chr in file_results and cell_type in file_results and resolution in file_results :
                            
                            results = file_results
                        
                        if chr in file_gold and cell_type in file_gold and resolution in file_gold :
                            
                            gold = file_gold
                    
                    lresults = genfromtxt(results, delimiter='\n')
                    lgold = genfromtxt(gold, delimiter='\n')
                    
                    #The two files have some shape difference, so we slice at the end the longest one
                    
                    if np.shape(lresults)[0] > np.shape(lgold)[0] :
                        lresults = np.copy(lresults[:np.shape(lgold)[0]])
                    
                    else :
                        lgold = np.copy(lgold[:np.shape(lresults)[0]])
                        
                    #We try similarity if compartment A = 0.0 or 1.0, because we don't know the corresponding compartment in gold standard    
                    
                    inv_results = np.ones(np.shape(lresults)[0],) - lresults # It changes 0 into 1 and 1 into 0, -1 into 2
                    inv_results = np.where(inv_results != 2.0, inv_results, -1.0) #We change back 2 into -1.0
                
                    compare_mat = np.equal(lresults,lgold) #Compare for each comp if it is equal to gold standard
                    
                    count = np.sum(compare_mat)
                    sim = count/min(np.shape(lresults)[0],np.shape(lgold)[0]) # Get the percentage
                    
                    compare_mat_inv = np.equal(inv_results,lgold) #Compare for each comp if it is equal to gold standard
                    count_int = np.sum(compare_mat_inv)
                    sim_inv = count_int/min(np.shape(lresults)[0],np.shape(lgold)[0]) # Get the percentage
                    
                    #We get the best percentage with comp A = 0.0 or 1.0
                    
                    list_similarity.append(max(sim,sim_inv))
                    
                    
                
                res_df[cell_type] = np.round(list_similarity,2) #We round to two figures to get a percentage
                res_df.index = list_chr #as index, we get the list of chromosomes
        
        #Write the file
        resfilename = resolution + "_similarity.csv"
        resfilename = os.path.join(path_parent,result_folder,resfilename)
        res_df.to_csv(resfilename)
            

def get_results_inter(folder):  
    #Get all the files in the results folder and the gold standard folder

    files_gold_standard = get_files_cluster.getfiles(os.path.join(path_parent,folder_gold_standards),"comp")
    files_results = get_files_cluster.getfiles(os.path.join(path_parent,folder),"comp")
    #Get the list of the chromosomes thanks to one of the folder
    
    list_chr = [str(i) for i in range(1,23)]
    list_chr.append("X")
             
    #For each cell type, we create a file of comparison for 100 kb
    for cell_type in os.listdir(os.path.join(path_parent,folder)):
        res_df = pd.DataFrame()

        if os.path.isdir(os.path.join(path_parent,folder,cell_type)):
            for chr1 in list_chr :
                list_similarity = []

                list_2 = list_chr[:]
                list_2.remove(chr1)
                for chr2 in list_2 :                
                
                    results = ""
                    gold = ""
                    
                    #Let's construct the good regex
                    
                    regex = "chr({}|{})_({}|{})_100kb_comp_{}".format(chr1,chr2,chr1,chr2,chr2)
                   
                    #Get the same results files and gold standards files
                    
                    for file_gold in files_gold_standard:
                        
                        if "chr" + chr2 + "_" in file_gold and cell_type in file_gold and "100kb" in file_gold :
                            
                            gold = file_gold
                     
                    for file_results in files_results:
                        matched = re.search(regex,file_results)

                        if matched and cell_type in file_results :
                            
                            results = file_results
                        
                    
                    print(results,gold)
                    
                    lresults = genfromtxt(results, delimiter='\n')
                    lgold = genfromtxt(gold, delimiter='\n')
                    
                    #The two files have some shape difference, so we slice at the end the longest one
                    
                    if np.shape(lresults)[0] > np.shape(lgold)[0] :
                        lresults = np.copy(lresults[:np.shape(lgold)[0]])
                    
                    else :
                        lgold = np.copy(lgold[:np.shape(lresults)[0]])
                        
                    #We try similarity if compartment A = 0.0 or 1.0, because we don't know the corresponding compartment in gold standard    
                    
                    inv_results = np.ones(np.shape(lresults)[0],) - lresults # It changes 0 into 1 and 1 into 0, -1 into 2
                    inv_results = np.where(inv_results != 2.0, inv_results, -1.0) #We change back 2 into -1.0
                
                    compare_mat = np.equal(lresults,lgold) #Compare for each comp if it is equal to gold standard
                    
                    count = np.sum(compare_mat)
                    sim = count/min(np.shape(lresults)[0],np.shape(lgold)[0]) # Get the percentage
                    
                    compare_mat_inv = np.equal(inv_results,lgold) #Compare for each comp if it is equal to gold standard
                    count_int = np.sum(compare_mat_inv)
                    sim_inv = count_int/min(np.shape(lresults)[0],np.shape(lgold)[0]) # Get the percentage
                    
                    #We get the best percentage with comp A = 0.0 or 1.0
                    
                    list_similarity.append(max(sim,sim_inv))
                
                
                list_similarity.insert(list_chr.index(chr1),-9999)
                res_df["chr" + chr1] = np.round(list_similarity,2) #We round to two figures to get a percentage
                list_index = ["chr" + c for c in list_chr]
                res_df.index = list_index #as index, we get the list of chromosomes
                
    
            #Write the file
            resfilename = cell_type + "_similarity_inter.csv"
            resfilename = os.path.join(path_parent,folder,resfilename)
            res_df.to_csv(resfilename)        
        
        