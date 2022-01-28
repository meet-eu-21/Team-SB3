#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 09:03:17 2022

@author: apauron
"""
import os
import get_files_cluster
import pandas as pd
from numpy import genfromtxt

### Get the parent folder of the working directory. Change it if you modify the name of the folders
path_parent = os.path.dirname(os.getcwd()) 
path_SB1 = os.path.join(path_parent,"Results_SB1_intra") #location of SB1 intrachromosomal results to convert


folder_results = "Results_Intra" 
path_SB3 = os.path.join(path_parent,folder_results) #location of SB3 intrachromosomal results to convert
list_chr = os.listdir(os.path.join(path_parent,folder_results,"HUVEC","25kb_resolution_intrachromosomal")) ## All the chromosomes

###################################################Convert SB1 results to SB3 results##########################################

def SB1toSB3(path):
    
    """
    A pipeline to convert SB1 generated compartments files into SB3 format.  
    
    Keyword arguments :
    path -- the path containing the folder in which there are the files containing SB1 results
  
    Returns :
    

    all the converted SB1 results in SB3 format in the "SB1_converted_SB3" folder
    
    """
    
    filestoconvert = get_files_cluster.getfiles(path,"") #get all the files in the path
    for file in filestoconvert :
        cell_type = file.split("/")[-1].split("_")[0]
        for resolution in ["25kb","100kb"] :
            if resolution in file :
                df_file = pd.read_csv(file,sep = ' ',header = None) #get the SB1 file
                df_file["chrname"] = df_file[0] + df_file[1].astype(str) #transform chr x to chrx
                df_file["comp"] = df_file[4] #get the comp number
                df_file = df_file[["chrname","comp"]] #because SB3 type file is only chrname and comp
                chr_values = pd.unique(df_file.chrname) #get the chr values
                grouped = df_file.groupby(df_file.chrname) #to split according to chr name
                for chr in chr_values :
                    split_df = grouped.get_group(chr)
                    split_df.comp = split_df.comp.replace([-1.0,0.0],[0.0,-1.0]) ## Change the format of filtered and B compartment bins
                    if not os.path.exists(os.path.join(path_parent,"SB1_converted_SB3",cell_type,resolution)): #Create folder if not exists
                            os.makedirs(os.path.join(path_parent,"SB1_converted_SB3",cell_type,resolution))
                    filename = os.path.join(path_parent,"SB1_converted_SB3",cell_type,resolution,chr + "_" + resolution +  "_comp.txt") 
                    split_df.comp.to_csv(filename,header = False, index = False) #create the files corresponding to our metric
                    

###################################################Convert SB3 results to SB1 results##########################################

def SB3toSB1(path):
    
    """
    A pipeline to convert SB3 generated compartments files into SB1 format.  
    
    Keyword arguments :
    path -- the path containing the folder in which there are the files containing SB1 results
  
    Returns :
    

    all the converted SB3 results in SB1 format in the "SB3_converted_SB1" folder
    
    """
    
    files_results = get_files_cluster.getfiles(path,"comp") #get files inside the path given

    for resolution in ["25kb","100kb"] : ## Because those are intrachromosomal results
        
        for cell_type in os.listdir(os.path.join(path_parent,folder_results)): ## adapt if not all cell types are present
            if os.path.isdir(os.path.join(path_parent,folder_results,cell_type)):
   
                list_df = []
                
                for chr in list_chr : ## List all the chromosomes
                    
                    for file_results in files_results :
                        
                        # find the good corresponding file to chr,cell_type and results
                        if chr in file_results and cell_type in file_results and resolution in file_results : 
                            
                            file_df = pd.DataFrame()
                            
                            # Transformation into a SB1 type file : chr x start end comp
                            
                            lresults = genfromtxt(file_results, delimiter='\n') 
                            file_df["comp"] = lresults 
                            file_df["chromosome"] = ["chr" for i in range(len(lresults))]
                            file_df["chrnum"] = [chr.replace("chr","") for i in range(len(lresults))]
                            
                            #According to resolution, create the start and end bins
                            
                            if resolution == "100kb" :
                                file_df["start"] = [100000.0*x for x in file_df.index.tolist()]
                            else :
                                file_df["start"] = [25000.0*x for x in file_df.index.tolist()]
                            if resolution == "100kb" :
                                file_df["end"] = [100000.0*(x+1) for x in file_df.index.tolist()]
                            else :
                                file_df["end"] = [25000.0*(x+1) for x in file_df.index.tolist()]
                            
                            #Append to a list the dataframe corresponding to the chromosome
                            
                            file_df_copy = file_df.copy()
                            file_df_copy = file_df_copy[["chromosome","chrnum","start","end","comp"]]
                            file_df_copy.comp[file_df.comp == 0.0] = -1.0
                            file_df_copy.comp[file_df.comp == -1.0] = 0.0
                            list_df.append(file_df_copy)
               
                #Concatenate all the dataframes with chromosomes of the same cell type
                
                res_df = pd.concat(list_df)
                res_df = res_df.sort_values(by = ["chrnum","start"])
                filename = os.path.join(path_parent,"SB3_converted_SB1",cell_type + "_" + resolution + "_COMPARTMENT" )
                res_df.to_csv(filename,header = False, index = False, sep = " ")
                        


