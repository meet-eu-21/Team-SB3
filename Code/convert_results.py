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

### Get the parent folder of the working directory
path_parent = os.path.dirname(os.getcwd())
path_SB1 = os.path.join(path_parent,"Results_SB1_intra")


folder_results = "Results_Intra"
path_SB3 = os.path.join(path_parent,folder_results)


def SB1toSB3(path):
    filestoconvert = get_files_cluster.getfiles(path,"")
    for file in filestoconvert :
        for resolution in ["25kb","100kb"] :
            if resolution in file :
                df_file = pd.read_csv(file,sep = ' ',header = None)
                df_file["chrname"] = df_file[0] + df_file[1].astype(str)
                df_file["comp"] = df_file[4]
                df_file = df_file[["chrname","comp"]]
                chr_values = pd.unique(df_file.chrname)
                grouped = df_file.groupby(df_file.chrname)
                for chr in chr_values :
                    split_df = grouped.get_group(chr)
                    newsplit = split_df.copy()
                    newsplit.comp[split_df.comp == 0.0] = -1.0
                    newsplit.comp[split_df.comp == -1.0] = 0.0
                    filename = os.path.join(path_parent,"SB1_converted_SB3",chr + "_" + resolution + "_comp.txt")
                    split_df.comp.to_csv(filename,header = False, index = False)
                    
list_chr = os.listdir(os.path.join(path_parent,folder_results,"HUVEC","25kb_resolution_intrachromosomal"))


def SB3toSB1(path):
    files_results = get_files_cluster.getfiles(path,"comp")

    for resolution in ["25kb","100kb"] :
        
        for cell_type in os.listdir(os.path.join(path_parent,folder_results)):
            
            list_df = []
            
            for chr in list_chr :
                
                for file_results in files_results :
                    
                    if chr in file_results and cell_type in file_results and resolution in file_results :
                        
                        file_df = pd.DataFrame()
                        
                        lresults = genfromtxt(file_results, delimiter='\n')
                        file_df["comp"] = lresults
                        file_df["chromosome"] = ["chr" for i in range(len(lresults))]
                        file_df["chrnum"] = [chr.replace("chr","") for i in range(len(lresults))]
                        if resolution == "100kb" :
                            file_df["start"] = [100000.0*x for x in file_df.index.tolist()]
                        else :
                            file_df["start"] = [25000.0*x for x in file_df.index.tolist()]
                        if resolution == "100kb" :
                            file_df["end"] = [100000.0*(x+1) for x in file_df.index.tolist()]
                        else :
                            file_df["end"] = [25000.0*(x+1) for x in file_df.index.tolist()]
                        
                        file_df_copy = file_df.copy()
                        file_df_copy = file_df_copy[["chromosome","chrnum","start","end","comp"]]
                        file_df_copy.comp[file_df.comp == 0.0] = -1.0
                        file_df_copy.comp[file_df.comp == -1.0] = 0.0
                        list_df.append(file_df_copy)
                        
            res_df = pd.concat(list_df)
            res_df = res_df.sort_values(by = ["chrnum","start"])
            filename = os.path.join(path_parent,"SB3_converted_SB1",cell_type + "_" + resolution + "_COMPARTMENT" )
            res_df.to_csv(filename,header = False, index = False, sep = " ")
                        


