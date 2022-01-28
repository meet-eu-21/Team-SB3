#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 10:23:17 2021

@author: apauron
"""

from os import walk
import os

############################### Get the list of all the files path in a folder, containing an expression such as .RAWobserved #############################"

def getfiles(mypath,expression) :
    
    
    """
    A small method returning all the files you want in a folder
    
    Keyword arguments :
    mypath -- the path containing the folder in which there are the files you want to list
    expression -- Gets only the files containing this expression, ex '.RAWobserved'


    
    Returns :
    
    a list rawobserveds of the files full paths
    
    
    """

    ## An empty list
    f = []
    for dirpath, dirnames, filenames in walk(mypath, topdown = False) : #walk allows to penetrate recursively folders until we reach files
            for name in filenames :
                f.append(os.path.join(dirpath,name)) #We get the full path of the file
    

    rawobserveds = []
    
    #Filter the files according to the expression
    
    for file in f :
        if expression in file :
            rawobserveds.append(file)
    
    return rawobserveds
