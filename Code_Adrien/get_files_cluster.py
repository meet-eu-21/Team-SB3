#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 10:23:17 2021

@author: apauron
"""

from os import walk
import os


def getfiles(mypath) :

    f = []
    for dirpath, dirnames, filenames in walk(mypath, topdown = False) :
            for name in filenames :
                f.append(os.path.join(dirpath,name))
    

    rawobserveds = []
    
    for file in f :
        if ".RAWobserved" in file :
            rawobserveds.append(file)
    
    return rawobserveds
