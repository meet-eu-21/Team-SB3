#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 09:58:58 2021

@author: apauron
"""
import wget
from pathlib import Path
import Compartments_SB3 ## To do the pipeline
import get_files ## To get the filename

url = 'http://www.lcqb.upmc.fr/meetu/dataforstudent/HiC/'
a = get_files.approfondissement_url(url)
# In[9]:


for file in a :
    file = file.replace("http://www.lcqb.upmc.fr/meetu/dataforstudent/HiC","")
    print(file)