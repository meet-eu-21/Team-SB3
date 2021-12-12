#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 12 10:23:17 2021

@author: apauron
"""

from os import listdir
from os.path import isfile, join


def getfiles(mypath) :
    return [f for f in listdir(mypath) if isfile(join(mypath, f))]

