#!/usr/bin/env python
# coding: utf-8

# In[1]:


import requests
from bs4 import BeautifulSoup, NavigableString 


# In[2]:


def find_directory(url):
    soup = BeautifulSoup(requests.get(url).text) 
    table_body=soup.find('body')
    rows = table_body.find_all('tr')
    liste_type = []
    liste_ajout_url = []
    a = None
    for row_number, row in enumerate(rows):
        cols=row.find_all('td')
        for number, x in enumerate(cols):
            b = x.find('a')
            if a != None and b != None:
                if a.get('src') != '/icons/back.gif':
                    if a.get('src') == '/icons/folder.gif':
                        liste_type.append('folder')
                    elif a.get('src') == '/icons/unknown.gif':
                        liste_type.append('file')
                    liste_ajout_url.append(b.get('href'))
            a = x.find('img')
    return liste_type, liste_ajout_url


# In[3]:


def terminer(liste):
    liste_repons = ['yes']
    for i in liste:
        if i[-12:] == '.RAWobserved':
            liste_repons.append('yes')
        else:
            liste_repons.append('no')
            break
    if len(set(liste_repons))==1:
        repons = 'stop'
    else:
        repons = 'continuer'
    if len(liste)==1:
        repons = 'continuer'
    return repons


# In[4]:


def approfondissement_url(url):
    liste_url = [url]
    while terminer(liste_url) == 'continuer':
        for element in liste_url:
            if element[-12:] != '.RAWobserved':
                a, b = find_directory(element)
                little_list = []
                for little_element in b:
                    little_list.append(element + little_element)
                liste_url += little_list
                liste_url.remove(element)
    return liste_url


# In[8]:


url = 'http://www.lcqb.upmc.fr/meetu/dataforstudent/HiC/'
a = approfondissement_url(url)


# In[9]:


