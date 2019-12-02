#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 09:38:57 2019

@author: ip10
"""
'''
A python script to take four files from galaxy, combine flanks from either side of an insertion to output a bedfile covering each region 
whilst esuring this bedfile is non-redundant

All four galaxy fasta files are required to be in the working directory

'''

from Bio import SeqIO
import re

# a function to find each flanks partner by taking the files and the location of the insertions
# input two files with flanking regions outputs

def finding_partner(file1, file2):
    G31 = []
    G32 = []
 
    for record in SeqIO.parse(file1, "fasta"):
        a = re.search(r"chr([\dXYUn]{,2}_\w{,10}_[a-z]+)_", record.description) 

        if a is None:
            d = re.search(r"chr[\dXY]{1,2}", record.description)
            if d is not None:
                G31.append(record.description)

       
    for record in SeqIO.parse(file2, "fasta"):
        a = re.search(r"chr([\dXYUn]{,2}_\w{,10}_[a-z]+)_", record.description) 
        if a is None:
            d = re.search(r"chr[\dXY]{1,2}", record.description)
            if d is not None:
                G32.append(record.description)

    return G31, G32

''' 
a function to write the non-redundant list to a bed file.
This function can be adapted whether you want the format of the bedfile to include 'chr' in the chromosome column 
by altering the brackets in the regex.
'''
    
def whatis(G,F):    
    count = 0
    with open('list_notchr.bed', 'a') as output:
        for gy in G:
            
            e = re.search(r"chr([\dXY]{,2})_(\d+)_(\d+)_([+-])\s(\w+)", gy)
            f = re.search(r"chr([\dXY]{,2})_(\d+)_(\d+)_([+-])\s(\w+)", F[count])
            yes = 0
            while yes == 0:
                if e.group(1) == f.group(1):
                    if int(e.group(2)) < int(f.group(2)):
                        if (int(f.group(2)) + 10000) > int(e.group(2)):
    
                            output.writelines(str(e.group(1))+"\t"+str(e.group(2))+"\t"+str(f.group(3))+"\t"+str(e.group(4)) +"\t"+str(e.group(5))+"\n")
                            count +=1
                            yes +=1
                    elif (int(e.group(2)) + 10000) > int(f.group(2)):
                        output.writelines(str(e.group(1))+"\t"+str(f.group(2))+"\t"+str(e.group(3))+"\t"+str(e.group(4)) +"\t"+str(e.group(5))+"\n")
                        count +=1
                        yes +=1
                else:
                    count +=1


hiy, jiy = finding_partner("/Users/ip10/Documents/Data/Galaxy34.fasta", "/Users/ip10/Documents/Data/Galaxy33.fasta")
ghy, jhu = finding_partner("/Users/ip10/Documents/Data/Galaxy31.fasta","/Users/ip10/Documents/Data/Galaxy32.fasta")

whatis(hiy,jiy)
whatis(ghy,jhu)

# a check to make sure start and finish locations are the correct way round in the bedfile created
with open('list_nonchr.bed') as check:
    umm = check.readlines()
    for lin in umm:
        e = re.search(r"([\dXY]{,2})\s(\d+)\s(\d+)\s([+-])\s(\w+)", lin)
        diff = int(e.group(3)) - int(e.group(2))
        print(diff)
        
