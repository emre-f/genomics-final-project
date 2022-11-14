import urllib3, gzip, shutil, os
import pandas as pd
import numpy as np
import glob, json

'''
Not used because prokka files themselves only had general COG terms which wasn't specific enough.
Therefore further data analysis was done. See read_go_annotations.py
'''

# directory = "./prokka_genomes/"

# genes = {} # key: gene COG name, value: { good_count: int, bad_count: int }
# prokkas = {} # key: prokka name, values: { good:(bool), genes:(array of strings) }

# for subfolder in glob.glob(directory + '*'): # for all subfolders (all resulting prokkas)
#     genome_name = subfolder[subfolder.index('\\') + 1 : subfolder.index('.fna')] # prokka_genomes/x.fna (getting the x)
#     split_name = genome_name.split('_')

#     good = split_name[0] == 'good' # good or bad?
#     bacteria_name = split_name[1]  # the bacteria name

#     prokkas[bacteria_name] = {
#         'good': good,
#         'genes': {} # keep track of counts
#     } 

#     tsvs = glob.glob(subfolder + '/*.tsv')
#     print(tsvs)
#     if len(tsvs) > 1:
#         print("ERROR: Multiple .tsv files found at " + subfolder)
#     elif len(tsvs) == 0:
#         print("ERROR: No .tsv files found at " + subfolder)
#     else: # should be only 1
#         data = pd.read_csv(tsvs[0], sep='\t')
#         cog_numbers = data['COG']
#         for cn in cog_numbers:
#             if not pd.isnull(cn):
                
#                 # add it to the prokkas dict
#                 prokkas[bacteria_name]['genes'][cn] = prokkas[bacteria_name]['genes'].get(cn, 0) + 1
    
#     # print(prokkas)
#     break