import urllib3, gzip, shutil, os
import pandas as pd
import numpy as np
import glob, json
import copy
from collections import OrderedDict

# For visualization
import seaborn as sns
import matplotlib.pyplot as plt
# plt.rcParams['font.size'] = '7'

def split_dict_to_multiple(input_dict, max_limit=200):
    """Splits dict into multiple dicts with given maximum size. 
    Returns a list of dictionaries."""
    chunks = []
    curr_dict ={}
    for k, v in input_dict.items():
        if len(curr_dict.keys()) < max_limit:
            curr_dict.update({k: v})
        else:
            chunks.append(copy.deepcopy(curr_dict))
            curr_dict = {k: v}
    # update last curr_dict
    chunks.append(curr_dict)
    return chunks

def take_top_GO_terms(count):
    target_file = "./GO_terms_raw.csv"
    data = pd.read_csv(target_file) 
    print("[0] GETTING GO TERM DATA FROM FILE")
    # GET THE ORDERED TERMS
    terms = {}
    for row in data[data.keys()[0]]:
        if pd.isnull(data['GO Term'][row]):
            print("ERROR: Null Go term")
            continue 

        terms[data['GO Term'][row]] = {
            'Total Score': data['Total Score'][row],
            'Unique Score': data['Unique Score'][row]
        }
    print("[1] COMPLETE")
    ordered_by_total = OrderedDict(sorted(terms.items(), key=lambda i: i[1]['Total Score'], reverse=True))
    ordered_by_total_keys = list(ordered_by_total.keys())

    ordered_by_unique = OrderedDict(sorted(terms.items(), key=lambda i: i[1]['Unique Score'], reverse=True))
    ordered_by_unique_keys = list(ordered_by_unique.keys())
    
    # FROM ORDERED TERMS GET THE TOP 'count' OF TERMS
    top_total_good_terms = []
    top_total_bad_terms = []

    top_unique_good_terms = []
    top_unique_bad_terms = []

    for i in range(count):
        top_total_good_terms.append(ordered_by_total_keys[i])
        top_total_bad_terms.append(ordered_by_total_keys[-(i + 1)])

        top_unique_good_terms.append(ordered_by_unique_keys[i])
        top_unique_bad_terms.append(ordered_by_unique_keys[-(i + 1)])
    print("[2] SORTED & SAVED DATA")
    return top_total_good_terms, top_total_bad_terms, top_unique_good_terms, top_unique_bad_terms
        
def get_all_bacteria_GO_counts():
    good_bacteria = {}
    bad_bacteria = {}

    # FOR EACH BACTERIA, AND EACH GO TERM LIST FILL OUT GO TERM COUNT
    directory = './go_annotations/sorted_proteomes/'

    for bact_file in glob.glob(directory + '*'):
        genome_name = bact_file[bact_file.index('\\') + 1: bact_file.index('.fna')] 
        split_name = genome_name.split('_')
        
        genome_name = split_name[1]
        short_name = 'P. ' + genome_name.split(' ', 1)[1] # P. instead of Pseudomonas
        print("GETTING GO TERM FREQ. FOR " + genome_name)
        if os.stat(bact_file).st_size == 0:
            print("ERROR: Empty file, skipping")
            continue 
        good_or_bad = split_name[0] == 'good' # good or bad?

        data = pd.read_csv(bact_file, sep='\t', header=None) 
        
        if good_or_bad:
            good_bacteria[short_name] = {} 

            for row in data[0]:
                count = int(row[ : row.rfind(' ')]) # count of GO term
                name = row[row.rfind(' ') + 1: ] # name of GO term
                good_bacteria[short_name][name] = count
        else:
            bad_bacteria[short_name] = {}

            for row in data[0]:
                count = int(row[ : row.rfind(' ')]) # count of GO term
                name = row[row.rfind(' ') + 1: ] # name of GO term
                bad_bacteria[short_name][name] = count

    return good_bacteria, bad_bacteria

def makeHeatmapFromCounts(term_lst, bact_lst, max_limit, color, file_name):
    # MIX & MATCH TO MAKE HEAT MAPS
    columns = list(term_lst)
    df = pd.DataFrame(columns=columns) # Will add GO terms as we go along
    bact_name_keys = list(bact_lst)
    for i in range(len(bact_name_keys)):
        name = bact_name_keys[i]
        bact_go_terms = bact_lst[name]
        counts = []

        for go_term in term_lst:
            if pd.isnull(go_term):
                print("ERROR: Null term")
                continue
            counts.append(bact_go_terms.get(go_term, 0))

        df.loc[i] = counts
    
    # PLOT DATA

    # TO SET SPECIFIC COLOR PATTERN
    # rdgn = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=50, sep=40, as_cmap=True)
    # plt.figure(figsize=(len(df.columns), len(df)))
    plt.rcParams.update({'axes.titlesize': 'large'})
    plt.figure(figsize=(200, 100))

    g = sns.heatmap(df, cmap=color, yticklabels=list(bact_lst.keys()), vmin = 0, vmax = max_limit)
    g.set_xticklabels(g.get_xticklabels(), rotation=90)
    plt.savefig('./figures/' + file_name, bbox_inches='tight')

    return df

# # TAKE THE TOP TERMS (Y-AXIS) (fun fact: 197 unique terms out of 200)
top_total_good_terms, top_total_bad_terms, top_unique_good_terms, top_unique_bad_terms = take_top_GO_terms(100)

# # TAKE ALL BACTERIAS AND THEIR GO TERM COUNTS
good_bacteria_with_term_counts, bad_bacteria_with_term_counts = get_all_bacteria_GO_counts()

''' NOTE: Running all at once usually hits memory errors, try to create the figures in chunks '''

''' 0. VISUALS BEGIN HERE '''
plt.rcParams['font.size'] = '24'
''' 1. ALL BACT TOGETHER. Outcome: total counts are useless'''
ALL_BACT_COUNTS = good_bacteria_with_term_counts.copy()
ALL_BACT_COUNTS.update(bad_bacteria_with_term_counts) # COMBINE LISTS
makeHeatmapFromCounts(top_total_good_terms, ALL_BACT_COUNTS, 150, 'Greens', 'all_bact_tot_good')
makeHeatmapFromCounts(top_total_bad_terms, ALL_BACT_COUNTS, 30, 'Reds', 'all_bact_tot_bad')
makeHeatmapFromCounts(top_unique_good_terms, ALL_BACT_COUNTS, 4, 'Greens', 'all_bact_uniq_good')
makeHeatmapFromCounts(top_unique_bad_terms, ALL_BACT_COUNTS, 4, 'Reds', 'all_bact_uniq_bad')

''' 2. UNIQUE COUNTS '''
plt.rcParams['font.size'] = '48'
''' 2.1. GOOD BACT '''
makeHeatmapFromCounts(top_unique_good_terms, good_bacteria_with_term_counts, 4, 'Greens', 'good_bact_uniq_good')
makeHeatmapFromCounts(top_unique_bad_terms, good_bacteria_with_term_counts, 4, 'Reds', 'good_bact_uniq_bad')

''' 2.2. BAD BACT '''
split_bad_dicts = split_dict_to_multiple(bad_bacteria_with_term_counts, 38)
for i in range(len(split_bad_dicts)):
    makeHeatmapFromCounts(top_unique_good_terms, split_bad_dicts[i], 4, 'Greens', 'bad_bact_uniq_good_' + repr(i))

for i in range(len(split_bad_dicts)):
    makeHeatmapFromCounts(top_unique_bad_terms, split_bad_dicts[i], 4, 'Reds', 'bad_bact_uniq_bad_' + repr(i))

''' 3. ALL GO TERMS '''

''' 3.1 GOOD BACT '''
makeHeatmapFromCounts(top_unique_good_terms + top_unique_bad_terms, good_bacteria_with_term_counts, 4, 'Greens', 'good_bact_uniq_all')

''' 3.2. BAD BACT '''
split_bad_dicts = split_dict_to_multiple(bad_bacteria_with_term_counts, 38)
for i in range(len(split_bad_dicts)):
    makeHeatmapFromCounts(top_unique_good_terms + top_unique_bad_terms, split_bad_dicts[i], 4, 'Reds', 'bad_bact_uniq_all_' + repr(i))


''' OLD STUFF [IGNORE] '''
# Note: None of them occur more than 10 times max (even ~5-6 in other lists)
# 
# makeHeatmapFromCounts(top_unique_good_terms, bad_bacteria_with_term_counts, 4, 'Greens')
# 
# makeHeatmapFromCounts(top_unique_bad_terms, bad_bacteria_with_term_counts, 4, 'Reds')

# TOTAL TERMS (not very descriptive)
# makeHeatmapFromCounts(top_total_good_terms, good_bacteria_with_term_counts, 150, 'Greens')
# makeHeatmapFromCounts(top_total_good_terms, bad_bacteria_with_term_counts, 400, 'Greens')
# makeHeatmapFromCounts(top_total_bad_terms, good_bacteria_with_term_counts, 100, 'Reds')
# makeHeatmapFromCounts(top_total_bad_terms, bad_bacteria_with_term_counts, 400, 'Reds')



