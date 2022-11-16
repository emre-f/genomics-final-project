import urllib3, gzip, shutil, os
import pandas as pd
import numpy as np
import glob, json
from collections import OrderedDict

# For visualization
import seaborn as sns
import matplotlib.pyplot as plt
# plt.rcParams['font.size'] = '7'

def take_top_GO_terms(count):
    target_file = "./GO_terms_raw.csv"
    data = pd.read_csv(target_file) 
    print("[0] GETTING GO TERM DATA FROM FILE")
    # GET THE ORDERED TERMS
    terms = {}
    for row in data[data.keys()[0]]:
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

def makeHeatmapFromCounts(term_lst, bact_lst, max_limit, color):
    # MIX & MATCH TO MAKE HEAT MAPS
    columns = list(term_lst)
    
    df = pd.DataFrame(columns=columns) # Will add GO terms as we go along
    bact_name_keys = list(bact_lst)
    for i in range(len(bact_name_keys)):
        name = bact_name_keys[i]
        bact_go_terms = bact_lst[name]
        counts = []

        for go_term in term_lst:
            counts.append(bact_go_terms.get(go_term, 0))

        df.loc[i] = counts
    
    # PLOT DATA

    # TO SET SPECIFIC COLOR PATTERN
    # rdgn = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=50, sep=40, as_cmap=True)
    plt.figure(figsize=(len(df), len(df.columns)))
    g = sns.heatmap(df, cmap=color, yticklabels=list(bact_lst.keys()), vmin = 0, vmax = max_limit)
    g.set_xticklabels(g.get_xticklabels(), rotation=90)
    plt.show()

    return df

# # TAKE THE TOP TERMS (Y-AXIS) (fun fact: 197 unique terms out of 200)
top_total_good_terms, top_total_bad_terms, top_unique_good_terms, top_unique_bad_terms = take_top_GO_terms(80)

# # TAKE ALL BACTERIAS AND THEIR GO TERM COUNTS
good_bacteria_with_term_counts, bad_bacteria_with_term_counts = get_all_bacteria_GO_counts()

# Note: None of them occur more than 10 times max (even ~5-6 in other lists)
makeHeatmapFromCounts(top_unique_good_terms, good_bacteria_with_term_counts, 4, 'Greens')
makeHeatmapFromCounts(top_unique_good_terms, bad_bacteria_with_term_counts, 4, 'Greens')
makeHeatmapFromCounts(top_unique_bad_terms, good_bacteria_with_term_counts, 4, 'Reds')
makeHeatmapFromCounts(top_unique_bad_terms, bad_bacteria_with_term_counts, 4, 'Reds')

makeHeatmapFromCounts(top_total_good_terms, good_bacteria_with_term_counts, 150, 'Greens')
makeHeatmapFromCounts(top_total_good_terms, bad_bacteria_with_term_counts, 400, 'Greens')
makeHeatmapFromCounts(top_total_bad_terms, good_bacteria_with_term_counts, 100, 'Reds')
makeHeatmapFromCounts(top_total_bad_terms, bad_bacteria_with_term_counts, 400, 'Reds')

# ALL OF THEM
good_bacteria_with_term_counts.update(bad_bacteria_with_term_counts) # COMBINE LISTS
makeHeatmapFromCounts(top_unique_good_terms, good_bacteria_with_term_counts, 4, 'Greens')
makeHeatmapFromCounts(top_unique_bad_terms, good_bacteria_with_term_counts, 4, 'Reds')
makeHeatmapFromCounts(top_total_good_terms, good_bacteria_with_term_counts, 150, 'Greens')
makeHeatmapFromCounts(top_total_bad_terms, good_bacteria_with_term_counts, 30, 'Reds')