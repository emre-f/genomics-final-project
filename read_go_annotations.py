import urllib3, gzip, shutil, os
import pandas as pd
import numpy as np
import glob, json
from collections import OrderedDict

directory = './go_annotations/sorted_proteomes/'

'''
COLLECT GO TERMS INTO A DICT

key: go term name
value: {
        key: good
        value: {
            unique_count
            total_count
        },
        key: bad
        value: {
            unique_count
            total_count
        }
    }
'''
go_terms = {}
GOOD_COUNT, BAD_COUNT = 0, 0
for go_file in glob.glob(directory + '*'):
    genome_name = go_file[go_file.index('\\') + 1: go_file.index('.fna')] # prokka_genomes/{name}.fna
    split_name = genome_name.split('_')

    good_or_bad = split_name[0] # good or bad?
    if good_or_bad == 'good':
        GOOD_COUNT += 1
    else:
        BAD_COUNT += 1

    bacteria_name = split_name[1]  # the bacteria name

    if os.stat(go_file).st_size == 0:
        print("ERROR: Empty file at " + go_file)
        continue 

    # Read the table file
    data = pd.read_csv(go_file, sep='\t', header=None) 
    
    for row in data[0]:
        count = int(row[ : row.rfind(' ')]) # count of GO term
        name = row[row.rfind(' ') + 1: ] # name of GO term

        if name not in go_terms:
            go_terms[name] = {
                'good': { 
                    'unique_count': 0,
                    'total_count': 0
                    },
                'bad': {
                    'unique_count': 0,
                    'total_count': 0
                    }
                }

        go_terms[name][good_or_bad]['unique_count'] += 1
        go_terms[name][good_or_bad]['total_count'] += count
    
# MAKE TOP-LISTS (top good, top unique good, top bad, top unique bad) #
# ordered_by_total_good = OrderedDict(sorted(go_terms.items(), key=lambda i: i[1]['good']['total_count'], reverse=True))
# ordered_by_total_bad = OrderedDict(sorted(go_terms.items(), key=lambda i: i[1]['bad']['total_count'], reverse=True))
# ordered_by_unique_good = OrderedDict(sorted(go_terms.items(), key=lambda i: i[1]['good']['unique_count'], reverse=True))
# ordered_by_unique_bad = OrderedDict(sorted(go_terms.items(), key=lambda i: i[1]['bad']['unique_count'], reverse=True))

# print(list(ordered_by_total_good.keys())[:5])
# print(list(ordered_by_total_bad.keys())[:5])
# print(list(ordered_by_unique_good.keys())[:5])
# print(list(ordered_by_unique_bad.keys())[:5])

'''
For every GO term in the go_terms dictionary, compare their 'good' and 'bad' counts

Because there are ~200 bad and ~40 good genomes, the counts have been "normalized" by dividing them
by the number of the total count. (eg. every bad genome is divided by ~200, and good ones by ~40)
'''
diff_total = {}
diff_unique = {}

GOOD_COUNT, BAD_COUNT = float(GOOD_COUNT), float(BAD_COUNT)

for gt in go_terms.keys():
    diff_total[gt] = go_terms[gt]['good']['total_count'] / GOOD_COUNT - go_terms[gt]['bad']['total_count'] / BAD_COUNT
    diff_unique[gt] = go_terms[gt]['good']['unique_count'] / GOOD_COUNT - go_terms[gt]['bad']['unique_count'] / BAD_COUNT

ordered_by_diff_total = OrderedDict(sorted(diff_total.items(), key=lambda i: i[1], reverse=True))
ordered_by_diff_unique = OrderedDict(sorted(diff_unique.items(), key=lambda i: i[1], reverse=True))

tot_keys_lst = list(ordered_by_diff_total.keys())
unique_keys_lst = list(ordered_by_diff_unique.keys())

def write_raw_data_to_table():
    df = pd.DataFrame(columns=['GO Term', 'Total Score', 'Unique Score'])
    keys = list(go_terms.keys())
    for i in range(len(keys)):
        df.loc[i] = [keys[i], '{0:.2f}'.format(diff_total[keys[i]]), '{0:.4f}'.format(diff_unique[keys[i]])]
    
    df.to_csv(r'./GO_terms_raw.csv', sep=",", encoding='utf-8', header='true')

def write_sorted_data_to_table():
    df = pd.DataFrame(columns=['TOT GOOD', 'TOT BAD', 'UNIQ GOOD', 'UNIQ BAD'])
    for i in range(len(go_terms)):
        
        # Get ith total occured GO name and its count
        total_top_name = tot_keys_lst[i]
        total_top_count = ordered_by_diff_total[tot_keys_lst[i]]

        # Get bottom ith total occured GO name and its count
        total_bottom_name = tot_keys_lst[-(i + 1)]
        total_bottom_count = ordered_by_diff_total[tot_keys_lst[-(i + 1)]]

        # Get ith total uniquely occured GO name and its count
        unique_top_name = unique_keys_lst[i]
        unique_top_count = ordered_by_diff_unique[unique_keys_lst[i]]

        # Get bottom ith total uniquely occured GO name and its count
        unique_bottom_name = unique_keys_lst[-(i + 1)]
        unique_bottom_count = ordered_by_diff_unique[unique_keys_lst[-(i + 1)]]

        # Save them into a final .csv file 
        df.loc[i] = [
            total_top_name + ' > {0:.2f}'.format(total_top_count),
            total_bottom_name +  ' > {0:.2f}'.format(abs(total_bottom_count)),
            unique_top_name +  ' > {0:.2f}'.format(unique_top_count),
            unique_bottom_name +  ' > {0:.2f}'.format(abs(unique_bottom_count))
        ]

    df.to_csv(r'./GO_term_frequencies.csv', sep=",", encoding='utf-8', header='true')




