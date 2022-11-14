import urllib3, gzip, shutil, os
import pandas as pd
import numpy as np
import glob, json
from collections import OrderedDict

target_file = './GO_terms_raw.csv'

def get_go_term(GO_ID):

    # CONSTRUCT WEB LINK
    base = 'ebi.ac.uk/QuickGO/services/ontology/go/terms/'

    http = urllib3.PoolManager()
    r = http.request('GET', base + GO_ID, preload_content=False) # MAKE THE REQUEST
    go_json = {}
    
    if r.status != 200: # failed to get file
        print('Failed to get file for ' + out_file_name + ', status ' + repr(r.status))
    else:
        while True:
            data = r.read()
            if not data:
                break
            go_json = json.loads(data)
    
    r.release_conn()

    go_id = go_json['results'][0]['id']

    if GO_ID != go_id:
        print('ERROR: Got different GO ID. Input: ' + GO_ID + ", Got Back: " + go_id)

    go_name = go_json['results'][0]['name']
    go_description = go_json['results'][0]['definition']['text']

    return go_name, go_description

def get_all_go_terms(target_file, write=False):
    data = pd.read_csv(target_file)

    detailed_go_terms = {}
    TOTAL_SCORE_CUTOFF, UNIQUE_SCORE_CUTOFF = 1, 0.001 # To not get 5000 terms

    for row in data[data.keys()[0]]:
        go_id = data['GO Term'][row]
        total_score = data['Total Score'][row]
        unique_score = data['Unique Score'][row]

        if pd.isnull(go_id):
            continue # no name

        if abs(total_score) <= TOTAL_SCORE_CUTOFF or abs(unique_score) <= UNIQUE_SCORE_CUTOFF:
            continue # not significant
            
        print("Getting term " + go_id)
        name, description = get_go_term(go_id)

        detailed_go_terms[go_id] = {
            'Total Score': total_score,
            'Unique Score': unique_score,
            'Name': name,
            'Description': description
        }
    
    if not write:
        return # no file to write to
    
    df = pd.DataFrame(columns=['GO Term', 'Total Score', 'Unique Score', 'Function Name', 'Function Description'])
    keys = list(detailed_go_terms.keys())
    for i in range(len(keys)):
        info = detailed_go_terms[keys[i]]
        df.loc[i] = [
            keys[i], 
            '{0:.2f}'.format(info['Total Score']), 
            '{0:.4f}'.format(info['Unique Score']),
            info['Name'],
            info['Description']
        ]
    
    df.to_csv(r'./GO_terms_detailed.csv', sep=",", encoding='utf-8', header='true')

    

get_all_go_terms(target_file, True)