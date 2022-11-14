import urllib3, gzip, shutil, os
import pandas as pd
import numpy as np

# Helper: Check if function is gzipped
def is_gzipped(filepath):
    if not os.path.isfile(filepath):
        return

    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

# Download file from url into out_dir directory with out_file_name name
def download_file(url, out_dir, out_file_name):
    http = urllib3.PoolManager()
    r = http.request('GET', url, preload_content=False) # MAKE THE REQUEST
    path = out_dir + out_file_name
    
    if r.status != 200: # failed to get file
        print('Failed to get file for ' + out_file_name + ', status ' + repr(r.status))
    else:
        with open(path, 'wb') as out:
            while True:
                data = r.read()
                if not data:
                    break
                out.write(data)

    r.release_conn()

# Unzip a gzip file if it is gzipped
def unzip_gzip(path, delete_old_file=True):
    if is_gzipped(path):
        unzipped_path = path[:path.index(".gz")]

        with gzip.open(path, 'r') as f_in, open(unzipped_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        
        if delete_old_file: 
            os.remove(path) # delete old file

        return unzipped_path
    
    return None

# Get genome from NCBI
def get_genome(ncbi_refseq_id, ncbi_gbk_id, custom_name=""):

    # CONSTRUCT WEB LINK
    base = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

    three_letter_code = ncbi_refseq_id[:ncbi_refseq_id.index('_')]
    just_refseq_id = ncbi_refseq_id[ncbi_refseq_id.index('_') + 1: ncbi_refseq_id.index('.')] # 006715895

    directory = three_letter_code + "/" + \
                just_refseq_id[0:3] + "/" + just_refseq_id[3:6] + "/" + just_refseq_id[6:9] + "/" + ncbi_refseq_id + "_" + ncbi_gbk_id
    
    # DOWNLOADING THE FILE
    file_name = ncbi_refseq_id + "_" + ncbi_gbk_id + "_genomic.fna.gz"
    url = base + "/" + directory + "/" + file_name
    
    # CLEAR URL FROM SPACES
    url = url.replace(' ', '_')

    download_file(url, "./downloads/", custom_name + ".fna.gz" if custom_name else file_name)

    # UNZIP FILE
    path = "./downloads/" + custom_name + ".fna.gz" if custom_name else file_name
    unzip_gzip(path)

# Gets the genomes for all genomes in the database
def download_all_genomes(limit=float('inf')):
    df = pd.read_excel('./469_final_db.xlsx') # a DataFrame

    counter = 0
    for index, row in df.iterrows():
        if counter >= limit:
            break

        if row['GP'] or row['BP']: # a good pseudomonas
            ncbi_refseq_id = row['NCBI REFSEQ ID']
            ncbi_assembly_id = row['NCBI ASSEMBLY ID']
            name = str(row['STRAIN'])
            if '[' in name: # remove the 'identifiers' such as [PPATH]
                name = name[:name.index('[') - 1] 

            # Check if data exists
            if ncbi_refseq_id == 'not_announced' or pd.isnull(ncbi_refseq_id) or \
            ncbi_assembly_id == 'not_announced' or pd.isnull(ncbi_assembly_id):
                continue

            # Good or bad bacteria?
            tag = "good_" if row['GP'] else "bad_"

            # If file already downloaded skip
            if os.path.isfile('./downloads/' + tag + name + '.fna'):
                continue

            get_genome(ncbi_refseq_id, ncbi_assembly_id, tag + name)
            counter += 1

# download_all_genomes()