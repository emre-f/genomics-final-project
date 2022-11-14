#!/bin/bash
# For every folder under /prokka/genomes/, finds the .gff file and creates gene ontology (GO) annotations
# Files are stored under /uniprots/ and /annotations/, both under /go_annotations/

for i in prokka/genomes/*
do
    cat prokka/genomes/"${i##*/}"/PROKKA.gff | grep -o "UniProtKB.*;" | awk -F'[:;=]' '{print $4" "$2}' >go_annotations/uniprots/"${i##*/}".txt
    python2.7 /usr/bin/uniprot2go.py -i go_annotations/uniprots/"${i##*/}.txt" -d /fsys1/data/uniprot2go/uniprot-vs-go-db.sl3 >go_annotations/annotations/"${i##*/}_go.annotations"
done