#!/bin/bash
# Finds every .gff file recursively, and names them PROKKA.gff

shopt -s globstar nullglob dotglob
for file in **/*.gff; do
    mv "$(dirname "${file}")/$(basename "${file%.*}").gff" "$(dirname "${file}")/PROKKA.gff"
done