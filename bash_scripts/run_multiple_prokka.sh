#!/bin/bash
# For every genome under /genomes/ runs prokka

for i in genomes/*
do
if [ -f prokka/genomes/"${i%}" ]; then
    echo "${i%} exists already."
else
    prokka "${i%}" --outdir prokka/"${i%}"
fi
done