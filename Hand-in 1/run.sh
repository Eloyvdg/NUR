#!/bin/bash

echo "Run handin template -- this is an example, make sure to edit this (e.g. rename template.tex to [YourLastName].tex)!"

echo "Clearing/creating the plotting directory"
if [ ! -d "plots" ]; then
  mkdir plots
fi
rm -rf plots/*

# Script that returns a plot
echo "Run the first script ..."
python3 Poisson.py > poisson.txt

echo "Run the second script ..."
python3 Vandermonde.py > vandermonde_output.txt

echo "Generating the pdf"

pdflatex template.tex
bibtex template.aux
pdflatex template.tex
pdflatex template.tex


