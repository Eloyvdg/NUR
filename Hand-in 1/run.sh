#!/bin/bash


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

echo "Run the third script ..."
python3 Timeit.py > timeit.txt

echo "Generating the pdf"

pdflatex genugten.tex
bibtex genugten.aux
pdflatex genugten.tex
pdflatex genugten.tex


