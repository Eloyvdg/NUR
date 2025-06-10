#!/bin/bash


echo "Clearing/creating the plotting directory"
if [ ! -d "plots" ]; then
  mkdir plots
fi
rm -rf plots/*

if [ ! -e galaxy_data.txt ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/galaxy_data.txt 
fi

# Script that returns a plot
echo "Run the first script ..."
python3 Ex1.py

echo "Run the second script ..."
python3 Ex2.py

echo "Run the third script ..."
python3 Ex3.py > results_3c.txt

echo "Generating the pdf"

pdflatex genugten.tex
bibtex genugten.aux
pdflatex genugten.tex
pdflatex genugten.tex
