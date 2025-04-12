#!/bin/bash


echo "Clearing/creating the plotting directory"
if [ ! -d "plots" ]; then
  mkdir plots
fi
rm -rf plots/*

# Script that returns a plot
echo "Run the first script ..."
python3 Ex1a.py > Ex1a.txt

echo "Run the second script ..."
python3 Ex1b.py > Ex1b.txt

echo "Run the third script ..."
python3 Ex1c.py > Ex1c.txt

echo "Run the fourth script ..."
python3 Ex1d.py > Ex1d.txt

echo "Generating the pdf"

pdflatex genugten.tex
bibtex genugten.aux
pdflatex genugten.tex
pdflatex genugten.tex
