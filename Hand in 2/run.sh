#!/bin/bash


echo "Clearing/creating the plotting directory"
if [ ! -d "plots" ]; then
  mkdir plots
fi
rm -rf plots/*

# Script that returns a plot
echo "Run the first script ..."
python3 Exercise1abc.py > exercise1a.txt

echo "Run the second script ..."
python3 Exercise1d.py > exercise1d.txt

echo "Run the third script ..."
python3 Exercise2a.py > exercise2a.txt

echo "Run the fourth script ..."
python3 Exercise2b.py > exercise2b.txt

echo "Generating the pdf"

pdflatex genugten.tex
bibtex genugten.aux
pdflatex genugten.tex
pdflatex genugten.tex
