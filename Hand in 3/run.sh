#!/bin/bash


echo "Clearing/creating the plotting directory"
if [ ! -d "plots" ]; then
  mkdir plots
fi
rm -rf plots/*

if [ ! -e satgals_m11.txt ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m11.txt
fi

if [ ! -e satgals_m12.txt ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m12.txt
fi

if [ ! -e satgals_m13.txt ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m13.txt
fi

if [ ! -e satgals_m14.txt ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m14.txt
fi

if [ ! -e satgals_m15.txt ]; then
  wget https://home.strw.leidenuniv.nl/~daalen/Handin_files/satgals_m15.txt
fi

# Script that returns a plot
echo "Run the first script ..."
python3 Ex1a.py > Ex1a.txt

echo "Run the second script ..."
python3 Ex1b.py > Ex1b.txt

echo "Run the third script ..."
python3 Ex1c.py > Ex1c.txt

#echo "Run the fourth script ..."
#python3 Ex1d.py > Ex1d.txt

echo "Generating the pdf"

pdflatex genugten.tex
bibtex genugten.aux
pdflatex genugten.tex
pdflatex genugten.tex
