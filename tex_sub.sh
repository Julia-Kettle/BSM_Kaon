#!/bin/bash

for BASIS in Lattice
do
for SCHEME in MOM
do
for KIN in E
do

cp results_TEMPLATE.tex tmp.tex
sed -e "s@BASIS@${BASIS}@g" -i.bak tmp.tex
sed -e "s@SCHEME@${SCHEME}@g" -i.bak tmp.tex
sed -e "s@KIN@${KIN}@g" -i.bak tmp.tex
mv tmp.tex results_${BASIS}_${SCHEME}_${KIN}.tex
pdflatex results_${BASIS}_${SCHEME}_${KIN}.tex 

done
done
done
