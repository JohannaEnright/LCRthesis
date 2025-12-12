#!/bin/bash

#takes filename as an argument
#will print out the number of pattern types observed in 
# say for example ~/proj1/periodicity_info/dna_repeats/S.cerevisiae/wholeseq/periodicity_dnawholeseq_1154222222_Saccharomyces_cerevisiae451.31.5_dna_proEntOutput


echo 1:; more $1 | grep -v "#" | cut -f 2 | grep -c 1
echo 2:; more $1 | grep -v "#" | cut -f 2 | grep -c 2
echo 3:; more $1 | grep -v "#" | cut -f 2 | grep -c 3
echo 4:; more $1 | grep -v "#" | cut -f 2 | grep -c 4
echo 5:; more $1 | grep -v "#" | cut -f 2 | grep -c 5
echo 6:; more $1 | grep -v "#" | cut -f 2 | grep -c 6
echo 7:; more $1 | grep -v "#" | cut -f 2 | grep -c 7
echo 8:; more $1 | grep -v "#" | cut -f 2 | grep -c 8
echo 9:; more $1 | grep -v "#" | cut -f 2 | grep -c 9
