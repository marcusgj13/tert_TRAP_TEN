# tert TRAP/TEN analysis

This repository contains python scripts used to isolate sequences of the TEN and TRAP domains from  
full-length telomerase reverse transcriptase (TERT) sequences. All scripts were intially implemented  
in python 3.5 but should be backwards compatible with python2.7.

Author: Marcus Gallagher-Jones
Institution: UCLA Department of chemistry and biochemistry
Email: marcusgj13@gmail.com

# Explanation  
  
The scripts directory contains the python scripts necessary to isolate the TRAP and TEN domains  
from full-length sequences in the SEQUENCES directory. Additional files are for the analysis of  
any secondary structure prediction or disorder prediction performed on isolated sequences.  

# Running scripts

findTEN.py and findTRAP.py require no inputs. Other scripts take either a fasta file as an input, 
a clustal results file or a text file that has a list of output files to be analysed (.horiz files  
for psipred and .diso files for disopred).
