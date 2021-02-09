# Primer Designer

Primer Designer is a python script for designing primers for genes using Primer3.

The full sequence version (FullDesigner.py) takes a FASTA file of genes and flanking regions, formatted as:
```text
>Gene_name
atcgGATCtgac
```
Where flanking regions are in lower case and the gene/ coding sequence is Upper case.

The flanking sequence version (FlankingDesigner.py) takes a FASTA file of flanking regions with a segment of gene,
formatted as:
```text
>Gene_name_LF
atcgGA

>Gene_name_RF
TCtgac
```
Where sequences are split into LF (left flanking region) and RF (right flanking region) and include a small portion
of the coding sequence (in test cases we used 200bp).

The output for both method is  a .csv file of left and right primers for each gene, along with additional information.

Additionally the scripts look for bsaI regions in both the gene/ coding sequence, and the primers. If there is one of 
these in the gene; the script will not generate primer for the specific gene, and will instead add the gene name to a
.txt file. For detection in primers, a column is added to the .csv denoting if bsaI is present in the primer or not.

(N.B. For the FlankingDesigner script this will only be looking in the included section of coding sequence, so if there 
is bsaI further in it will not be detected.)

## Installation

Make sure that you have biopython and primer3-py installed first
```bash
pip install biopython
pip install primer3-py
```

Download the FullDesigner.py, FlankingDesigner.py and PrimerFilters.py files and keep them in the same directory.

## Usage

To run the script from the command line, with an input FASTA file:
```bash
python main.py -i inputfile.fasta
```
or use:
```bash
python main.py -i inputfile.fasta -n 8
```
to change the number of primer pairs generated (n can be any integer, and the default is 5).

### Author
Joshua M Ball (joshua.ball@earlham.ac.uk)