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

Additionally, the scripts look for bsaI regions, and regions for other enzymes of interest, in the theoretical product. For detection in the product, a column is 
added to the .csv denoting if bsaI is present or not.

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
python FullDesigner.py -i inputfile.fasta
```
or use:
```bash
python FullDesigner.py -i inputfile.fasta -n 8
```
to change the number of primer pairs generated (n can be any integer, and the default is 5).

Users can also select the length of the primer product that is generated. The default is between 200 and 500 bp.
However, Primer3 will generally create primers to the shorter end of this range.

To choose the lengths users can use the following arguments:
```bash
python FullDesigner.py -i inputfile.fasta -l=500 -u=550
```
Where -l is for the lower length range and -u for the upper. 

### Author
Joshua M Ball (joshua.ball@earlham.ac.uk)
