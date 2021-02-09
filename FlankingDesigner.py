__author__ = 'Joshua Ball (joshua.ball@earlham.ac.uk)'
__version__ = '1.0.0'

import argparse
import primer3
from Bio import SeqIO
import csv
from datetime import datetime
from PrimerFilters import *

# argument parser for command line
parser = argparse.ArgumentParser(description='Primer Designer')
parser.add_argument('-i', '--input', help='Input file name', required=True)
parser.add_argument('-n', '--number', help='Number of Primer pairs to return', default=5)
args = parser.parse_args()

# create the .csv file and enter headers. Also assign date and time stamped names for CSV and TXT files to be made.
csvfilename = 'LF-RF Primers '+datetime.today().strftime('%y-%m-%d %H.%M.%S')+'.csv'
txtfilename = 'LF-RF CDS with bsaI '+datetime.today().strftime('%y-%m-%d %H.%M.%S')+'.txt'
headings = ['Gene', 'Primer', 'bsaI in Primer', 'Flank', 'Pair Penalty', 'Left Penalty', 'Right Penalty', 'Primer Forward', 'Primer Reverse',
        'Left (Start, Length)', 'Right (Start, Length)', 'Left TM', 'Right TM', 'Left GC%', 'Right GC%', 'Left Self Any TH',
        'Right Self Any TH', 'Left Self End TH', 'Right Self End TH', 'Left Hairpin TH', 'Right Hairpin TH',
        'Left End Stability', 'Right End Stability', 'Pair Compl Any TH', 'Pair Compl End TH', 'Pair Product Size']
with open(csvfilename, 'w') as f:
    writer = csv.DictWriter(f, fieldnames=headings)
    writer.writeheader()
    f.close()

# using BioPython to get info from fasta file.
# 'Whole Chromosome CDS list 1000bp.fasta' 'Test_Genes.fasta'
myfast = SeqIO.parse(str(args.input), 'fasta') #str(args.input)
for seq_record in myfast:
    ID = seq_record.id
    SEQ = str(seq_record.seq)
    CDS_len = sum(1 for c in SEQ if c.isupper())
    CDS_start = SEQ.find(find_CDS_start(SEQ))
    #print(ID)
    #print(SEQ)
    #print(find_CDS_start(SEQ))
    #print("CDS Start = ", CDS_start)
    #print("CDS Length = ", CDS_len)

    if bsa_in_CDS(SEQ):
        bsaincds = []
        bsaincds.append(ID)
        with open(txtfilename, 'a') as f2:
            for item in bsaincds:
                f2.write('%s\n' % item)
                f2.close()

    else:
        # Use Primer 3 and get results for Left Flanking Region.
        num = int(args.number) # number of primer pairs to generate.

        primerlist = (primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': ID,
                'SEQUENCE_TEMPLATE': SEQ,
            },
            {
                'PRIMER_TASK': 'generic',
                'PRIMER_PICK_LEFT_PRIMER': 1,
                'PRIMER_PICK_RIGHT_PRIMER': 1,
                'PRIMER_NUM_RETURN': num,
                'PRIMER_OPT_SIZE': 20,
                'PRIMER_MIN_SIZE': 18,
                'PRIMER_MAX_SIZE': 30,
                'PRIMER_PRODUCT_SIZE_RANGE': [200, 500],
                'PRIMER_OPT_TM': 55,
                'PRIMER_MIN_TM': 50,
                'PRIMER_MAX_TM': 60,
                'PRIMER_EXPLAIN_FLAG': 1,
                'PRIMER_MAX_END_STABILITY': 6.0,
                'PRIMER_MIN_GC': 45.0,
                'PRIMER_OPT_GC_PERCENT': 50.0,
                'PRIMER_MAX_GC': 80.0
            }
        ))
        #print(primerlist)
        # primerlist is the initial output of primer3

        newlist = cleanupdata(primerlist)
        #print(newlist)
        # newlist is a filtered version of primerlist, number of primers returned and explains have been removed

        # Separate the data based on the primer number in each key! This results in a dictionary of sub-dictionaries where a
        # sub-dictionary is a primer pair
        reg = re.compile(r'\d+') # regex that finds the number (primer number) in the keys (e.g. PRIMER_PAIR_0_PENALTY)
        splitprimers = {} # This will be the dictionary that all sub-dictionaries go into
        for k, v in newlist.items():
            m=reg.search(k)
            if m:
                numb = m.group()
                splitprimers.setdefault(numb, {})[k]=v
            else:
                print('Invalid key:', k, v)

        #print(splitprimers)

        # This adds 2 new keys and values to each primer-pair dictionary. 1 for the Primer pair number, and 1 for the Gene ID
        # This also then turns the dictionary of sub-dictionaries into a LIST of dictionaries called primers, this is important
        # for further steps.

        primers = []
        for j in splitprimers:
            splitprimers[j]['Primer'] = str(int(j)+1)
            splitprimers[j]['Gene'] = ID
            splitprimers[j]['Flank'] = 'Left'
            primers.append(splitprimers[j])

        #print(primers)

        # Check for bsaI and filter primers to remove any that contain bsaI

        fprimers = []
        check = re.compile(r'PRIMER_LEFT_\d+$')
        check2 = re.compile(r'PRIMER_RIGHT_\d+$')
        for i in primers:
            for j in i:
                if check.search(j):
                    primerstart = i[j][0]
                elif check2.search(j):
                    primerend = i[j][0]+1
            primer = SEQ[primerstart:primerend]
            if bsa_in_primer(primer):
                i['bsaI in Primer'] = 'Yes'
                fprimers.append(i)
            else:
                i['bsaI in Primer'] = 'No'
                fprimers.append(i)

        #print(fprimers)

        fprimers = replacekey(fprimers)
        #print(fprimers)

        # Output .txt or .csv
        headings = ['Gene', 'Primer', 'bsaI in Primer', 'Flank', 'Pair Penalty', 'Left Penalty', 'Right Penalty', 'Primer Forward',
                    'Primer Reverse',
                    'Left (Start, Length)', 'Right (Start, Length)', 'Left TM', 'Right TM', 'Left GC%', 'Right GC%',
                    'Left Self Any TH',
                    'Right Self Any TH', 'Left Self End TH', 'Right Self End TH', 'Left Hairpin TH', 'Right Hairpin TH',
                    'Left End Stability', 'Right End Stability', 'Pair Compl Any TH', 'Pair Compl End TH',
                    'Pair Product Size']
        with open(csvfilename, 'a') as f:
            writer = csv.DictWriter(f, fieldnames=headings)
            writer.writerows(fprimers)
            f.close()