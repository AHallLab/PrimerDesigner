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
csvfilename = 'Primers '+datetime.today().strftime('%y-%m-%d %H.%M.%S')+'.csv'
txtfilename = 'CDS with bsaI '+datetime.today().strftime('%y-%m-%d %H.%M.%S')+'.txt'
headings = ['Gene', 'Primer', 'bsaI in Primer', 'Flank', 'Pair Penalty', 'Left Penalty', 'Right Penalty', 'Primer Forward', 'Primer Reverse',
        'Left (Start, Length)', 'Right (Start, Length)', 'Left TM', 'Right TM', 'Left GC%', 'Right GC%', 'Left Self Any TH',
        'Right Self Any TH', 'Left Self End TH', 'Right Self End TH', 'Left Hairpin TH', 'Right Hairpin TH',
        'Left End Stability', 'Right End Stability', 'Pair Compl Any TH', 'Pair Compl End TH', 'Pair Product Size']
with open(csvfilename, 'w') as f:
    writer = csv.DictWriter(f, fieldnames=headings)
    writer.writeheader()
    f.close()

# using BioPython to get info from fasta file.
# 'Whole Chromosome CDS list 1000bp.fasta' 'Test_fastas/Test_Genes.fasta'
myfast = SeqIO.parse(str(args.input), 'fasta')
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

        # To be used in the start and length for excluded regions
        exclude_start_left = CDS_start+200 # Start point for excluded Left primer region.
        # Depending on CDS length either use final position or -200 for the length element of excluded Right primer
        # region. Depending if CDS Length is < 200
        if CDS_len < 200:
            exclude_start_right = CDS_len
        else:
            exclude_start_right = CDS_len-200
        #print(CDS_start+(exclude_start_right))

        leftprimerlist = (primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': ID,
                'SEQUENCE_TEMPLATE': SEQ,
                'SEQUENCE_EXCLUDED_REGION': [exclude_start_left, len(SEQ)-exclude_start_left]
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
        #print(leftprimerlist)
        # leftprimerlist is the initial output of primer3

        leftnewlist = cleanupdata(leftprimerlist)
        #print(leftnewlist)
        # leftnewlist is a filtered version of leftprimerlist, number of primers returned and explains have been removed

        # Separate the data based on the primer number in each key! This results in a dictionary of sub-dictionaries where a
        # sub-dictionary is a primer pair
        reg = re.compile(r'\d+') # regex that finds the number (primer number) in the keys (e.g. PRIMER_PAIR_0_PENALTY)
        leftsplitprimers = {} # This will be the dictionary that all sub-dictionaries go into
        for k, v in leftnewlist.items():
            m=reg.search(k)
            if m:
                numb = m.group()
                leftsplitprimers.setdefault(numb, {})[k]=v
            else:
                print('Invalid key:', k, v)

        #print(leftsplitprimers)

        # This adds 2 new keys and values to each primer-pair dictionary. 1 for the Primer pair number, and 1 for the Gene ID
        # This also then turns the dictionary of sub-dictionaries into a LIST of dictionaries called primers, this is important
        # for further steps.

        leftprimers = []
        for j in leftsplitprimers:
            leftsplitprimers[j]['Primer'] = str(int(j)+1)
            leftsplitprimers[j]['Gene'] = ID
            leftsplitprimers[j]['Flank'] = 'Left'
            leftprimers.append(leftsplitprimers[j])

        #print(leftprimers)

        # Check for bsaI and filter primers and append key and value of 'Yes' or 'No'
        leftfprimers = []
        check = re.compile(r'PRIMER_LEFT_\d+$')
        check2 = re.compile(r'PRIMER_RIGHT_\d+$')
        for i in leftprimers:
            for j in i:
                if check.search(j):
                    leftprimerstart = i[j][0]
                elif check2.search(j):
                    leftprimerend = i[j][0]+1
                else:
                    pass
            leftprimer = SEQ[leftprimerstart:leftprimerend]
            if bsa_in_primer(leftprimer):
                i['bsaI in Primer'] = 'Yes'
                leftfprimers.append(i)
            else:
                i['bsaI in Primer'] = 'No'
                leftfprimers.append(i)

        #print(leftfprimers)

        leftfprimers = replacekey(leftfprimers)
        # print(leftfprimers)

        # Repeat the same process for the Right Flanking Region
        rightprimerlist = (primer3.bindings.designPrimers(
            {
                'SEQUENCE_ID': ID,
                'SEQUENCE_TEMPLATE': SEQ,
                'SEQUENCE_EXCLUDED_REGION': [1, CDS_start + (exclude_start_right)]
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
        # print(rightprimerlist)

        rightnewlist = cleanupdata(rightprimerlist)
        # print(rightnewlist)

        reg = re.compile(r'\d+')
        rightsplitprimers = {}
        for k, v in rightnewlist.items():
            m = reg.search(k)
            if m:
                numb = m.group()
                rightsplitprimers.setdefault(numb, {})[k] = v
            else:
                print('Invalid key:', k, v)

        # print(rightsplitprimers)

        rightprimers = []
        for j in rightsplitprimers:
            rightsplitprimers[j]['Primer'] = str(int(j) + 1)
            rightsplitprimers[j]['Gene'] = ID
            rightsplitprimers[j]['Flank'] = 'Right'
            rightprimers.append(rightsplitprimers[j])

        # print(rightprimers)

        rightfprimers = []
        check = re.compile(r'PRIMER_LEFT_\d+$')
        check2 = re.compile(r'PRIMER_RIGHT_\d+$')
        for i in rightprimers:
            for j in i:
                if check.search(j):
                    rightprimerstart = i[j][0]
                elif check2.search(j):
                    rightprimerend = i[j][0] + 1
                else:
                    pass
            rightprimer = SEQ[rightprimerstart:rightprimerend]
            if bsa_in_primer(rightprimer):
                i['bsaI in Primer'] = 'Yes'
                rightfprimers.append(i)
            else:
                i['bsaI in Primer'] = 'No'
                rightfprimers.append(i)

        rightfprimers = replacekey(rightfprimers)
        # print(rightfprimers)

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
            writer.writerows(leftfprimers)
            writer.writerows(rightfprimers)
            f.close()