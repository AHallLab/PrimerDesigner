import re

### FUNCTIONS
# function to search the input fasta for the first uppercase letter. Input fasta seq has the CDS in all uppercase and
# flanking regions in lower case.
def find_CDS_start(sequence):
    for i in range(len(sequence)):
        if sequence[i].isupper():
            return sequence[i]

# function that looks for bsaI in the coding sequence.
def bsa_in_CDS(sequence):
    if re.search(r'GGTCTC', sequence) or re.search(r'GAGACC', sequence):
        return True
    else:
        return False

# function that looks for bsa in the primers.
def bsa_in_primer(primer):
    if re.search(r'ggtctc', primer) or re.search(r'gagacc', primer):
        return True
    else:
        return False

# function to remove some of the clunkier bits of info from primer3 output
def cleanupdata(list_of_primers):
    newlist = list_of_primers
    toremove = ['PRIMER_LEFT_EXPLAIN', 'PRIMER_RIGHT_EXPLAIN', 'PRIMER_PAIR_EXPLAIN', 'PRIMER_LEFT_NUM_RETURNED',
                'PRIMER_RIGHT_NUM_RETURNED', 'PRIMER_INTERNAL_NUM_RETURNED', 'PRIMER_PAIR_NUM_RETURNED']
    for key in toremove:
        del newlist[key]
    return newlist

# function to replace the Primer3 keys with our own so that each dictionary uses same keys, for easy writing to csv
def replacekey(finallist):
    replacekeys = ['Pair Penalty', 'Left Penalty', 'Right Penalty', 'Primer Forward', 'Primer Reverse',
                   'Left (Start, Length)', 'Right (Start, Length)', 'Left TM', 'Right TM', 'Left GC%', 'Right GC%',
                   'Left Self Any TH', 'Right Self Any TH', 'Left Self End TH', 'Right Self End TH', 'Left Hairpin TH',
                   'Right Hairpin TH', 'Left End Stability', 'Right End Stability', 'Pair Compl Any TH',
                   'Pair Compl End TH', 'Pair Product Size', 'Primer', 'Gene', 'Flank', 'bsaI in Primer']
    for i in range(len(finallist)):
        curr_keys = list(finallist[i].copy().keys())
        for j in range(len(replacekeys)):
            finallist[i][replacekeys[j]] = finallist[i].pop(str(curr_keys[j]))
    return finallist
