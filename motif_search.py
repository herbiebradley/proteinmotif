from __future__ import print_function

import urllib
import re

from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO

def get_match_list(protein_string):
    print(protein_string)
    # Compile regex from regex string, deal with multiline strings:
    regex = re.compile(r"[^JORUZXP][^JORUZXPKRHW][VLSWFNQ][ILTYWFN][FIY][^JORUZXPKRH]", flags=re.MULTILINE)
    # Create iterator of matches including start and end positions:
    matches = regex.finditer(protein_string)
    # Create match list of the repeating form [match, start position, end position, ... ]:
    match_list = []
    for match in matches:
        match_list.append(match.group())
        # Get match positions, convert to Bioinformatics counting:
        match_list.append(match.start() + 1)
        match_list.append(match.end() + 1)
    print(match_list)
    return match_list

if __name__ == "__main__":

    accesion_code = input("Type accession code: ")
    try:
        with ExPASy.get_sprot_raw(accesion_code) as handle:
            record = SwissProt.read(handle)
        protein_string = record.sequence
        match_list = get_match_list(protein_string)
    except urllib.error.HTTPError:
        print("No protein record found, please enter valid accession code.")
