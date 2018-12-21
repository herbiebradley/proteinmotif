from __future__ import print_function

import urllib
import regex as re # Third party regex module

from Bio import ExPASy
from Bio import SwissProt
from Bio import SeqIO

def get_matches(protein_string):
    # Compile regex from regex string, deal with multiline strings:
    amyloid_regex = r"[^JORUZXP][^JORUZXPKRHW][VLSWFNQ][ILTYWFN][FIY][^JORUZXPKRH]"
    compiled_regex = re.compile(amyloid_regex, flags=re.MULTILINE)
    # Third party regex module allows easy searching for overlapping matches.
    # Create iterator of matches including start and end positions:
    matches = compiled_regex.finditer(protein_string, overlapped=True)
    # Create match list of the repeating form [match, start position, end position, ... ]:
    match_list = []
    for match in matches:
        # Get match positions, convert to Bioinformatics counting:
        match_start = match.start() + 1
        match_end = match.end() + 1
        # Create list of lists for matches:
        match_list.append([match.group(), match_start, match_end])
    # Enumerate list of lists so we can get indexes in table:
    return enumerate(match_list)

def match_list_from_uniprot(accession_code):
    # Validation check for valid UniProt accession code:
    results, error = None, None
    try:
        with ExPASy.get_sprot_raw(accession_code) as handle:
            record = SwissProt.read(handle)
        protein_string = record.sequence
        results = get_matches(protein_string)
    except urllib.error.HTTPError:
        error = "No protein record found, please enter valid accession code."
    return results, error
