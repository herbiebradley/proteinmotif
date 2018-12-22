import urllib
import csv
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
    return match_list

def match_list_from_uniprot(accession_code):
    # Validation check for valid UniProt accession code:
    try:
        # Query SwissProt with the accession code:
        with ExPASy.get_sprot_raw(accession_code) as handle:
            record = SwissProt.read(handle)
        # Get the raw protein sequence from the SwissProt record:
        protein_string = record.sequence
        # Get list of lists of matches:
        match_list = get_matches(protein_string)
        # Generate CSV file for download:
        generate_csv(match_list)
        # Enumerate list of lists so we can get indexes in table:
        results = enumerate(match_list)
    # Exception is thrown if accession_code is not in SwissProt database:
    except urllib.error.HTTPError:
        results = None
    return results

def generate_csv(match_list):
    filename = 'results.csv'
    with open(filename, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
        writer.writerow(['Match Sequence', 'Match Start', 'Match End']) # First row
        writer.writerows(match_list) # Write each list in match_list as it's own row
