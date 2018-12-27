import urllib
import csv
import regex as re # Third party regex module

from Bio import ExPASy
from Bio import SwissProt

def _get_matches(protein_string):
    """This function searches the protein for matches of the amyloid regex pattern
    and returns a list of lists with all the matches, including overlapping matches.
    This function is intended to be module-level private.

    Args:
        protein_string (String): The raw protein string as single letter codes.

    Returns:
        List: A list of lists. The inner lists each have 3 items, and represent
    a single match of the amyloid motif pattern in the protein.
    """
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

def _generate_csv(match_list):
    """Generates a .csv file from the amyloid pattern matches. The file is stored
    in the base project directory as 'results.csv'. This function is intended to
    be module-level private.

    Args:
        match_list (List): The list of lists containing the matches.
    """
    filename = 'results.csv'
    with open(filename, 'w', newline='') as outcsv:
        writer = csv.writer(outcsv, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
        writer.writerow(['Match Sequence', 'Match Start', 'Match End']) # First row
        writer.writerows(match_list) # Write each list in match_list as separate row

def match_list_from_uniprot(accession_code):
    """This function takes in a UniProt accession number from the form and returns
    None if the number is invalid, or an enumerate object containing all matches
    of the amyloid motif pattern in the protein corresponding to the accession
    number. This is the only function in the module intended to be public.

    Args:
        accession_code (String): The UniProt accession number that was entered
    into the text box in the form.

    Returns:
        None, if accession number does not represent a valid protein.
        Otherwise, returns an enumerated list of lists in which each inner list
    corresponds to a different match of the amyloid motif pattern in the protein.
    """
    # Validation check for valid UniProt accession code:
    try:
        # Query SwissProt with the accession code:
        with ExPASy.get_sprot_raw(accession_code) as handle:
            record = SwissProt.read(handle)
        # Get the raw protein sequence from the SwissProt record:
        protein_string = record.sequence
        # Get list of lists of matches:
        match_list = _get_matches(protein_string)
        # Generate CSV file for download:
        _generate_csv(match_list)
        # Enumerate list of lists so we can get indexes in table:
        results = enumerate(match_list)
    # Exception is thrown if accession_code is not in SwissProt database:
    except urllib.error.HTTPError:
        results = None
    return results
