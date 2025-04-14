from Bio import Entrez, SeqIO
from pathlib import Path
import re
import csv
import time


species_lut = {'human': 'Homo sapiens',           'monkey': 'Macaca mulatta',
               'mouse': 'Mus musculus',           'rat':    'Rattus norvegicus',
               'dog':   'Canis familiaris',       'rabbit': 'Oryctolagus cuniculus'}

Entrez.email = "lalondematthew@hotmail.com"
parent_path = Path("C:\\Users\\mlalonde\\OneDrive - Kymera Therapeutics\\Desktop\\IRF5 Project\\Alignments")
protein = 'IRF5'
genus_species = species_lut['monkey']
timestamp = time.strftime("%Y%m%d-%H%M%S")
file_prefix = f"{genus_species.replace(' ', '-')}_{protein}"


# Function to fetch protein sequences and save to a FASTA file
def download_irf5_isoforms():
    output_file = parent_path / Path(f"{file_prefix}_isoforms_{timestamp}.fasta")
    # Search for IRF5 protein isoforms in the NCBI database
    search_term = f"{protein}[Gene Name] AND {genus_species}[Organism] AND RefSeq[filter]"
    handle = Entrez.esearch(db="protein", term=search_term, retmax=20)  # Adjust retmax if needed
    record = Entrez.read(handle)
    handle.close()

    protein_ids = record["IdList"]  # Get list of protein IDs

    if not protein_ids:
        print(f"No protein sequences found for {protein}.")
        return

    # Fetch sequences using IDs
    handle = Entrez.efetch(db="protein", id=protein_ids, rettype="fasta", retmode="text")
    sequences = handle.read()
    handle.close()

    # Save to FASTA file
    with open(output_file, "w") as f:
        f.write(sequences)

    return output_file


def remove_duplicate_sequences(input_fasta):
    """
    Reads a FASTA file, removes identical sequences (regardless of sequence names),
    and writes unique sequences to a new FASTA file.

    :param input_fasta: Path to the input FASTA file.
    """
    output_file = parent_path / Path(f"{file_prefix}_unique_isoforms_{timestamp}.fasta")
    unique_seqs = {}

    # Read the input FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_str = str(record.seq)  # Convert sequence to string for comparison

        # Store the first occurrence of each unique sequence
        if seq_str not in unique_seqs:
            unique_seqs[seq_str] = record

    # Write unique sequences to output FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(unique_seqs.values(), output_handle, "fasta")

    return output_file


def tryptic_digest(sequence):
    """
    Simulates trypsin digestion of a protein sequence.
    Trypsin cleaves at K (Lysine) and R (Arginine), unless followed by P (Proline).

    :param sequence: Protein sequence as a string.
    :return: List of tryptic peptides (filtered to length >6).
    """
    # Regular expression for trypsin cleavage: Split at K or R, but not if followed by P
    peptides = re.split(r'(?<=[KR])(?!P)', sequence)

    # Filter peptides to only include those >6 residues
    return [pep for pep in peptides if len(pep) > 6]


def digest_fasta_and_tabulate(input_fasta):
    """
    Reads a FASTA file, performs tryptic digestion on each protein sequence,
    and outputs a CSV file with unique peptides and their corresponding isoforms.

    :param input_fasta: Path to input FASTA file.
    :param output_csv: Path to save the peptide-isoform table (CSV).
    """
    output_file = parent_path / Path(f"{file_prefix}_tryptic_pp_composition_{timestamp}.csv")
    peptide_dict = {}

    # Read input FASTA and perform digestion
    for record in SeqIO.parse(input_fasta, "fasta"):
        protein_id = record.id
        peptides = tryptic_digest(str(record.seq))

        for peptide in peptides:
            if peptide not in peptide_dict:
                peptide_dict[peptide] = set()  # Use a set to avoid duplicate entries
            peptide_dict[peptide].add(protein_id)

    # Write results to CSV file
    with open(output_file, "w", newline="") as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(["Peptide", "Protein Isoforms"])  # Header

        for peptide, isoforms in sorted(peptide_dict.items()):
            csv_writer.writerow([peptide, "; ".join(sorted(isoforms))])

    return output_file


def main():
    if __name__ == "__main__":
        fas1 = download_irf5_isoforms()
        fas2 = remove_duplicate_sequences(fas1)
        csv1 = digest_fasta_and_tabulate(fas2)
        print(csv1)

main()
