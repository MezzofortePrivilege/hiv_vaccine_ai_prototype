import argparse
import pandas as pd
from Bio import SeqIO
from mhcflurry import Class1PresentationPredictor

# --- Configuration ---
# A list of common MHC Class I alleles to predict binding for.
# This list can be expanded based on target populations.
# You can find more alleles at: http://www.iedb.org/allele.php
COMMON_ALLELES = [
    "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-A*11:01",
    "HLA-A*24:02", "HLA-B*07:02", "HLA-B*08:01", "HLA-B*44:02"
]

# Peptides with a predicted binding affinity < 500 nM are typically
# considered potential binders. We'll use a stricter threshold.
BINDING_THRESHOLD_NM = 250

def extract_reference_sequence(fasta_file):
    """
    Extracts the first sequence from the FASTA file to use as a reference.
    This is typically the HXB2 reference sequence in Los Alamos DB files.

    Args:
        fasta_file (str): Path to the FASTA alignment file.

    Returns:
        str: The reference protein sequence, with gaps removed.
    """
    try:
        first_record = next(SeqIO.parse(fasta_file, "fasta"))
        print(f"Using reference sequence: {first_record.id}")
        return str(first_record.seq).replace('-', '')
    except FileNotFoundError:
        print(f"Error: The file '{fasta_file}' was not found.")
        return None
    except Exception as e:
        print(f"An error occurred while reading the FASTA file: {e}")
        return None

def parse_regions(region_strings):
    """
    Parses region strings like "142-159" into a list of (start, end) tuples.

    Args:
        region_strings (list): A list of strings, e.g., ["142-159", "161-179"].

    Returns:
        list: A list of tuples, e.g., [(141, 159), (160, 179)]. Uses 0-based index.
    """
    parsed = []
    for r in region_strings:
        try:
            start, end = map(int, r.split('-'))
            # Convert from 1-based inclusive to 0-based exclusive for slicing
            parsed.append((start - 1, end))
        except ValueError:
            print(f"Warning: Could not parse region '{r}'. Skipping.")
    return parsed

def generate_peptides(sequence, min_len=8, max_len=11):
    """
    Generates all possible peptide fragments (k-mers) of specified lengths
    from a given sequence.

    Args:
        sequence (str): The protein sequence segment.
        min_len (int): Minimum peptide length.
        max_len (int): Maximum peptide length.

    Returns:
        set: A set of unique peptide strings.
    """
    peptides = set()
    for length in range(min_len, max_len + 1):
        for i in range(len(sequence) - length + 1):
            peptides.add(sequence[i:i+length])
    return peptides

def main():
    """
    Main function to run epitope prediction.
    """
    parser = argparse.ArgumentParser(description="Predict MHC Class I binding epitopes in conserved regions.")
    parser.add_argument('--file', type=str, required=True, help="Path to the FASTA alignment file.")
    parser.add_argument('--regions', nargs='+', required=True, help='List of conserved regions, e.g., "142-159" "161-179".')

    args = parser.parse_args()

    # 1. Load the MHCflurry predictor
    print("Loading MHCflurry predictor... (This may take a moment on first run)")
    try:
        predictor = Class1PresentationPredictor.load()
    except Exception as e:
        print(f"Error loading mhcflurry model: {e}")
        print("Please ensure you have run 'mhcflurry-downloads fetch' after installation.")
        return

    # 2. Extract the reference sequence
    ref_sequence = extract_reference_sequence(args.file)
    if not ref_sequence:
        return

    # 3. Parse the user-provided regions
    conserved_regions = parse_regions(args.regions)

    # 4. Generate all candidate peptides from the conserved regions
    all_peptides = set()
    for start, end in conserved_regions:
        if start < len(ref_sequence) and end <= len(ref_sequence):
            region_sequence = ref_sequence[start:end]
            peptides_from_region = generate_peptides(region_sequence)
            all_peptides.update(peptides_from_region)
        else:
            print(f"Warning: Region {start+1}-{end} is out of bounds for the reference sequence (length {len(ref_sequence)}). Skipping.")

    if not all_peptides:
        print("No valid peptides generated from the provided regions.")
        return

    print(f"\nGenerated {len(all_peptides)} unique peptides from conserved regions for analysis.")

    # 5. Make predictions
    print(f"Predicting binding for {COMMON_ALLELES}...")
    predictions = predictor.predict(peptides=list(all_peptides), alleles=COMMON_ALLELES)

    # 6. Filter for strong binders and display results
    strong_binders = predictions[predictions['presentation_score'] > 0.5] # mhcflurry uses a score from 0-1
    
    print(f"\n--- Results ---")
    print(f"Found {len(strong_binders)} potential epitopes with a presentation score > 0.5")

    if not strong_binders.empty:
        # Sort for better presentation
        strong_binders_sorted = strong_binders.sort_values(by="presentation_score", ascending=False)
        
        # Display the top 20 hits for brevity
        print("\nTop 20 Predicted Epitopes:")
        print(strong_binders_sorted.head(20).to_string())
    else:
        print("No strong binders found with the current settings.")


if __name__ == "__main__":
    main()
