import argparse
import pandas as pd
import os
from Bio import AlignIO
from mhcflurry import Class1AffinityPredictor

# --- Configuration ---
# Expanded list of common MHC Class I alleles for broader coverage.
COMPREHENSIVE_ALLELES = [
    "HLA-A*01:01", "HLA-A*02:01", "HLA-A*02:03", "HLA-A*02:06", "HLA-A*03:01",
    "HLA-A*11:01", "HLA-A*23:01", "HLA-A*24:02", "HLA-A*26:01", "HLA-A*30:01",
    "HLA-A*30:02", "HLA-A*31:01", "HLA-A*32:01", "HLA-A*68:01", "HLA-A*68:02",
    "HLA-B*07:02", "HLA-B*08:01", "HLA-B*15:01", "HLA-B*35:01", "HLA-B*40:01",
    "HLA-B*44:02", "HLA-B*44:03", "HLA-B*51:01", "HLA-B*53:01", "HLA-B*57:01",
    "HLA-B*58:01", "HLA-C*03:04", "HLA-C*04:01", "HLA-C*05:01", "HLA-C*06:02",
    "HLA-C*07:01", "HLA-C*07:02"
]
BINDING_THRESHOLD_NM = 500 # Peptides < 500 nM are considered potential binders.
OUTPUT_DIR = "results"

def get_ref_and_map(fasta_file):
    """
    Extracts the first sequence from the alignment, creates a gapless version,
    and returns a map of alignment coordinates to sequence coordinates.
    """
    try:
        alignment = AlignIO.read(fasta_file, "fasta")
        ref_record = alignment[0]
        print(f"\nProcessing file: {os.path.basename(fasta_file)}")
        print(f"Using reference sequence: {ref_record.id}")

        coord_map = {}
        seq_idx = 0
        for align_idx, char in enumerate(ref_record.seq):
            if char != '-':
                coord_map[align_idx] = seq_idx
                seq_idx += 1
        
        ref_sequence = str(ref_record.seq).replace('-', '')
        return ref_sequence, coord_map

    except FileNotFoundError:
        print(f"Error: The file '{fasta_file}' was not found.")
        return None, None
    except Exception as e:
        print(f"An error occurred while reading the FASTA file '{fasta_file}': {e}")
        return None, None

def parse_regions(region_strings, coord_map):
    """
    Parses region strings and converts alignment coordinates to sequence coordinates.
    """
    parsed = []
    for r in region_strings:
        try:
            align_start, align_end = map(int, r.split('-'))
            
            seq_start_idx = None
            for i in range(align_start - 1, align_end):
                if i in coord_map:
                    seq_start_idx = coord_map[i]
                    break
            
            seq_end_idx = None
            for i in range(align_end - 1, align_start - 2, -1):
                if i in coord_map:
                    seq_end_idx = coord_map[i]
                    break

            if seq_start_idx is not None and seq_end_idx is not None:
                parsed.append((seq_start_idx, seq_end_idx + 1))
            else:
                print(f"Warning: Region '{r}' contains only gaps in the reference. Skipping.")

        except ValueError:
            print(f"Warning: Could not parse region '{r}'. Skipping.")
    return parsed

def generate_peptides(sequence, min_len=8, max_len=11):
    """
    Generates all possible peptide fragments (k-mers).
    """
    peptides = set()
    for length in range(min_len, max_len + 1):
        for i in range(len(sequence) - length + 1):
            peptides.add(sequence[i:i+length])
    return peptides

def main():
    parser = argparse.ArgumentParser(description="Predict MHC Class I binding epitopes in conserved regions for all PRO files in a directory.")
    parser.add_argument('--directory', type=str, required=True, help="Path to the directory containing FASTA alignment files.")
    parser.add_argument('--regions', nargs='+', required=True, help='List of conserved regions, e.g., "142-159" "161-179".')
    args = parser.parse_args()

    print("Loading MHCflurry predictor...")
    try:
        predictor = Class1AffinityPredictor.load()
    except Exception as e:
        print(f"Error loading mhcflurry model: {e}\nPlease run 'mhcflurry-downloads fetch'.")
        return

    # Find all protein fasta files in the specified directory
    try:
        all_files = os.listdir(args.directory)
        pro_files = [os.path.join(args.directory, f) for f in all_files if f.endswith('_PRO.fasta')]
    except FileNotFoundError:
        print(f"Error: Directory not found at '{args.directory}'")
        return

    if not pro_files:
        print(f"No files ending with '_PRO.fasta' found in '{args.directory}'")
        return

    print(f"Found {len(pro_files)} protein alignment files to process.")

    # Loop through each file and perform the analysis
    for fasta_file in pro_files:
        ref_sequence, coord_map = get_ref_and_map(fasta_file)
        if not ref_sequence:
            continue

        conserved_regions = parse_regions(args.regions, coord_map)

        all_peptides = set()
        for start, end in conserved_regions:
            if start < len(ref_sequence) and end <= len(ref_sequence):
                region_sequence = ref_sequence[start:end]
                peptides_from_region = generate_peptides(region_sequence)
                all_peptides.update(peptides_from_region)
            else:
                print(f"Warning: Mapped region {start+1}-{end} is out of bounds. Skipping.")

        if not all_peptides:
            print("No valid peptides generated from the provided regions for this file.")
            continue

        print(f"Generated {len(all_peptides)} unique peptides for analysis.")
        print(f"Predicting binding for {len(COMPREHENSIVE_ALLELES)} common alleles...")
        
        all_predictions = []
        peptide_list = list(all_peptides)
        for allele in COMPREHENSIVE_ALLELES:
            df = predictor.predict_to_dataframe(peptides=peptide_list, allele=allele)
            all_predictions.append(df)

        predictions = pd.concat(all_predictions)
        
        strong_binders = predictions[predictions['prediction'] < BINDING_THRESHOLD_NM]
        
        print(f"\n--- Results for {os.path.basename(fasta_file)} ---")
        print(f"Found {len(strong_binders)} potential epitopes with binding affinity < {BINDING_THRESHOLD_NM} nM")

        if not strong_binders.empty:
            strong_binders_sorted = strong_binders.sort_values(by="prediction", ascending=True)
            
            os.makedirs(OUTPUT_DIR, exist_ok=True)
            
            # Create a unique output filename for each input file
            base_name = os.path.basename(fasta_file).replace('.fasta', '')
            output_filename = os.path.join(OUTPUT_DIR, f"predicted_epitopes_{base_name}.csv")
            
            strong_binders_sorted.to_csv(output_filename, index=False)
            print(f"All strong binders saved to '{output_filename}'")

            print("\nTop 20 Predicted Epitopes (sorted by best affinity):")
            print(strong_binders_sorted.head(20).to_string(
                columns=["peptide", "allele", "prediction"],
                header=["Peptide", "MHC Allele", "Affinity (nM)"],
                index=False
            ))
        else:
            print("No strong binders found with the current settings for this file.")

if __name__ == "__main__":
    main()
