import argparse
import numpy as np
from Bio import AlignIO
from math import log

def calculate_shannon_entropy(alignment):
    """
    Calculates the Shannon entropy for each column in a multiple sequence alignment.

    Args:
        alignment (MultipleSeqAlignment): A Biopython alignment object.

    Returns:
        list: A list of entropy values for each position in the alignment.
    """
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()
    entropies = []

    print(f"Analyzing alignment with {num_sequences} sequences and length {alignment_length}...")

    for i in range(alignment_length):
        column = alignment[:, i]
        frequencies = {}
        for char in column:
            if char != '-':
                frequencies[char] = frequencies.get(char, 0) + 1
        
        num_valid_chars = sum(frequencies.values())
        if num_valid_chars == 0:
            entropies.append(0)
            continue

        entropy = 0.0
        for char in frequencies:
            p_i = frequencies[char] / num_valid_chars
            if p_i > 0:
                entropy -= p_i * log(p_i, 2)
        
        entropies.append(entropy)

    print("Entropy calculation complete.")
    return entropies

def find_conserved_regions(entropies, threshold, min_length=15):
    """
    Identifies contiguous regions where the entropy is below a certain threshold.
    Default min_length is now 15, suitable for MHC-II epitopes.

    Args:
        entropies (list): A list of entropy values.
        threshold (float): The entropy threshold to define conservation.
        min_length (int): The minimum length for a region to be considered.

    Returns:
        list: A list of tuples, where each tuple represents the start and end
              of a conserved region (e.g., [(10, 25), (45, 55)]).
    """
    conserved_regions = []
    in_region = False
    start = 0

    for i, entropy in enumerate(entropies):
        if entropy < threshold and not in_region:
            in_region = True
            start = i
        elif entropy >= threshold and in_region:
            in_region = False
            end = i - 1
            if (end - start + 1) >= min_length:
                conserved_regions.append((start + 1, end + 1))

    if in_region:
        end = len(entropies) - 1
        if (end - start + 1) >= min_length:
            conserved_regions.append((start + 1, end + 1))

    return conserved_regions

def main():
    """
    Main function to parse arguments and run the analysis.
    """
    parser = argparse.ArgumentParser(description="Analyze HIV sequence conservation.")
    parser.add_argument('--file', type=str, required=True, help="Path to the FASTA alignment file.")
    parser.add_argument('--threshold', type=float, default=0.2, help="Entropy threshold for defining a conserved region. Lower is more conserved.")
    parser.add_argument('--min_length', type=int, default=15, help="Minimum length of a conserved region to report.")
    
    args = parser.parse_args()

    try:
        alignment = AlignIO.read(args.file, "fasta")
    except FileNotFoundError:
        print(f"Error: The file '{args.file}' was not found.")
        return
    except Exception as e:
        print(f"An error occurred while reading the alignment file: {e}")
        return

    entropies = calculate_shannon_entropy(alignment)
    conserved_regions = find_conserved_regions(entropies, args.threshold, args.min_length)

    print(f"\n--- Results ---")
    print(f"Found {len(conserved_regions)} conserved regions with entropy < {args.threshold} and minimum length {args.min_length}:")
    if conserved_regions:
        for start, end in conserved_regions:
            print(f"  - Region: {start}-{end} (Length: {end - start + 1})")
    else:
        print("No conserved regions found at this threshold. Try increasing it.")

if __name__ == "__main__":
    main()
