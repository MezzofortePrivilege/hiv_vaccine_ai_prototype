import argparse
import os
from Bio import AlignIO, SeqIO # Correctly import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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

def parse_regions_to_seq_coords(region_strings, coord_map):
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
                parsed.append((seq_start_idx, seq_end_idx + 1, r)) # Also keep original region string
            else:
                print(f"Warning: Region '{r}' contains only gaps in the reference. Skipping.")

        except ValueError:
            print(f"Warning: Could not parse region '{r}'. Skipping.")
    return parsed

def main():
    parser = argparse.ArgumentParser(description="Prepare conserved region sequences for NetMHCIIpan analysis.")
    parser.add_argument('--directory', type=str, required=True, help="Path to the directory containing FASTA alignment files.")
    parser.add_argument('--regions', nargs='+', required=True, help='List of conserved regions, e.g., "142-159" "161-179".')
    args = parser.parse_args()

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
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Loop through each file and perform the analysis
    for fasta_file in pro_files:
        ref_sequence, coord_map = get_ref_and_map(fasta_file)
        if not ref_sequence:
            continue

        conserved_regions_coords = parse_regions_to_seq_coords(args.regions, coord_map)
        
        output_records = []
        for start, end, region_str in conserved_regions_coords:
            if start < len(ref_sequence) and end <= len(ref_sequence):
                region_sequence = ref_sequence[start:end]
                
                # Create a SeqRecord for the FASTA output
                record_id = f"ConservedRegion_{region_str}"
                record = SeqRecord(Seq(region_sequence), id=record_id, description=f"Source: {os.path.basename(fasta_file)}")
                output_records.append(record)
            else:
                print(f"Warning: Mapped region {start+1}-{end} is out of bounds. Skipping.")

        if not output_records:
            print("No valid regions found to extract for this file.")
            continue

        # Create a unique output filename for each input file
        base_name = os.path.basename(fasta_file).replace('.fasta', '')
        output_filename = os.path.join(OUTPUT_DIR, f"mhc2_candidates_from_{base_name}.fasta")

        # Write the sequences to the output FASTA file
        SeqIO.write(output_records, output_filename, "fasta")
        print(f"Successfully extracted {len(output_records)} conserved sequences to '{output_filename}'")
        print("--> Next step: Upload this file to the NetMHCIIpan web server.")

if __name__ == "__main__":
    main()
