from Bio import SeqIO
from pathlib import Path
import time


# Project paths
BASE_DIR = Path(__file__).resolve().parent.parent
FASTA_FILE = BASE_DIR / "data" / "sample.fasta"


def load_fasta(file_path):
    """
    Read DNA sequences from a FASTA file.
    """

    sequences = {}

    for record in SeqIO.parse(str(file_path), "fasta"):

        sequence = str(record.seq).upper()

        # Keep only valid DNA bases
        sequence = "".join(
            base for base in sequence
            if base in "ATCGN"
        )

        sequences[record.id] = sequence

    return sequences


def find_matches(genome, target, max_mismatches=2):
    """
    Search for target DNA using a sliding window.

    Returns matching positions and mismatch counts.
    """

    matches = []

    target_length = len(target)

    if target_length == 0:
        raise ValueError("Target sequence cannot be empty.")

    if target_length > len(genome):
        return matches

    for i in range(len(genome) - target_length + 1):

        window = genome[i:i + target_length]

        mismatches = 0

        for j in range(target_length):

            if target[j] != window[j]:
                mismatches += 1

        if mismatches <= max_mismatches:

            matches.append({
                "position": i + 1,
                "sequence": window,
                "mismatches": mismatches
            })

    return matches


def find_matches_in_chunks(
    genome,
    target,
    max_mismatches=2,
    chunk_size=100000
):
    """
    Search a genome in manageable chunks.
    Overlap prevents missing matches at chunk boundaries.
    """

    all_matches = []
    target_length = len(target)

    if target_length == 0:
        raise ValueError("Target cannot be empty.")

    for chunk_start in range(0, len(genome), chunk_size):

        # Include overlap so boundary-spanning matches are detected
        chunk = genome[
            chunk_start:chunk_start + chunk_size + target_length - 1
        ]

        # Only process windows starting inside this chunk's main region
        valid_starts = min(
            chunk_size,
            len(genome) - chunk_start - target_length + 1
        )

        for i in range(valid_starts):

            window = chunk[i:i + target_length]

            mismatches = sum(
                target[j] != window[j]
                for j in range(target_length)
            )

            if mismatches <= max_mismatches:

                all_matches.append({
                    "position": chunk_start + i + 1,
                    "sequence": window,
                    "mismatches": mismatches
                })

    return all_matches



def display_results(all_matches):

    print("\n" + "=" * 65)
    print("             GENEWEAVER - ALIGNMENT RESULTS")
    print("=" * 65)

    total_matches = 0

    for chromosome, matches in all_matches.items():

        print(f"\nSequence: {chromosome}")

        for match in matches:

            total_matches += 1

            print(
                f"Position: {match['position']} | "
                f"Sequence: {match['sequence']} | "
                f"Mismatches: {match['mismatches']}"
            )

    print("\n" + "-" * 65)
    print(f"Total potential matches: {total_matches}")
    print("=" * 65)


def main():

    print("\nStarting GeneWeaver CPU Engine...")

    target = input(
        "Enter target DNA sequence (A, T, C, G): "
    ).strip().upper()

    if not target or any(base not in "ATCG" for base in target):
        print("Error: Enter a non-empty DNA sequence containing only A, T, C, G.")
        return

    max_mismatches = 2

    if not FASTA_FILE.exists():
        print(f"FASTA file not found: {FASTA_FILE}")
        return

    print("\nLoading FASTA file...")

    start_time = time.perf_counter()

    sequences = load_fasta(FASTA_FILE)

    all_matches = {}

    for chromosome, genome in sequences.items():
        matches = find_matches_in_chunks(
            genome,
            target,
            max_mismatches,
            chunk_size=100000
        )

        all_matches[chromosome] = matches

    end_time = time.perf_counter()

    display_results(all_matches)

    print(f"\nExecution time: {end_time - start_time:.6f} seconds")
    print("GeneWeaver CPU scan completed!")


if __name__ == "__main__":
    main()
