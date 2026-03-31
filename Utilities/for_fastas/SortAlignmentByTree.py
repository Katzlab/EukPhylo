#Author, date: Ozan and Claude, March 31 2026
#Motivation: Sorts a multiple sequence alignment (FASTA) so that the order of sequences matches the top-to-bottom tip order of an input phylogenetic tree (Newick or Nexus format).
#Dependencies: Python3, BioPython
#Inputs: A fasta file to be sorted and a tree file to be used as a reference
#Outputs: A sorted fasta file 
#Example: python3 SortAlignmentByTree.py -a alignment.fasta -t tree_file.nwk [-f <tree_format>] [-o <output.fasta>] [--missing ignore|warn|error]
#Note: The arguments in brackets are optional. 
#Note: The script will compare sequence IDs in the tree and fasta files. In case of mismatches, it prints a warning by default. This behaviour can be modified by changing the value of "--missing".
#Note: If fasta file contains sequences that are missing from the tree file, and if "ignore" or "warn" is specified, the script will append these at the end of the alignment.

#Dependencies
import argparse
import os
import sys

from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Tree


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def get_args():
    """Parse and return command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            "Sort a multiple sequence alignment (FASTA) to match the "
            "top-to-bottom tip order in a phylogenetic tree file "
            "(Newick or Nexus)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python SortAlignmentByTree.py -a seqs.fasta -t tree.nwk\n"
            "  python SortAlignmentByTree.py -a seqs.fasta -t tree.nex "
            "-f nexus -o sorted.fasta\n"
            "  python SortAlignmentByTree.py -a seqs.fasta -t tree.nwk "
            "--missing warn\n"
        ),
    )

    parser.add_argument(
        "-a", "--alignment",
        required=True,
        help="Path to the input alignment file (FASTA format).",
    )
    parser.add_argument(
        "-t", "--tree",
        required=True,
        help="Path to the input tree file (Newick or Nexus).",
    )
    parser.add_argument(
        "-f", "--format",
        dest="tree_format",
        default=None,
        choices=["newick", "nexus"],
        help=(
            "Tree file format. If omitted, the format is auto-detected from "
            "the file extension."
        ),
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help=(
            "Path and filename for the sorted output FASTA"
        ),
    )
    parser.add_argument(
        "--missing",
        default="warn",
        choices=["ignore", "warn", "error"],
        help=(
            "Behaviour when a tree tip has no matching sequence (or vice "
            "versa). 'ignore' silently skips mismatches; 'warn' prints a "
            "warning (default); 'error' aborts."
        ),
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Functions for detecting tree format and getting tip order
# ---------------------------------------------------------------------------

def detect_tree_format(tree_path: str) -> str:

    newick_exts = {".nwk", ".tree", ".newick", ".tre"}
    nexus_exts  = {".nex", ".nexus", ".nxs"}

    ext = os.path.splitext(tree_path)[-1].lower()

    if ext in newick_exts:
        return "newick"
    if ext in nexus_exts:
        return "nexus"

    raise ValueError(
        f"Cannot infer tree format from extension '{ext}'. "
        "Please supply -f newick or -f nexus explicitly."
    )


def get_tip_order(tree: Tree) -> list:

    return [tip.name for tip in tree.get_terminals()]


# ---------------------------------------------------------------------------
# Function for mapping sequence IDs in the fasta file
# ---------------------------------------------------------------------------

def load_alignment(aln_path: str) -> dict:

    records = {}
    for rec in SeqIO.parse(aln_path, "fasta"):
        records[rec.id] = rec
    return records


# ---------------------------------------------------------------------------
# Function for sorting the fasta file
# ---------------------------------------------------------------------------

def sort_alignment(tip_order: list, records: dict, missing: str) -> list:

    sorted_records = []
    seen_in_tree = set()

    for tip in tip_order:
        seen_in_tree.add(tip)
        if tip in records:
            sorted_records.append(records[tip])
        else:
            msg = (
                f"Tree tip '{tip}' has no matching sequence in the alignment."
            )
            if missing == "error":
                sys.exit(f"[ERROR] {msg}")
            elif missing == "warn":
                print(f"[WARNING] {msg}", file=sys.stderr)
            # 'ignore' → do nothing

    # Report sequences in the alignment that are absent from the tree
    orphan_seqs = [sid for sid in records if sid not in seen_in_tree]
    if orphan_seqs:
        msg = (
            f"{len(orphan_seqs)} sequence(s) in the alignment have no "
            f"matching tip in the tree and will be appended at the end: "
            f"{', '.join(orphan_seqs)}"
        )
        if missing == "error":
            sys.exit(f"[ERROR] {msg}")
        elif missing == "warn":
            print(f"[WARNING] {msg}", file=sys.stderr)
        for sid in orphan_seqs:
            sorted_records.append(records[sid])

    return sorted_records


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def build_output_path(aln_path: str) -> str:

    stem, ext = os.path.splitext(aln_path)
    return f"{stem}_sorted{ext}"


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

def main():
    args = get_args()

    # Resolve tree format
    tree_format = args.tree_format
    if tree_format is None:
        try:
            tree_format = detect_tree_format(args.tree)
        except ValueError as exc:
            sys.exit(f"[ERROR] {exc}")

    # Load tree
    print(f"[INFO] Reading tree from '{args.tree}' (format: {tree_format}).")
    try:
        tree = Phylo.read(args.tree, tree_format)
    except Exception as exc:
        sys.exit(f"[ERROR] Could not parse tree file: {exc}")

    tip_order = get_tip_order(tree)
    print(f"[INFO] Tree has {len(tip_order)} terminal tips.")

    # Load alignment
    print(f"[INFO] Reading alignment from '{args.alignment}'.")
    try:
        records = load_alignment(args.alignment)
    except FileNotFoundError:
        sys.exit(f"[ERROR] Alignment file not found: '{args.alignment}'")
    except Exception as exc:
        sys.exit(f"[ERROR] Could not parse alignment: {exc}")

    print(f"[INFO] Alignment contains {len(records)} sequences.")

    # Sort fasta
    sorted_records = sort_alignment(tip_order, records, args.missing)
    print(
        f"[INFO] Sorted alignment contains {len(sorted_records)} sequences."
    )

    # Write output 
    output_path = args.output if args.output else build_output_path(args.alignment)

    try:
        written = SeqIO.write(sorted_records, output_path, "fasta")
    except Exception as exc:
        sys.exit(f"[ERROR] Could not write output file: {exc}")

    print(f"[INFO] Wrote {written} sequences to '{output_path}'.")


if __name__ == "__main__":
    main()
