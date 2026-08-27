#!/usr/bin/env python3

"""Trim query overhangs from a two-sequence FASTA alignment."""

import argparse
import sys
from contextlib import ExitStack
from typing import TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def trim_query_overhangs(records: list[SeqRecord]) -> SeqRecord:
    """Remove terminal columns in which the reference sequence is gapped."""
    if len(records) != 2:
        raise ValueError(
            f"Expected exactly two aligned sequences, but found {len(records)}"
        )

    reference, query = records
    if len(reference) != len(query):
        raise ValueError(
            "Reference and query must have the same aligned length "
            f"({len(reference)} != {len(query)})"
        )

    reference_sequence = str(reference.seq)
    gap_characters = {"-"}
    if not reference_sequence or all(
        base in gap_characters for base in reference_sequence
    ):
        raise ValueError("Reference sequence contains only gaps")

    start = 0
    while reference_sequence[start] in gap_characters:
        start += 1

    end = len(reference_sequence)
    while reference_sequence[end - 1] in gap_characters:
        end -= 1

    query_seq: str = str(query.seq)
    query.seq = Seq("".join([b for b in query_seq[start:end] if b.isalpha()])) # join back and remove gaps

    print(f"ref start: {start}, ref end: {end}", file = sys.stderr)
    return query


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Trim query overhangs from a two-sequence MAFFT FASTA alignment. "
            "The first sequence is treated as the reference."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="input aligned FASTA (default: standard input)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="-",
        help="output trimmed aligned FASTA (default: standard output)",
    )
    return parser.parse_args()


def open_input(path: str, stack: ExitStack) -> TextIO:
    return sys.stdin if path == "-" else stack.enter_context(open(path))


def open_output(path: str, stack: ExitStack) -> TextIO:
    return sys.stdout if path == "-" else stack.enter_context(open(path, "w"))


def main() -> None:
    args = parse_args()
    try:
        with ExitStack() as stack:
            input_handle = open_input(args.input, stack)
            records = list(SeqIO.parse(input_handle, "fasta"))
            trimmed_query = trim_query_overhangs(records)
            output_handle = open_output(args.output, stack)
            SeqIO.write(trimmed_query, output_handle, "fasta")
    except (OSError, ValueError) as error:
        raise SystemExit(f"Error: {error}") from error


if __name__ == "__main__":
    main()
