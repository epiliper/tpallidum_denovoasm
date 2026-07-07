#!/usr/bin/env python3

import regex
from Bio import SeqIO
from Bio.Seq import Seq

def get_terminal_seqs(ref_fasta, term_len):
    record = next(SeqIO.parse(ref_fasta, "fasta"))
    seq = str(record.seq).upper()
    if len(seq) < term_len * 2:
        raise ValueError(
            f"Reference sequence ({len(seq)} bp) is shorter than "
            f"2x terminal length ({term_len * 2} bp)."
        )
    n_term = seq[:term_len]
    c_term = seq[-term_len:]
    return n_term, c_term

def split_reads_at_termini(reads_fasta: str, n_term: str, c_term: str,
                            max_edits: int, out_prefix: str, check_rc: bool = True):
    n_pattern = f"({regex.escape(n_term)}){{e<={max_edits}}}"
    c_pattern = f"({regex.escape(c_term)}){{e<={max_edits}}}"

    out_path = f"{out_prefix}_split.fasta"
    n_records = both_split = 0

    with open(out_path, "w") as out_f:
        for record in SeqIO.parse(reads_fasta, "fasta"):
            n_records += 1
            fwd = str(record.seq).upper()
            candidates = [("+", fwd)]
            if check_rc:
                candidates.append(("-", str(Seq(fwd).reverse_complement())))

            # Use whichever strand actually has a hit (prefer + on ties/no hits)
            seq, strand_used = fwd, "+"
            for strand_label, s in candidates:
                if regex.search(n_pattern, s, flags=regex.IGNORECASE) or \
                   regex.search(c_pattern, s, flags=regex.IGNORECASE):
                    seq, strand_used = s, strand_label
                    break

            match_5 = regex.search(n_pattern, seq, flags=regex.IGNORECASE)
            match_3 = regex.search(c_pattern, seq, flags=regex.IGNORECASE)

            start5, _end5 = match_5.span() if match_5 else (None, None)
            _start3, end3 = match_3.span() if match_3 else (None, None)

            pieces = []
            if match_5:
                pieces.append((f"{out_prefix}_{record.id}|5trim|strand={strand_used}", seq[start5:]))
            if match_3:
                pieces.append((f"{out_prefix}_{record.id}|3trim|strand={strand_used}", seq[:end3]))

            # if we found both 5' and 3' then split
            if match_5 and match_3:
                for name, s in pieces:
                    if s:
                        out_f.write(f">{name}\n{s}\n")

            # just write original if there isn't a hit
            else:
                 out_f.write(f">{record.id}\n{record.seq}\n")

            if match_5 and match_3:
                both_split += 1

    print(f"  both termini (junction) split: {both_split}")
    print(f"Output written to {out_path}")

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("--contigs", required = True, type = str, help = "contigs")
    parser.add_argument("--ref", required = True, type = str, help = "reference fasta")
    parser.add_argument("--terminal_length", required = False, type = int, help = "length of ends extracted from ref used in regex search", default = 40)
    parser.add_argument("--max_edit_distance", required = False, default = 5, type = int, help = "max edit distance of match (shared by both 5' and 3)")
    parser.add_argument("--output_prefix", required = True, type = str, help = "output prefix for emitted files")
    args = parser.parse_args()

    n_term, c_term = get_terminal_seqs(args.ref, args.terminal_length)
    print(n_term, c_term)

    split_reads_at_termini(args.contigs, n_term, c_term, args.max_edit_distance, args.output_prefix, True)
