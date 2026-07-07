#!/usr/bin/env python3

from argparse import ArgumentParser
from csv import DictReader

# column name defines
TEMPLATE = "template"
SCORE = "score"
EXPECTED = "expected"
TEMPLATE_LENGTH = "template_length"
TEMPLATE_IDENTITY = "template_identity"
TEMPLATE_COVERAGE = "template_coverage"
QUERY_IDENTITY = "query_identity"
QUERY_COVERAGE = "query_coverage"
DEPTH = "depth"
Q_VAL = "q_value"
P_VAL = "p-value"

FIELDNAMES = [TEMPLATE, SCORE, EXPECTED, TEMPLATE_LENGTH, TEMPLATE_IDENTITY, TEMPLATE_COVERAGE, QUERY_IDENTITY, QUERY_COVERAGE, DEPTH, Q_VAL, P_VAL]


def select_best_kma_reference(kma_file: str) -> str:
    with open(kma_file, "r") as inf:
        reader = DictReader(inf, delimiter = "\t", fieldnames = FIELDNAMES)
        rows = [row for row in reader]
        rows = rows[1:]

        if not rows:
            raise ValueError(f"KMA results file {kma_file} is empty!")

        rows.sort(key = lambda x: (float(x[SCORE]), float(x[TEMPLATE_COVERAGE]), float(x[TEMPLATE_IDENTITY]), float(x[DEPTH]), float(x[P_VAL])), reverse = True)
        best_ref = rows[0]
        print(f"best reference for file {kma_file}: {best_ref}")
        return best_ref[TEMPLATE] # return reference name

def find_ref_in_fasta(fasta: str, ref_name: str) -> str:
    found = False
    ret = ""

    with open(fasta, "r") as inf:
        for l in inf:
            if found:
                ret += l # add sequence
                return ret

            if l.startswith(">"):
                if ref_name in l:
                    found = True
                    ret += l # add header
                    continue

    raise ValueError(f"ref {ref_name} is not in the fasta file {fasta}!")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", required = True, help = "kma res file", type = str)
    parser.add_argument("--prefix", required = True, help = "prefix for output file", type = str)
    parser.add_argument("-r", required = True, help = "reference fasta that contains references listed in the KMA file", type = str)

    args = parser.parse_args()

    best_ref = select_best_kma_reference(args.i)
    ref_fasta = find_ref_in_fasta(args.r, best_ref)

    with open(f"{args.prefix}_best_ref.fa", "w") as outf:
        outf.write(ref_fasta)

