#!/usr/bin/env python3

from argparse import ArgumentParser
from csv import DictReader

# column name defines
TEMPLATE = "# Template"
SCORE = "Score"
EXPECTED = "Expected"
TEMPLATE_LENGTH = "Template_length"
TEMPLATE_IDENTITY = "Template_identity"
TEMPLATE_COVERAGE = "Template_Coverage"
QUERY_IDENTITY = "Query_Identity"
QUERY_COVERAGE = "Query_Coverage"
DEPTH = "Depth"
Q_VAL = "q_value"
P_VAL = "p-value"


def select_best_kma_reference(kma_file: str) -> str:
    with open(kma_file, "w") as inf:
        reader = DictReader(inf, delimiter = "\t")
        rows = [row for row in reader]

        if not rows:
            raise ValueError(f"KMA results file {kma_file} is empty!")

        rows.sort(key = lambda x: (x[SCORE], x[TEMPLATE_COVERAGE], x[TEMPLATE_IDENTITY], x[DEPTH], x[P_VAL]), reverse = True)
        best_ref = rows[0]
        print(f"best reference for file {kma_file}: {best_ref}")
        return best_ref[TEMPLATE] # return reference name

def find_ref_in_fasta(fasta: str, ref_name: str) -> str:
    found = False
    ret = ""

    with open(fasta, "w") as inf:
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

    with open(f"{args.prefix}_best_ref.fa") as outf:
        outf.write(ref_fasta)

