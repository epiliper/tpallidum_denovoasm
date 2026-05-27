#!/usr/bin/env python3

import subprocess

def polish(query: str, ref: str, prefix:str, output_fasta: str):
    qheader = seq = ""

    refheader = ""

    if not "denovo" in prefix:
        prefix += "_denovo"

    if not prefix.startswith(">"):
        prefix = ">" + prefix

    with open(query) as qif:
        qheader = qif.readline().strip();
        seq = qif.readline().strip();

        seq = seq.lstrip("Nn")
        seq = seq.rstrip("Nn")

    assert(qheader[0] == ">"), f"FASTA header of query file {query} is invalid"
    assert(seq), f"FASTA file {query} doesn't seem to have sequence after leading/trailing N removal, invalid format?"

    mafft_in = "mafft_in"
    mafft_out = prefix.replace(">", "") + "_aligned.fa"

    # create alignment input
    with open(mafft_in, "w") as alnin:

        _ = alnin.write(f"{prefix}\n")
        _ = alnin.write(f"{seq}\n")

        with open(ref) as inref:
            for line in inref:
                _line = line.strip()
                if _line.startswith(">"):
                    assert(refheader == ""), f"reference file {ref} has multiple records. Not allowed."
                    refheader = _line

                _ = alnin.write(f"{_line}\n")

    assert(refheader), f"no fasta headers detected in {ref}"



    # align with mafft
    mafft_cmd = ["mafft", "--reorder", mafft_in]

    with open(mafft_out, "w") as out:
        _ = subprocess.check_call(mafft_cmd, stdout = out)

    seqmap = {prefix: "", refheader: ""}

    curref = ""

    # mafft has a little "quirk :)" where it can split records across multiple lines, so we need to concat
    with open(mafft_out) as aln:
        for line in aln:
            if line.startswith(">"):
                s = line.strip()
                assert(s in seqmap), f"unknown ref: {line}"
                curref = s

            else:
                seqmap[curref] += line.strip()

    qseq, refseq = list(seqmap[prefix]), list(seqmap[refheader])

    # replace missing sequence ends of query relative to reference
    for i in range(len(qseq)):
        if qseq[i] != "-":
            break;
        else:
            qseq[i] = "N"

    for i in range(len(qseq) - 1, -1, -1):
        if qseq[i] != "-":
            break;
        else:
            qseq[i] = "N"

    # trim overhanging query ends to ends of reference.
    loverhang = roverhang = 0

    for i in range(len(refseq)):
        if refseq[i] == "-":
            loverhang += 1
        else: 
            break

    for i in range(len(refseq) - 1, -1, -1):
        if refseq[i] == "-":
            roverhang += 1
        else: 
            break

    qseq = qseq[loverhang:]
    qseq = qseq[:len(qseq) - roverhang]
    refseq = refseq[loverhang:]
    refseq = refseq[:len(refseq) - roverhang]

    def convert(c: str) -> str:
        if c == "n" or c == "-":
            return "N"
        return c

    # replace interior gaps or n with N
    qseq = [convert(c) for c in qseq]

    # sanity check
    assert(len(qseq) == len(refseq)), f"length mismatch, polishing failed (this is an internal error): query = {len(qseq)}, ref = {len(refseq)}"

    # write 
    with open(output_fasta, "w") as outf:
        _ = outf.write(f"{prefix}_polished\n")
        _ = outf.write("".join(qseq))

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    _ = parser.add_argument("--query", type = str, required = True, help = "query fasta to polish")
    _ = parser.add_argument("--reference", type = str, required = True, help = "reference fasta")
    _ = parser.add_argument("--header_name", type = str, required = True, help = "text to rename header")
    _ = parser.add_argument("--output_fasta", type = str, required = True, help = "output fasta")

    args = parser.parse_args()

    polish(args.query, args.reference, args.header_name, args.output_fasta)
