from glob import glob
from collections import defaultdict
from csv import DictWriter
import os

def get_mate_idx(filename: str):
    r1, r2 = filename.rfind("1"), filename.rfind("2")
    if r1 & r2: return 0
    return 1 if r1 else 2

def find_pairs(reads: set[str], sample_suffix: str) -> dict[str, list[str]]:
    matemap: dict[str, list[str]] = defaultdict(list[str])
    for r in reads:
        prefix = os.path.basename(r)
        prefix = prefix.split(sample_suffix)[0].strip(" _")
        assert(prefix)

        matemap[prefix].append(r)

    for v in matemap.values():
        v.sort(key = lambda x: get_mate_idx(x))
        assert(len(v) <= 2)

    return matemap

def create_samplesheet():
    trna_files = set([os.path.join(args.dir, f) for f in glob(f"**/*{args.tprefix}*[fastq|fq]*", root_dir = args.dir, recursive = True)])
    rrna_files = set([os.path.join(args.dir, f) for f in glob(f"**/*{args.rprefix}*[fastq|fq]*", root_dir = args.dir, recursive = True)])
    gna_files = set([os.path.join(args.dir, f) for f in glob(f"**/*{args.gprefix}*[fastq|fq]*", root_dir = args.dir, recursive = True)])

    intersect = set.intersection(trna_files, rrna_files, gna_files)

    if len(intersect) > 0:
        raise ValueError(f"patterns for nucleotide types ({args.tprefix} {args.rprefix} {args.gprefix}) are not distinctive: these files are associated with more than one nucleic acid type: {intersect}")

    trna_files = find_pairs(trna_files, args.tprefix)
    rrna_files = find_pairs(rrna_files, args.rprefix)
    gna_files = find_pairs(gna_files, args.gprefix)

    sample_df: dict[str, dict[str, str]] = defaultdict(dict)

    for sample, values in trna_files.items():
        r1 = values[0]
        try:
            r2 = values[1] 
        except IndexError:
            r2 = ""

        sample_df[sample]["trna_fastq_1"] = r1
        sample_df[sample]["trna_fastq_2"] = r2

    for sample, values in rrna_files.items():
        r1 = values[0]
        try:
            r2 = values[1] 
        except IndexError:
            r2 = ""

        sample_df[sample]["rrna_fastq_1"] = r1
        sample_df[sample]["rrna_fastq_2"] = r2

    for sample, values in gna_files.items():
        r1 = values[0]
        try:
            r2 = values[1] 
        except IndexError:
            r2 = ""

        sample_df[sample]["gna_fastq_1"] = r1
        sample_df[sample]["gna_fastq_2"] = r2

    fields = ["sample", "gna_fastq_1", "gna_fastq_2", "rrna_fastq_1", "rrna_fastq_2", "trna_fastq_1", "trna_fastq_2"]

    # remove any samples missing a nucleic acid type
    to_del = []
    for k, v in sample_df.items():
        if "gna_fastq_1" not in v:
            to_del.append(k)

        elif "rrna_fastq_1" not in v:
            to_del.append(k)

        elif "trna_fastq_1" not in v:
            to_del.append(k)

    if to_del:
        print(f"found {len(to_del)} samples with missing fastqs for at least one nucleotide type. These will be removed from the final samplesheet.")

    for k in to_del:
        print(f"removing {k}")
        del sample_df[k]

    if not sample_df:
        raise ValueError("All samples identified had at least one nucleic acid type FASTQ missing; the samplesheet is empty.")


    with open(args.output, "w") as outf:
        writer = DictWriter(outf, fields)
        writer.writeheader()
        for k in sample_df:
            writer.writerow({field: sample_df[k].get(field) or k for field in fields})

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("dir", help = "directory that contains all subdirs for nucleic acid types")

    parser.add_argument("--tprefix", required = True, help = "prefix that distinguishes tRNA filenames", type = str)
    parser.add_argument("--rprefix", required = True, help = "prefix that distinguishes rRNA filenames")
    parser.add_argument("--gprefix", required = True, help = "prefix that distinguishes gNA filenames")
    parser.add_argument("--output", help = "name of output samplesheet", default = "samplesheet.csv")

    args = parser.parse_args()

    create_samplesheet()

