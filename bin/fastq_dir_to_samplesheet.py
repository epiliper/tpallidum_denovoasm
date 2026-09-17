from glob import glob
from collections import defaultdict
from csv import DictWriter
import os
import re

FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


def is_fastq(filename: str) -> bool:
    return filename.lower().endswith(FASTQ_SUFFIXES)

def find_fastqs(directory: str, sample_prefix: str) -> set[str]:
    if not directory.startswith("s3://"):
        return {
            os.path.join(directory, filename)
            for filename in glob(
                f"**/*{sample_prefix}*", root_dir=directory, recursive=True
            )
            if is_fastq(filename)
        }

    try:
        import boto3
    except ImportError as error:
        raise ImportError("boto3 is required when the input directory is an S3 URI") from error

    bucket, _, prefix = directory[5:].partition("/")
    if not bucket:
        raise ValueError(f"Invalid S3 directory: {directory}")

    prefix = prefix.strip("/")
    if prefix:
        prefix += "/"

    files = set()
    paginator = boto3.client("s3").get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for item in page.get("Contents", []):
            key = item["Key"]
            filename = key.rsplit("/", 1)[-1]
            if sample_prefix in filename and is_fastq(filename):
                files.add(f"s3://{bucket}/{key}")

    return files

def get_mate_idx(filename: str):
    filename = os.path.basename(filename)
    matches = re.findall(r"(?:R|^|[._-])([12])(?=[._-]|$)", filename, re.IGNORECASE)
    return int(matches[-1]) if matches else 0

def find_pairs(reads: set[str], sample_suffix: str) -> dict[str, list[str]]:
    matemap: dict[str, list[str]] = defaultdict(list[str])
    for r in reads:
        prefix = os.path.basename(r)
        prefix = prefix.split(sample_suffix)[0].strip(" _")
        assert(prefix)

        matemap[prefix].append(r)

    for v in matemap.values():
        v.sort(key = lambda x: get_mate_idx(x))
        assert len(v) <= 2, f"Invalid mate bucket: {v}"

    return matemap

def create_samplesheet():
    fields = ["sample", "gna_fastq_1", "gna_fastq_2", "rrna_fastq_1", "rrna_fastq_2", "trna_fastq_1", "trna_fastq_2"]
    sample_rows = []
    missing_samples = []

    for directory in args.dir:
        trna_files = find_fastqs(directory, args.tprefix)
        rrna_files = find_fastqs(directory, args.rprefix)
        gna_files = find_fastqs(directory, args.gprefix)
        intersect = set.intersection(trna_files, rrna_files, gna_files)

        if intersect:
            raise ValueError(f"patterns for nucleotide types ({args.tprefix} {args.rprefix} {args.gprefix}) are not distinctive in {directory}: these files are associated with more than one nucleic acid type: {intersect}")

        pairs = (
            ("trna", find_pairs(trna_files, args.tprefix)),
            ("rrna", find_pairs(rrna_files, args.rprefix)),
            ("gna", find_pairs(gna_files, args.gprefix)),
        )
        samples: dict[str, dict[str, str]] = defaultdict(dict)
        for nucleic_type, paired_files in pairs:
            for sample, values in paired_files.items():
                samples[sample][f"{nucleic_type}_fastq_1"] = values[0]
                samples[sample][f"{nucleic_type}_fastq_2"] = (
                    values[1] if len(values) > 1 else ""
                )

        for sample, values in samples.items():
            if not all(
                f"{nucleic_type}_fastq_1" in values
                for nucleic_type in ("gna", "rrna", "trna")
            ):
                missing_samples.append((directory, sample))
                continue
            sample_rows.append({"sample": sample, **values})

    if missing_samples:
        print(f"found {len(missing_samples)} samples with missing fastqs for at least one nucleotide type. These will be removed from the final samplesheet.")
        for directory, sample in missing_samples:
            print(f"removing {sample} from {directory}")

    if not sample_rows:
        raise ValueError("All samples identified had at least one nucleic acid type FASTQ missing; the samplesheet is empty.")


    with open(args.output, "w") as outf:
        writer = DictWriter(outf, fields)
        writer.writeheader()
        writer.writerows(sample_rows)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("dir", nargs="+", help = "local directories or S3 URIs that contain subdirs for nucleic acid types")

    parser.add_argument("--tprefix", required = True, help = "prefix that distinguishes tRNA filenames", type = str)
    parser.add_argument("--rprefix", required = True, help = "prefix that distinguishes rRNA filenames")
    parser.add_argument("--gprefix", required = True, help = "prefix that distinguishes gNA filenames")
    parser.add_argument("--output", help = "name of output samplesheet", default = "samplesheet.csv")

    args = parser.parse_args()

    create_samplesheet()
