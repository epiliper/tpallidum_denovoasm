#!/usr/bin/env python3 

from enum import StrEnum
from csv import DictReader, DictWriter
from Bio import SeqIO, SeqRecord
import pysam
from typing import Any

# The columns of the instruction sheet we expect. These are meant to be read into a MaskInstruction class. See mask_instruction.csv in assets/ for an example.
class MaskInstructionSheetColumn(StrEnum):
    REGION = "region"
    """Locus name"""

    LINKAGE = "linkage_required"
    """should be an interval (in genomic coordinates, 1-based) delimited by a dash (-): e.g. 1-5"""

    MASK_IF_CONTIGUOUS = "mask_if_contiguous"
    """genomic coordinate interval"""

    MASK_IF_NOT_CONTIGUOUS = "mask_if_not_contiguous"
    """genomic coordinate interval"""

    MASK_MARGIN = "mask_outside"
    """should be Y/N or y/n"""


class MaskInstruction:
    region: str
    linkage_required: tuple[int, int] | None
    mask_if_contiguous: tuple[int, int] | None
    mask_if_not_contiguous: tuple[int, int] | None

    def __init__(self, name: str, linkage_required: tuple[int, int] | None, mask_if_contiguous: tuple[int, int] | None, mask_if_not_contiguous: tuple[int, int] | None):
        self.region = name
        self.linkage_required = linkage_required
        self.mask_if_contiguous = mask_if_contiguous
        self.mask_if_not_contiguous = mask_if_not_contiguous

# convert genomic interval to zero-based
def parse_numeric_interval(val: str) -> tuple[int, int] | None:
    val = val.strip()
    if val == "-": return

    split = val.split("-")
    if len(split) != 2: raise ValueError(f"Invalid coordinate range format: {val}")

    try:
        start = max(int(split[0]) - 1, 0)
    except ValueError:
        raise ValueError("Invalid coordinate start: {}", split[0])

    try:
        end = max(int(split[1]) - 1, 0)
    except ValueError:
        raise ValueError("Invalid coordinate start: {}", split[1])

    assert start < end, f"Invalid numeric interval {start}-{end}, 'start' must be < 'end'"

    return (start, end)


def mask_region(input: list[str], start: int, end: int, region_name: str) -> list[str]:
    assert start < end, f"bad input to masking: {start} should be < than {end}"
    print(f"Masking region {start} - {end} ({region_name})")
    return input[: start] + (["N"] * (end - start + 1)) + input[end + 1 :]


def parse_instruction_sheet(path: str) -> list[MaskInstruction]:
    ret: list[MaskInstruction] = []
    if not path.endswith(".csv"):
        raise ValueError("input path should be a CSV file.")


    with open(path, "r") as inf:
        reader = DictReader(inf)

        # sanity check stuff
        for row in reader:
            region = row[MaskInstructionSheetColumn.REGION].strip()
            linkage_required = parse_numeric_interval(row[MaskInstructionSheetColumn.LINKAGE])
            mask_if_contiguous = parse_numeric_interval(row[MaskInstructionSheetColumn.MASK_IF_CONTIGUOUS])
            mask_if_not_contiguous = parse_numeric_interval(row[MaskInstructionSheetColumn.MASK_IF_NOT_CONTIGUOUS])

            ret.append(MaskInstruction(region, linkage_required, mask_if_contiguous, mask_if_not_contiguous))

    return ret

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser()
    _ = parser.add_argument("--instructions", required = True, type = str, help = "sheet with masking instructions")
    _ = parser.add_argument("--bam", required = True, type = str, help = "bam alignment of contigs against reference")
    _ = parser.add_argument("--fasta", required = True, type = str, help = "input fasta for masking")
    _ = parser.add_argument("--output_prefix", required = True, type = str, help = "output prefix")
    args = parser.parse_args()

    instructions = parse_instruction_sheet(args.instructions)

    mode = "rb" if args.bam.endswith(".bam") else "r"
    bamreader = pysam.AlignmentFile(args.bam, mode)
    assert len(bamreader.references) == 1, "input BAM file must have only one reference."

    fastaref = next(SeqIO.parse(args.fasta, "fasta"))
    fastaseq = list(str(fastaref.seq))
    bamref = bamreader.get_reference_name(0)

    reports: list[dict[Any, Any]] = []

    for ins in instructions:
        assert ins.mask_if_not_contiguous, f"loci should have an interval specified in column '{MaskInstructionSheetColumn.MASK_IF_NOT_CONTIGUOUS}'"

        start, end = ins.linkage_required if ins.linkage_required else ins.mask_if_not_contiguous
        records = bamreader.fetch(bamref, start, end)

        # if we don't specify linkage required, we assume that we mask the entire interval specified in mask_if_not_contiguous
        if not ins.mask_if_contiguous and not ins.linkage_required:
            start, end = ins.mask_if_not_contiguous
            fastaseq = mask_region(fastaseq, start, end, ins.region)
            reports.append({"region": ins.region, "mask": True, "coordinates": f"{start} - {end}"})
            continue

        # otherwise, see if we have a contig covering the entire region
        assert ins.linkage_required, "linkage_required field missing"
        start, end = ins.linkage_required

        linkage_found = False
        for record in records:
            assert record.reference_end, f"Record {record.query_name} does not have a reference end, anomalous data?"

            if record.is_unmapped or record.is_supplementary or record.is_secondary: continue

            if (record.reference_start <= start) and (record.reference_end >= end):
                linkage_found = True
                break;

        # no single contig found for this region
        if not linkage_found:
            print(f"Masking region {start} - {end} ({ins.region})")
            start, end = ins.mask_if_not_contiguous
            fastaseq = mask_region(fastaseq, start, end, ins.region)
            reports.append({"region": ins.region, "mask": True, "coordinates": f"{start} - {end}"})

        elif ins.mask_if_contiguous: 
            start, end  = ins.mask_if_contiguous
            fastaseq = mask_region(fastaseq, start, end, ins.region)
            reports.append({"region": ins.region, "mask": True, "coordinates": f"{start} - {end}"})

        else:
            reports.append({"region": ins.region, "mask": False, "coordinates": f"{start} - {end}"})

    bamreader.close()

    assert len(fastaseq) == len(str(fastaref.seq)),  f"length mismatch: {len(fastaseq)} vs {len(fastaref.seq)}"

    # write modified fasta
    with open(f"{args.output_prefix}_genemasked.fasta", "w") as outf:
        _ = outf.write(f">{fastaref.id}_genemasked\n{"".join(fastaseq)}\n")

    # write mask report
    with open(f"{args.output_prefix}_genemask.tsv", "w") as outf:
        writer = DictWriter(outf, fieldnames = reports[0].keys(), delimiter = "\t")
        writer.writeheader()

        for row in reports:
            writer.writerow(row)
