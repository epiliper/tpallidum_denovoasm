#!/usr/bin/env Rscript

# Adapted for Nextflow Aug 2020 for T. pallidum

# HSV script but works with any viral sequence: This script makes a new reference sequence from de novo assembled scaffolds
# Pavitra Roychoudhury
# Sep 2017

# Built to be called from hsv_wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

# 2026-May-1: EP removed the code that generates the BAM; an upstream step does so already

rm(list = ls())
sessionInfo()
library(Rsamtools)
library(GenomicAlignments)
library(Biostrings)
library(RCurl)
library(parallel)
# Get args from command line
args <- (commandArgs(TRUE))
if (length(args) == 0) {
    print("No arguments supplied.")
} else {
		sampname <- args[[1]]
    bamname <- args[[2]]
    reffname <- args[[3]]
    min_contig_len <- as.numeric(args[[4]])
}

con_seq_final = ""

# note: ported from another python script by EP
report_contig_stats <- function(query_seqs, num_fastas, scaf_len, output_file) {
  #' Reports N50, number of total and filtered alignments, and scaffold length to a TSV file.
  #'
  #' @param query_seqs  A GAlignments object
  #' @param num_fastas  Number of seqs prior to any filtering (integer)
  #' @param scaf_len    Length of scaffold assembled with the seqs (integer)
  #' @param output_file Path to the output TSV file (string)
  #' @return NULL (invisibly)

  if (length(query_seqs) == 0) {
    longest <- num_filtered <- shortest <- n50 <- scaf_len <- -1L
  } else {
    lengths <- sort(GenomicAlignments::qwidth(query_seqs), decreasing = TRUE)

    longest      <- lengths[1]
    shortest     <- lengths[length(lengths)]
    num_filtered <- length(lengths)

    get_n50 <- function(lengths) {
      #' Compute N50 from a vector of contig lengths sorted in descending order
      total_length <- sum(lengths)
      cumu_length  <- 0L
      for (l in lengths) {
        cumu_length <- cumu_length + l
        if (cumu_length >= total_length / 2) return(l)
      }
      return(-1L)
    }

    n50 <- get_n50(lengths)
    message(sprintf("Writing contig stats to report file: %s", output_file))
  }

  header <- paste(
    "num_total_contigs", "num_filtered_contigs", "N50",
    "shortest_contig_len", "longest_contig_len", "scaffold_length",
    sep = "\t"
  )
  values <- paste(num_fastas, num_filtered, n50, shortest, longest, scaf_len, sep = "\t")

  writeLines(c(header, values), con = output_file)
  invisible(NULL)
}

# Make a new reference from scaffolds
make_ref_from_assembly <- function(bamfname, reffname) {
    require(Rsamtools)
    require(GenomicAlignments)
    require(parallel)
    ncores <- detectCores()

    # Read reference sequence
    ref_seq <- readDNAStringSet(reffname)

    if (!is.na(bamfname) & class(try(scanBamHeader(bamfname), silent = T)) != "try-error") {
        # Index bam if required
        if (!file.exists(paste(bamfname, ".bai", sep = ""))) {
            baifname <- indexBam(bamfname)
        } else {
            baifname <- paste(bamfname, ".bai", sep = "")
        }

        # Import bam file
        params <- ScanBamParam(
            flag = scanBamFlag(isUnmappedQuery = FALSE),
            what = c("qname", "rname", "strand", "pos", "qwidth", "mapq", "cigar", "seq")
        )
        gal <- readGAlignments(bamfname, index = baifname, param = params)

        # Remove any contigs with width less than specified in args
        original_len = length(gal)

        gal <- gal[width(gal) > min_contig_len]

        # First lay contigs on reference space--this removes insertions and produces a seq of the same length as ref
        qseq_on_ref <- sequenceLayer(mcols(gal)$seq, cigar(gal), from = "query", to = "reference")
        qseq_on_ref_aligned <- stackStrings(qseq_on_ref, 1, max(mcols(gal)$pos + qwidth(gal) - 1, width(ref_seq)),
            shift = mcols(gal)$pos - 1, Lpadding.letter = "N", Rpadding.letter = "N"
        )

        # Make a consensus matrix and get a consensus sequence from the aligned scaffolds
        cm <- consensusMatrix(qseq_on_ref_aligned, as.prob = T, shift = 0)[c("A", "C", "G", "T", "N", "-"), ]
        # cm[c('N','-'),]<-0;
        cm["N", ] <- 0
        cm <- apply(cm, 2, function(x) if (all(x == 0)) {
            return(x)
        } else {
            return(x / sum(x))
        })
        cm["N", colSums(cm) == 0] <- 1
        con_seq <- DNAStringSet(gsub("\\?", "N", consensusString(cm, threshold = 0.25)))
        con_seq <- DNAStringSet(gsub("\\+", "N", con_seq))


        # Now fill in the Ns with the reference
        temp <- as.matrix(con_seq)
        temp[temp == "N"] <- as.matrix(ref_seq)[temp == "N"]
        con_seq <- DNAStringSet(paste0(temp, collapse = ""))
        names(con_seq) <- sub(".bam", "_consensus", basename(bamfname))

        # Look for insertions in bam cigar string
        cigs_ref <- cigarRangesAlongReferenceSpace(cigar(gal),
            with.ops = F, ops = "I",
            reduce.ranges = T, drop.empty.ranges = F,
            pos = mcols(gal)$pos
        )
        cigs_query <- cigarRangesAlongQuerySpace(cigar(gal),
            ops = "I", with.ops = F,
            reduce.ranges = T, drop.empty.ranges = F
        )
        all_ins <- mclapply(c(1:length(cigs_query)), function(i) {
            extractAt(mcols(gal)$seq[i], cigs_query[[i]])[[1]]
        })

        # Merge all insertions
        all_ins_merged <- do.call("rbind", mclapply(c(1:length(cigs_ref)), function(i) {
            return(data.frame(
                start_ref = start(cigs_ref[[i]]), end_ref = end(cigs_ref[[i]]),
                start_query = start(cigs_query[[i]]), end_query = end(cigs_query[[i]]),
                ins_seq = all_ins[[i]], width_ins = width(all_ins[[i]])
            ))
        },
        mc.cores = ncores
        ))
        all_ins_merged <- all_ins_merged[order(all_ins_merged$end_ref), ]

        # write.csv(all_ins_merged,'./testing/all_ins.csv',row.names=F);

        # TO DO: Check for overlaps--should be minimal since scaffolds don't usually overlap that much
        if (any(table(all_ins_merged$start_ref) > 1)) {
            print("Overlapping insertions")
            # not the best way, but just pick the first for now
            all_ins_merged <- all_ins_merged[!duplicated(all_ins_merged[, c("start_ref", "end_ref")]), ]
        }

        # Now the beauty part of inserting the strings back in
        # Split ref seq by the insert positions
        if (nrow(all_ins_merged) != 0) {
            new_strs <- DNAStringSet(rep("", nrow(all_ins_merged) + 1))
            for (i in 1:nrow(all_ins_merged)) {
                if (i == 1) {
                    new_strs[i] <- paste0(
                        extractAt(con_seq, IRanges(start = 1, end = all_ins_merged$end_ref[i]))[[1]],
                        all_ins_merged$ins_seq[i]
                    )
                } else {
                    new_strs[i] <- paste0(
                        extractAt(con_seq, IRanges(
                            start = all_ins_merged$start_ref[i - 1],
                            end = all_ins_merged$end_ref[i]
                        ))[[1]],
                        all_ins_merged$ins_seq[i]
                    )
                }
            }

            # Last bit
            new_strs[i + 1] <- paste0(extractAt(con_seq, IRanges(
                start = all_ins_merged$start_ref[i],
                end = width(con_seq)
            ))[[1]])
            temp_str <- paste0(as.character(new_strs), collapse = "")

            # Remove gaps to get final sequence
            con_seq_final <- DNAStringSet(gsub("-", "", temp_str))

            # No insertions
        } else {
            con_seq_final <- con_seq
        }

        names(con_seq_final) <- sub(".bam", "_consensus", basename(bamfname))

        writeXStringSet(
            con_seq_final,
            paste0("./", sampname, "_initial_consensus.fasta")
        )


        # Delete bai file
        file.remove(baifname)
    } else {
        print("Bam file could not be opened.")
        return(NA)
    }

		report_contig_stats(gal, original_len, length(con_seq_final), paste0(sampname, "_contig_stats.tsv"))
		return (names(con_seq_final))
}


# Make a new reference scaffold
newref <- make_ref_from_assembly(bamname, reffname)

if (is.na(newref)) print("Failed to generate consensus from scaffolds")
