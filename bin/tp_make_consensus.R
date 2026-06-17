#!/usr/bin/env Rscript
# Adapted for Nextflow Aug 2020 for T. pallidum

# This script imports bam files and makes a consensus sequence
# Pavitra Roychoudhury
# Adapted from hsv_generate_consensus.R on 6-Mar-19
# Built to be called from wgs_pipeline.sh with input arguments specifying input filename
# Requires wgs_functions.R which contains several utility scripts plus multiple R packages listed below

rm(list = ls())
sessionInfo()
library(Rsamtools)
library(GenomicAlignments)
library(ShortRead)
library(Biostrings)
library(RCurl)
# Get args from command line
args <- (commandArgs(TRUE))
if (length(args) == 0) {
  print("No arguments supplied.")
} else {
  sampname <- args[[1]]
  bamfname <- args[[2]]
  ref <- args[[3]]
  min_depth <- as.integer(args[[4]])
}

# Files, directories, target site
merged_bam_folder <- "./"
mapped_reads_folder <- "./"
con_seqs_dir <- "./"

n_mapped_reads <- function(bamfname) {
  require(Rsamtools)
  indexBam(bamfname)
  if (file.exists(bamfname) & class(try(scanBamHeader(bamfname), silent = T)) != "try-error") {
    return(idxstatsBam(bamfname)$mapped)
  } else {
    return(NA)
  }
}

# Takes in a bam file, produces consensus sequence
generate_consensus <- function(bamfname, min_depth) {
  require(Rsamtools)
  require(GenomicAlignments)
  require(Biostrings)
  require(parallel)
  ncores <- detectCores()

  print(paste0("Using minimum consensus depth of ", min_depth))

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
    qseq_on_ref <- sequenceLayer(mcols(gal)$seq, cigar(gal), from = "query", to = "reference")
    cm <- consensusMatrix(qseq_on_ref, as.prob = F, shift = start(gal) - 1, width = seqlengths(gal))[c("A", "C", "G", "T", "N", "-"), ]
    poor_cov <- which(colSums(cm) < min_depth)
    cm <- apply(cm, 2, function(x) x / sum(x))
    cm[, poor_cov] <- 0
    cm["N", poor_cov] <- 1

    tmp_str <- strsplit(consensusString(cm, ambiguityMap = "?", threshold = 0.75), "")[[1]]
    ambig_sites <- which(tmp_str == "?")
    ambig_bases <- unlist(lapply(ambig_sites, function(i) {
      mixedbase <- paste(names(cm[, i])[cm[, i] > 0], collapse = "")
      if (mixedbase %in% IUPAC_CODE_MAP) {
        return(names(IUPAC_CODE_MAP)[IUPAC_CODE_MAP == mixedbase])
      } else {
        return("N")
      }
    }))
    tmp_str[ambig_sites] <- ambig_bases
    con_seq <- DNAStringSet(paste0(tmp_str, collapse = ""))
    names(con_seq) <- sub(".bam", "_consensus", basename(bamfname))
    rm(tmp_str)

    # Remove gaps and leading and trailing Ns to get final sequence
    con_seq_trimmed <- DNAStringSet(gsub("N*N$", "", gsub("^N*", "", as.character(con_seq))))
    con_seq_final <- DNAStringSet(gsub("-", "", as.character(con_seq_trimmed)))
    names(con_seq_final) <- sub(".bam", "_consensus", basename(bamfname))

    # Delete bai file
    file.remove(baifname)

    return(con_seq_final)
  } else {
  	print("Error generating consensus")
  	quit(status = 1)
  }
}

clean_consensus_tp <- function(sampname, bamfname, ref, min_depth) {
  # con_seqs <- lapply(bam, generate_consensus)
  consensus = generate_consensus(bamfname, min_depth)
  writeXStringSet(consensus, file = paste("./", sampname, "_consensus.fasta", sep = ""), format = "fasta")
}

# Make consensus sequence--returns TRUE if this worked
conseq <- clean_consensus_tp(sampname, bamfname, ref, min_depth)
