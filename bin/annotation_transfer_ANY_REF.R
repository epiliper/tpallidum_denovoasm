#!/usr/bin/env Rscript

library(tidyverse)
library(phylotools)
library(readxl)
library(Biostrings)

`%!in%` <- negate(`%in%`)

# dir.create("/Users/Nicole/Desktop/devira/transfer_annotations/test_output3/")

args <- commandArgs(TRUE)
alignment <- read.fasta(args[[1]])
ORFs <- read_excel(args[[2]])
target <- read.fasta(args[[3]])
outprefix <- args[[4]]

# read in alignment of final target sequence vs its best reference (NOT MASKED YET - or only masked in ARP repeat, 470 repeats, intra-rRNA trna)
# alignment <- read.fasta("/Users/Nicole/Desktop/devira/transfer_annotations/test_alignment1.fasta")
# ORFs <- read_excel(("/Users/Nicole/Desktop/devira/transfer_annotations/NC_021508_coordinates_FINAL.xlsx")) # %>% select(1:3)

# UPDATE AS APPROPRIATE FOR FINAL EXTENSION
# target <- read.fasta("/Users/Nicole/Desktop/devira/transfer_annotations/CP073506.1_denovo_polished_changes_Jul07.fasta")
# sample <- gsub("\\.1_denovo_polished_changes_Jul07", "", target[1,1])
sample <- outprefix
target[1, 1] <- sample


#######
split <- as.data.frame(strsplit(alignment$seq.text, split = ""))
colnames(split) <- sub(" .*", "", alignment$seq.name) # get everything before first space only
# colnames(split) <- c("ref", "target")

# which column in the alignment is the reference? Can't assume #1
all_refs <- c("CP007548", "NC_016842", "NC_016843", "NC_018722", "NC_021179", "NC_021490", "NC_021508", "NZ_CP032303")
# which reference was used?
ref_index <- which(names(split) %in% all_refs)
# Get column names
chosen_ref <- names(split)[ref_index]

# put ref in first position
split <- split %>% relocate(all_of(chosen_ref))
# relabel so uniform:
colnames(split) <- c("ref", "target")
split$gaps_ref <- ifelse(split$ref == "-", 0, 1)
split$ref_position <- cumsum(split$gaps_ref)

split$ref_position <- ifelse(split$gaps_ref == 0, NA, split$ref_position)

split$gaps_target <- ifelse(split$target == "-", 0, 1)
split$target_position <- cumsum(split$gaps_target)
split$target_position <- ifelse(split$gaps_target == 0, NA, split$target_position)

lookup_table_min <- data.frame(ref_min = split$ref_position, target_min = split$target_position)
colnames(lookup_table_min) <- c(paste0(chosen_ref, "_min"), "target_min")
lookup_table_max <- data.frame(ref_max = split$ref_position, target_max = split$target_position)
colnames(lookup_table_max) <- c(paste0(chosen_ref, "_max"), "target_max")

lookup_table_ref_min <- read_csv(paste0(chosen_ref, "_lookup.csv")) %>% select(1:2)
lookup_table_ref_max <- read_csv(paste0(chosen_ref, "_lookup.csv")) %>% select(3:4)

new_orfs <- ORFs %>%
  select(1:3, Direction) %>%
  left_join(lookup_table_ref_min) %>%
  left_join(lookup_table_ref_max) %>%
  left_join(lookup_table_min) %>%
  left_join(lookup_table_max)
new_orfs$chrom_sequence <- str_sub(target$seq.text, new_orfs$target_min, new_orfs$target_max)
sum(is.na(new_orfs$chrom_sequence))

new_orfs$chrom_sequence <- ifelse(is.na(new_orfs$chrom_sequence), "NNN", new_orfs$chrom_sequence)
new_orfs$gene_seq <- ifelse(new_orfs$Direction == "reverse", as.character(reverseComplement(DNAStringSet(new_orfs$chrom_sequence))), new_orfs$chrom_sequence)

bacterial_code <- getGeneticCode("11")
new_orfs$AA <- as.character(translate(DNAStringSet(new_orfs$gene_seq), genetic.code = bacterial_code, if.fuzzy.codon = "solve"))
new_orfs$target <- sample
new_orfs <- new_orfs %>% relocate(target)

write_csv(new_orfs, paste0(sample, "_extracted_ORFs.csv"))
