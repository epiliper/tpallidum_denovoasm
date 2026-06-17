#!/usr/bin/env python3
import subprocess
import argparse
import sys
import os
import regex
from itertools import chain
import numpy as np
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
import pandas as pd

#This script is ONLY to mask for upload to NCBI.
#intra-rRNA tRNA, arp repeats, tp0470 repeats, and tprk V regions are masked

# Matches a read to a specified string of nucleotides, "primer".
def fuzzy_match(read_seq, primer,num_mismatches):
	# Finds exact match first.
	exact_match = regex.search(primer,read_seq)
	if exact_match:
		return exact_match[0].rstrip()
	# If can't find exact match, searches for best match with less than specified number of mismatches.
	else:
		if(num_mismatches==1):
			fuzzy_match = regex.search(r"(?b)("+primer + "){s<=1}", read_seq)
		else:
			fuzzy_match = regex.search(r"(?b)("+primer + "){s<=3}", read_seq)
		if fuzzy_match:
			return fuzzy_match.group()
		else:
			return 0
			#if(len(fuzzy_match)>1):
				#print("!!!!!!!!! Multiple matches found!!!!!!!!!")
				#print(fuzzy_match.group())
				#print(fuzzy_match[0])
				#print(fuzzy_match[1])
				#return 0
			#return fuzzy_match[0]

# Masks repetitive/variable regions based on upstream and downstream elements.
def mask_variable_region(fasta_seq,upstream,downstream,gene_name,num_mismatches):
	# Finds the upstream portion of the gene.
	gene_upstream = fuzzy_match(fasta_seq,upstream,num_mismatches)

	#print(upstream)
	#print(fasta_seq)
	#print(gene_upstream)
	#print(num_mismatches)

	if (gene_upstream == 0):
		#print(gene_name,"not found.")
		return fasta_seq

	# Finds the start position of the gene of interest.
	gene_start = str.index(fasta_seq, gene_upstream) + len(gene_upstream)

	# Searches for the downstream portion of the gene starting from the start of the gene to a certain number of bases away, depending on the gene,
	# to further limit multiple matches.
	if(gene_name == "arp gene"):
		fuzzy_gene_end = gene_start + 3000
	elif(gene_name == "tp0470 gene"):
		fuzzy_gene_end = gene_start + 3000
	else:
		fuzzy_gene_end = gene_start + 10000
	gene_downstream = fuzzy_match(fasta_seq[gene_start:fuzzy_gene_end],downstream,num_mismatches)

	if (gene_downstream == 0):
		#print(gene_name,"not found.")
		return fasta_seq

	# Finds the end position of the gene
	gene_end = str.index(fasta_seq,gene_downstream)

	# Grabs the sequence between gene start and end positions
	gene_seq = fasta_seq[gene_start:gene_end]

	#print(gene_name,"found.")#: ",gene_seq,". Masking...")

	# Replace the found genes with the same length of Ns
	return fasta_seq.replace(gene_seq,'N'*len(gene_seq))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Masking variable regions in for Tpallidum WGS.')
	parser.add_argument('-f', '--fasta', help='Fasta file to mask variable regions for.')

	# Checks for argument sanity.
	try:
		args = parser.parse_args()
	except:
		parser.print_help()
		sys.exit(1)

	fasta_file = args.fasta

	masked_sequence = ""

	sample_name = fasta_file.split(".fasta")[0]
	#print("Masking genes in ",fasta_file,"...")
	masked_fasta_file = open(sample_name+"_masked_NCBI.fasta","w+")

	for line_num, line in enumerate(open(fasta_file)):
		# Writes fasta header
		if(line_num==0):
			masked_fasta_file.write(line)
		# Reads fasta sequence and masks variable/repetitive regions
		else:

			masked_sequence = masked_sequence + line

	masked_sequence = masked_sequence.replace("\n","")

	# intra-rrna trna1 gene
	seq_to_replace = masked_sequence[230000:235000]
	replaced_gene = mask_variable_region(seq_to_replace,"TCTCCCCTTCCCTTTTGAAAA","CTATTATTCTTTATGTCCCTT","intra-rrna trna1",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	# intra-rrna trna2 gene
	seq_to_replace = masked_sequence[279000:285000]
	replaced_gene = mask_variable_region(seq_to_replace,"TCTCCCCTTCCCTTTTGAAAA","CTATTATTCTTTATGTCCCTT","intra-rrna trna2",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	# arp gene
	seq_to_replace = masked_sequence[457255:469137]
	replaced_gene = mask_variable_region(seq_to_replace,"TTTGGTTTCCCCTTTGTCTC","AGGTCGCTTCTCAGCATACG","arp gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	##tp0470 gene
	seq_to_replace = masked_sequence[493859:504150]
	replaced_gene = mask_variable_region(seq_to_replace,"GCGCGCTTGAGAGCTTCAAA","GCGCTGCAGCCACTGCTCAA","tp0470 gene",1)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	# Masking tprK variable regions!
	# Some seem to require more than 1 mismatch.
	seq_to_replace = masked_sequence[972000:984000]
			
	#replaced_gene = mask_variable_region(seq_to_replace,"ATCAGTAGTAGTCTTAAATCC","CCAGGCCAGCTCCGCATA","tprK V1",2)
	replaced_gene = mask_variable_region(seq_to_replace,"TCAGTAGTAGTCTTAAATCC","CCCAGGCCAGCTCCGCATAG","tprK V1",2)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)
			
	seq_to_replace = masked_sequence[972000:984000]
	#replaced_gene = mask_variable_region(seq_to_replace,"AATATCTCCCCCCAATCCATA","GTCGGTGTTAGACGCAAA","tprK V2",2)
	replaced_gene = mask_variable_region(seq_to_replace,"AGCCGAACAAAATATCTCCC","GTCGGTGTTAGACGCAAACT","tprK V2",2)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

	seq_to_replace = masked_sequence[972000:984000]
	#replaced_gene = mask_variable_region(seq_to_replace,"TCATACTCACCTTAGCCCCGAC","GGAGTTGCCGGTGAGCTC","tprK V3",3)
	replaced_gene = mask_variable_region(seq_to_replace,"GCACACAGACCCCAAAGCTT","GTGGAGTTGCCGGTGAGCTC","tprK V3",3)

	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)
			
	seq_to_replace = masked_sequence[972000:984000]
	#replaced_gene = mask_variable_region(seq_to_replace,"AACAACGCATCTGCGCC","AGCAGCCAGAGCACACA","tprK V4",2)
	replaced_gene = mask_variable_region(seq_to_replace,"ACCCCAACGTCAACAACGCA","TAGCAGCCAGAGCACACAGA","tprK V4",2)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)
			
	seq_to_replace = masked_sequence[972000:984000]
	#replaced_gene = mask_variable_region(seq_to_replace,"TCGAGCTTAATATAGGCAGC","CGATGCGAAATATCCTCC","tprK V5",2)
	replaced_gene = mask_variable_region(seq_to_replace,"GACCCCTTGGTTTCGAGCTT","TCCTCCCGCCGAGAACCAAC","tprK V5",2)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)
			
	seq_to_replace = masked_sequence[972000:984000]
	#replaced_gene = mask_variable_region(seq_to_replace,"TTCCATACACCGGGAA","CATGTACGTACGCACATC","tprK V6",2)
	replaced_gene = mask_variable_region(seq_to_replace,"CCCCAGACTTTTCCATACAC","CATGTACGTACGCACATCAA","tprK V6",2)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)
			
	seq_to_replace = masked_sequence[972000:984000]
	#replaced_gene = mask_variable_region(seq_to_replace,"CGGACTGACCACTACCCCACACTC","CAAGTTTGCATACACTTT","tprK V7",2)
	replaced_gene = mask_variable_region(seq_to_replace,"ACTGACCACTACCCCACACT","CAAGTTTGCATACACTTTAA","tprK V7",2)
	masked_sequence = masked_sequence.replace(seq_to_replace,replaced_gene)

    # replace ambiguous bases with N
	nonambig = set(["A", "G", "C", "T", "N"])
	filter = lambda b: b if b.upper() in nonambig else "N" 
	masked_sequence = "".join([filter(b) for b in masked_sequence])


	masked_fasta_file.write(masked_sequence)
