#!/bin/python

import collections
from collections import Counter
from sys import argv 
script, filename, filename2 = argv 
import re

p = open(filename) #bed formatted coordinate file
g = open(filename2) #celegans genome (fa)

def revcom(DNA):
	DNAb = list(DNA)
	DNAr = DNAb[::-1]
	DNArc = []
	for base in DNAr:
		if base == 'A':
			DNArc.append('T')
		elif base == 'C':
			DNArc.append('G')
		elif base == 'G':
			DNArc.append('C')
		elif base == 'T':
			DNArc.append('A')
		else:
			print "non-DNA character in string!"
	return(DNArc)


def fastaparse(file):
	fstring = ''
	for line in file:
		line = line.strip()
		header = re.match('^>',line)
		if header:
			fstring = fstring + '%' + line + '%'
		else:
			fstring = fstring + line
	fstring = fstring + '%'

	SeqIDs = []
	headers = re.findall('>(.*?)%', fstring)
	for SeqID in headers:
		SeqIDs.append(SeqID)

	Seqs = []
	nheaders = re.findall('%([^>]+)%', fstring)
	for Seq in nheaders:
		Seqs.append(Seq.upper())

	return (SeqIDs,Seqs)
ce_SeqIDs, ce_Seqs = fastaparse(g)


sed "s/draw($prev_number;n_)/draw($number;n_)/g" file.txt > tmp

def coordparse(filename):
	fstring = ''
	for line in filename:
		line = line.strip()
		#coords = re.search('([I|V|X]+)\s+(\d+)\s+(\d+)\s+[A-Z,]+\s([+-])', line)
		coords = re.search('([I|V|X]+)\s+(\d+)\s+(\d+)\s+([+-])', line)
		if coords:
		    chrom = coords.group(1)
		    c1 = int(coords.group(2))
		    c2 = int(coords.group(3))
		    strand = coords.group(4)
		else: 
			pass
		index = ce_SeqIDs.index(chrom) #index of relevant chr
		seq = ce_Seqs[index]
		forward = re.match('\+',strand) #is it on forward strand?
		if forward:
			RNA_seq = seq[c1:c2]
			#RNA_seq_us = seq[c1-51:c2] #60bp upstream + RNA_seq
			#RNA_seq_us = seq[c1-61:c1-1] #60bp upstream
			out_line = ">Chr: {chrom} {c1} {c2} +\n{RNA_seq}\n".format(**locals())
			print out_line
		else:
			RNA_seq = seq[c1:c2]
			#RNA_seq_us = seq[c1-1:c2+50] #60bp upstream + RNA_seq
			#RNA_seq_us = seq[c2:c2+60] #60bp upstream
			RNA_seqrc = ''.join(revcom(RNA_seq))
			#RNA_seq_usrc = ''.join(revcom(piRNA_seq_us))
			out_line = ">Chr: {chrom} {c1} {c2} -\n{RNA_seq}\n".format(**locals())
			print out_line
		
run = coordparse(p)

#def coordparse(filename, outfile):
#with open(outfile, 'w+') as f:
#f.write