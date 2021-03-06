#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Théo Jamay et Rebecca Goulancourt"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Théo Jamay et Rebecca Goulancourt"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Théo Jamay et Rebecca Goulancourt"
__email__ = "jamay.theo@gmail.com et rgoulancourt@yahoo.com"
__status__ = "Developpement"


def isfile(path):
	"""Check if path is an existing file.
	  :Parameters:
			  path: Path to the file
	"""
	if not os.path.isfile(path):
		if os.path.isdir(path):
			msg = "{0} is a directory".format(path)
		else:
			msg = "{0} does not exist.".format(path)
		raise argparse.ArgumentTypeError(msg)
	return path


def get_arguments():
	"""Retrieves the arguments of the program.
	  Returns: An object that contains the arguments
	"""
	# Parsing arguments
	parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
									 .format(sys.argv[0]))
	parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
						help="Amplicon is a compressed fasta file (.fasta.gz)")
	parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default=400,
						help="Minimum sequence length for dereplication (default 400)")
	parser.add_argument('-m', '-mincount', dest='mincount', type=int, default=10,
						help="Minimum count for dereplication  (default 10)")
	parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default=100,
						help="Chunk size for dereplication  (default 100)")
	parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default=8,
						help="kmer size for dereplication  (default 10)")
	parser.add_argument('-o', '-output_file', dest='output_file', type=str,
						default="OTU.fasta", help="Output file")
	return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
	"""Read sequences from fasta file.
	Parameters
	----------
	amplicon_file : str
			Name of file that contains sequences in FASTA format.
	minseqlen : Longueur minimum des séquences
	Returns
	-------
	générateur
			sequence
	"""
	with gzip.open(amplicon_file, "rt") as amp_file:
		seq = ""
		for line in amp_file:
			if line.startswith(">"):
				if len(seq) >= minseqlen:
					yield seq
				seq = ""
			else:
				seq += line.strip()
		yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
	"""Compte occurence sequences
	Parameters
	----------
	amplicon_file : str
			Name of file that contains sequences in FASTA format.
	minseqlen : int
			Longueur minimum des séquences
	mincount : int
			Comptage minimum des séquences
	Returns
	-------
	générateur
			sequence,occurence
	"""
	dict_occu = {}
	for seq in read_fasta(amplicon_file, minseqlen):
		if seq not in dict_occu:
			dict_occu[seq] = 1
		else:
			dict_occu[seq] += 1
	for key, occu in sorted(dict_occu.items(), key=lambda t: t[1], reverse=True):
		if occu >= mincount:
			yield [key, occu]


def get_unique(ids):
	""""""
	return {}.fromkeys(ids).keys()


def common(lst1, lst2):
	""""""
	return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
	""""""
	len_seq = len(sequence)
	if len_seq < chunk_size * 4:
		raise ValueError("Sequence length ({}) is too short to be splitted in 4"
						 " chunk of size {}".format(len_seq, chunk_size))
	return [sequence[i:i+chunk_size]
			for i in range(0, len_seq, chunk_size)
			if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
	"""Cut sequence into kmers"""
	for i in range(0, len(sequence) - kmer_size + 1):
		yield sequence[i:i+kmer_size]


def get_identity(alignment_list):
	"""Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
	Retourne le pourcentage d'identite entre les deux."""
	id_nu = 0
	for i in range(len(alignment_list[0])):
		if alignment_list[0][i] == alignment_list[1][i]:
			id_nu += 1
	return round(100.0 * id_nu / len(alignment_list[0]), 2)


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
	"""kmer present dans sequence non chimèrique
	Parameters
	----------
	kmer_dict : dict
			dictionnaire kmer
	sequence : str
			sequence non chimerique
	id_seq : int
			identifiant sequence non chimerique
	kmer_size : int
			Longueur des “kmer”
	Returns
	-------
	dict
			dictionnaire kmer
	"""
	for kmer in cut_kmer(sequence, kmer_size):
		if kmer not in kmer_dict:
			kmer_dict[kmer] = [id_seq]
		elif kmer in kmer_dict and id_seq not in kmer_dict[kmer]:
			kmer_dict[kmer].append(id_seq)
	return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
	"""recherche 2 sequences parents
	Parameters
	----------
	kmer_dict : dict
			dictionnaire kmer
	sequence : str
			sequence candidate
	kmer_size : int
			Longueur des “kmer”
	Returns
	-------
	tuple
			index des 2 sequences parentes
	"""
	list_id = []
	for kmer in cut_kmer(sequence, kmer_size):
		if kmer in kmer_dict:
			for value in kmer_dict[kmer]:
				list_id.append(value)
	list_best = Counter(list_id).most_common(2)
	if list_best[0][1] >= 2:
		return list_best[0][0], list_best[1][0]
	return -1


def detect_chimera(perc_identity_matrix):
	"""Determine si une sequence est chimérique ou pas
	Parameters
	----------
	perc_identity_matrix : matrice
			matrice des pourcentages id entre candidat et les 2 parents
	Returns
	-------
	Booléen
			True si c'est une sequence chimérique
	"""
	list_std = []
	par1 = False
	par2 = False
	for data in perc_identity_matrix:
		list_std.append(round(statistics.stdev(data), 2))
		if data[0] > data[1]:
			par1 = True
		elif data[0] < data[1]:
			par2 = True
	if statistics.mean(list_std) > 5:
		if par1 and par2:
			return True
	return False


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""fonction global pour déterminer si une sequence est chimérique ou pas
	Parameters
	----------
	amplicon_file : str
			Name of file that contains sequences in FASTA format.
	minseqlen : int
			Longueur minimum des séquences
	mincount : int
			Comptage minimum des séquences
	chunk_size : int
			Taille des partitions de séquence
	kmer_size : int
			Longueur des “kmer”
	Returns
	-------
	generateur
			sequence,occurence
	"""
	kmer_dict = {}
	seqs = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
	yield seqs[0]
	yield seqs[1]
	kmer_dict = get_unique_kmer(kmer_dict, seqs[0][0], 0, kmer_size)
	kmer_dict = get_unique_kmer(kmer_dict, seqs[1][0], 1, kmer_size)
	for i in range(2, len(seqs)):
		chunks_cand = get_chunks(seqs[i][0], chunk_size)
		perc_identity_matrix = []
		for j in range(len(chunks_cand)):
			best = search_mates(kmer_dict, chunks_cand[j], kmer_size)
			if best != -1 :
				chunk_par1 = get_chunks(seqs[best[0]][0], chunk_size)
				chunk_par2 = get_chunks(seqs[best[1]][0], chunk_size)
				list_tmp = []
				list_tmp.append(get_identity(nw.global_align(chunks_cand[j], chunk_par1[j],
					gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))))
				list_tmp.append(get_identity(nw.global_align(chunks_cand[j], chunk_par2[j],
					gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))))
				perc_identity_matrix.append(list_tmp)
		if detect_chimera(perc_identity_matrix) is False:
			kmer_dict = get_unique_kmer(kmer_dict, seqs[i], i, kmer_size)
			yield seqs[i]


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
	"""regroupement glouton des sequences otu
	Parameters
	----------
	amplicon_file : str
			Name of file that contains sequences in FASTA format.
	minseqlen : int
			Longueur minimum des séquences
	mincount : int
			Comptage minimum des séquences
	chunk_size : int
			Taille des partitions de séquence
	kmer_size : int
			Longueur des “kmer”
	Returns
	-------
	list
			[sequence otu,occurence]
	"""
	seqs = list(chimera_removal(amplicon_file, minseqlen,
				mincount, chunk_size, kmer_size))
	list_otu = [seqs[0]]
	for seq in seqs[1:]:
		flag = True
		for seq_otu in list_otu:
			if get_identity(nw.global_align(seq[0], seq_otu[0], gap_open=-1, gap_extend=-1,
					matrix=os.path.abspath(os.path.join(os.path.dirname(__file__), "MATCH")))) >= 97:
				flag = False
				break
		if flag :
			list_otu.append(seq)
	return list_otu


def fill(text, width=80):
	"""Split text with a line return to respect fasta format"""
	return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
	"""fichier output contenant la liste des otu au format fasta
	Parameters
	----------
	OTU_list : list
			[sequence otu,occurence] : liste seq otu avec leur occurence
	output_file : str
			nom du fichier de sortie
	Returns
	-------
	list
			[sequence otu,occurence]
	"""
	with open(output_file, "w") as otu_file:
		for i in range(len(OTU_list)):
			otu_file.write(
				">OTU_{} occurrence:{}\n".format(i+1, OTU_list[i][1]))
			otu_file.write(fill(OTU_list[i][0])+"\n")


# ==============================================================
# Main program
# ==============================================================
def main():
	"""
	Main program function
	"""
	# Get arguments
	args = get_arguments()
	# Votre programme ici
	OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
	write_OTU(OTU_list, args.output_file)


if __name__ == '__main__':
	main()
