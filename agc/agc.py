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

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file,"rt") as amp_file :
        seq = ""
        for line in amp_file :
            if line.startswith(">") :
                if len(seq) >= minseqlen :
                    yield seq
                seq = ""
            else :    
                seq += line.strip()
        yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    dict_occu = {}
    for seq in read_fasta(amplicon_file, minseqlen) :
        if seq not in dict_occu :
            dict_occu[seq] = 1
        else :
            dict_occu[seq] += 1
    for key,occu in sorted(dict_occu.items(), key=lambda t:t[1],reverse=True) :
        if occu >= mincount :
            yield [key,occu]


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
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
    for kmer in cut_kmer(sequence,kmer_size) :
        if kmer not in kmer_dict :
            kmer_dict[kmer] = [id_seq]
        elif kmer in kmer_dict and id_seq not in kmer_dict[kmer] :
            kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    list_id = []
    for kmer in cut_kmer(sequence,kmer_size) :
        if kmer in kmer_dict :
            for value in kmer_dict[kmer] :
                list_id.append(value)
    list_best = Counter(list_id).most_common(2)
    return list_best[0][0],list_best[1][0]


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    kmer_dict = {}
    list_otu = abundance_greedy_clustering(amplicon_file, minseqlen, mincount,chunk_size,kmer_size)
    list_chim = [list_otu[0],list_otu[1]]
    for seq_otu in list_otu :
        for i in range(len(list_chim)) :
            kmer_dict = get_unique_kmer(kmer_dict,list_chim[i][0],i,8)
        chunks = get_chunks(seq_otu[0], 37)
        for chunk in chunks :
            list_test = search_mates(kmer_dict,chunk,kmer_size)
            #if list_test[0][1] >= 2 and list_test[1][1] >= 2 :



def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    seqs=list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    list_otu = [[seqs[0][0],seqs[0][1]]]
    for seq in seqs[1:] :
        for seq_otu in list_otu :
            if get_identity(nw.global_align(seq[0], seq_otu[0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))) < 97 :
                list_otu.append([seq[0],seq[1]])
                break
    return list_otu


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    with open(output_file,"w") as otu_file :
        for i in range(len(OTU_list)) :
            otu_file.write(">OTU_{}occurrence:{}\n".format(i,OTU_list[i][1]))
            otu_file.write(fill(OTU_list[i][0])+"\n")



#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Votre programme ici
    #read_fasta(args.amplicon_file, args.minseqlen)
    # for seq in read_fasta(args.amplicon_file, args.minseqlen) :
    #     print(seq)
    # for seq in dereplication_fulllength(args.amplicon_file, args.minseqlen, args.mincount) :
    #     print(seq)
    # OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    # output_file = "OTU.fasta"
    # write_OTU(OTU_list, output_file)
    chimera_removal(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)

if __name__ == '__main__':
    main()