#!/usr/bin/env python
"""silva_to_qiime2: a script to parse silva release data for silva into Qiime2 usable formats

"""
# Adam Rivers 2019, adapted from from code by by Mike Robeson Feb 2013


import subprocess
import tempfile
import argparse
import os
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def myparser():
    parser = argparse.ArgumentParser(description='silva_to_Qiime2: prepare Silva databases for Qiime')
    parser.add_argument('--infile', '-i', type=str, required=True,
                        help='A silva fasta file.')
    parser.add_argument('--taxafile', '-t', type=str, help="The intermediate taxonomy file ", required=True)
    parser.add_argument('--outfasta', '-o', type=str, help="The trimmed Fastq file, if it \
                        ends in 'gz' it will be gzipped", required=True)
    parser.add_argument('--clusterid', '-c', type=float, help="the percent identity to cluster at", default=0.99)
    parser.add_argument('--threads', '-n', type=int, help="The number of threads to use for clustering", default=2)
    return parser

# Taxonomy parsing functions

def _get_list_and_length(taxonomy_string):
    """Convert taxonomy string to list and length val
    """
    taxalist = taxonomy_string.split(';')
    length = len(taxalist)
    return taxalist, length

def _organize_taxonomy_info(taxonomy_dict):
    """Takes {seqID:taxonomy} and returns {seqID:(['D_1__'],[taxonomy_level])}

    """
    namesdict = {}
    i = 0
    for seq_id, taxonomy in taxonomy_dict.items():
        fixed_string = taxonomy.replace(';;', ';')
        taxalist, length = _get_list_and_length(fixed_string)
        namesdict[seq_id] = (taxalist, length)
        if length > i:
            i = length
    # D means taxonomic depth 0 being root / top
    depth_list = ['D_' + str(n) + '__' for n in range(i)]
    taxdict = {}
    for seq_id, taxonomy_tup in namesdict.items():
        new_tax_list = taxonomy_tup[0] + ['']*(i-taxonomy_tup[1])
        depth_tax_tuple = zip(depth_list, new_tax_list)
        taxdict[seq_id] = depth_tax_tuple
    return taxdict

def _make_tax_string(taxonomy_tup_list):
    ttl = ''.join(taxonomy_tup_list)
    return ttl

def _make_taxonomy_string(taxonomy_dict):
    flatlist = [seq_id+ '\t' + ';'.join(map(_make_tax_string, taxonomy_tup_list))\
          for seq_id, taxonomy_tup_list in taxonomy_dict.items()]
    fl_str = '\n'.join(flatlist)
    return fl_str


# Main silva class
class Silvadata:

    def __init__(self, infile):
        self.infile = infile
        self.taxa_dict = {}
        self.tempdir = tempfile.mkdtemp()
        self.bioseqs = []
        self.tempseqs = None
        self.rep_file = None
        self.uc_file = None

    def parse_fasta(self):
        with open(self.infile, 'r') as ifile:
            records = SeqIO.parse(ifile, 'fasta')
            for rec in records:
                seqid = rec.id
                taxdata = rec.description.split(" ")[1]
                self.taxa_dict[seqid] = taxdata
                dnaseq = rec.seq.back_transcribe()
                self.bioseqs.append(SeqRecord(dnaseq, id=seqid))
        with open(os.path.join(self.tempdir, "temp.fasta"), 'w') as ofile:
            SeqIO.write(self.bioseqs, ofile, "fasta")
            self.tempseqs = (os.path.join(self.tempdir, "temp.fasta"))

    def cluster(self, threads=1, clusterid=0.99):
        try:
            self.rep_file = os.path.join(self.tempdir, "rep.fasta")
            self.uc_file = os.path.join(self.tempdir, "rep.uc")
            parameters = ["vsearch",
                          "--cluster_size", self.tempseqs,
                          "--centroids", self.rep_file,
                          "--uc", self.uc_file,
                          "--strand", "both",
                          "--id", str(clusterid),
                          "--threads", str(threads)]
            p2 = subprocess.run(parameters, stderr=subprocess.PIPE)
            print(p2.stderr.decode('utf-8'))
            p2.check_returncode()
        except subprocess.CalledProcessError as e1:
            print("Could not perform clustering with Vsearch. Error from Vsearch was:\n {}".format(p2.stderr.decode('utf-8')))
            raise e1
        except FileNotFoundError as f1:
            print("Vsearch was not found, make sure Vsearch is installed and executable")
            raise f1


    def create_taxonomy(self):
        organized_taxa = _organize_taxonomy_info(self.taxa_dict)
        taxastring = _make_taxonomy_string(organized_taxa)
        with open(os.path.join(self.tempdir, "taxonomy.txt"), "w") as ofile2:
            ofile2.write(taxastring)

    def return_data(self, taxafile, outfasta):
        shutil.copy(os.path.join(self.tempdir, "taxonomy.txt"), taxafile)
        shutil.copy(self.rep_file, outfasta)

    def remove_tempdir(self):
        shutil.rmtree(self.tempdir)

def main(args=None):
    parser = myparser()
    if not args:
        args = parser.parse_args()
    data = Silvadata(args.infile)
    data.parse_fasta()
    data.cluster(threads=args.threads, clusterid=args.clusterid)
    data.create_taxonomy()
    data.return_data(args.taxafile, args.outfasta)
    data.remove_tempdir()


if __name__ == '__main__':
    main()
