#!/usr/bin/env python3

__author__ = "Sergey Knyazev"
__email__ = "sergey.n.knyazev@gmail.com"
__date__ = "10/2/18"

from Bio import SeqIO
import gzip
import random
import argparse


class LessReadsThanNeeded(Exception):
    pass


class Subsampler(object):
    def __init__(self, config):
        self.__fastq1 = config["fastq1"]
        self.__fastq2 = config["fastq2"]
        self.__fastq1_out = config["fastq1_out"]
        self.__fastq2_out = config["fastq2_out"]
        self.__n_samples = config["n_samples"]

    def run(self):
        n_reads = 0
        with gzip.open(self.__fastq1, 'rt') as f:
            for _ in SeqIO.parse(f, "fastq"):
                n_reads += 1
        read_indices = set(random.sample(set(range(n_reads)), self.__n_samples))
        with gzip.open(self.__fastq1, 'rt') as if1, \
                gzip.open(self.__fastq2, 'rt') as if2, \
                gzip.open(self.__fastq1_out, 'wt') as of1, \
                gzip.open(self.__fastq2_out, 'wt') as of2:
            SeqIO.write((r[1] for r in enumerate(SeqIO.parse(if1, "fastq")) if r[0] in read_indices), of1, "fastq")
            SeqIO.write((r[1] for r in enumerate(SeqIO.parse(if2, "fastq")) if r[0] in read_indices), of2, "fastq")


def parse_args():
    arg_parser = argparse.ArgumentParser("Reduces the number of reads in the fastq files.")
    arg_parser.add_argument("--n-samples", type=int, help="number of reads in output files")
    arg_parser.add_argument("--fastq1", help="input forward fastq")
    arg_parser.add_argument("--fastq2", help="input backward fastq")
    arg_parser.add_argument("--fastq1_out", help="output forward fastq")
    arg_parser.add_argument("--fastq2_out", help="output backward fastq")
    return arg_parser.parse_args()


def main():
    config = vars(parse_args())
    Subsampler(config).run()


if __name__ == "__main__":
    main()
