#!/usr/bin/env ipython3
"""
Author: Sergey Knyazev
Email: sergey.n.knyazev@gmail.com
Created: 10/26/2017
"""
import os, errno
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import glob
import shutil
import subprocess
import argparse
import random


class ReadGenerator(object):
    def __init__(self, config):
        self.__config = config
        if "read_length" not in self.__config:
            self.__config["read_length"] = 250
        if "coverage" not in self.__config:
            self.__config["coverage"] = 10000
        if "insert_size" not in self.__config:
            self.__config["insert_size"] = 100
        if "insert_stdev" not in self.__config:
            self.__config["insert_stdev"] = 300
        if "mate_frag" not in self.__config:
            self.__config["mate_frag"] = 500
        if "mate_frag_stdev" not in self.__config:
            self.__config["mate_frag_stdev"] = 50
        if "mate_pulldown_error_p" not in self.__config:
            self.__config["mate_pulldown_error_p"] = 0.3
        if "duplicate_probability" not in self.__config:
            self.__config["duplicate_probability"] = 0.01
        if "out_dir" not in self.__config:
            self.__config["out_dir"] = os.path.curdir
        if "out_prefix" not in self.__config:
            self.__config["out_prefix"] = "sim_reads"
        self.__config["simseq_home"] = os.path.abspath(self.__config["simseq_home"])
        if "simseq_profile" not in self.__config:
            self.__config["simseq_profile"] = os.path.join(self.__config["simseq_home"], 'profiles/miseq_250bp.txt')
        self.__config["out_dir"] = os.path.abspath(self.__config["out_dir"])
        if "work_dir" not in self.__config:
            self.__config["work_dir"] = os.path.join(self.__config["out_dir"], "_1_tmp_1_")
        else:
            self.__config["work_dir"] = os.path.abspath(self.__config["work_dir"])
        self.__config["sample_fasta"] = os.path.abspath(self.__config["sample_fasta"])
        if "primers_fasta" in self.__config:
            self.__config["primers_fasta"] = os.path.abspath(self.__config["primers_fasta"])
        self.__config["picard_home"] = os.path.abspath(self.__config["picard_home"])

        self.__simseq_jar = os.path.join(self.__config["simseq_home"], 'SimSeqNBProject/store/SimSeq.jar')

    def generate_reads(self):
        os.makedirs(self.__config["work_dir"])
        if "primers_fasta" not in self.__config:
            self.__config["primers_fasta"] = os.path.join(self.__config["work_dir"], "primers.fasta")
            self.__generate_primers()
        cur_dir = os.getcwd()

        os.chdir(self.__config["work_dir"])

        self.__separate_sample()
        self.__generate_sam()
        self.__sam2fasta()
        self.__join_fastas()

        os.chdir(cur_dir)
        shutil.rmtree(self.__config["work_dir"])

    def __separate_sample(self):
        primers = list(SeqIO.parse(self.__config["primers_fasta"], 'fasta'))
        seqs = list(SeqIO.parse(self.__config["sample_fasta"], 'fasta'))
        for i in range(len(seqs)):
            seqs[i].seq = Seq(str(primers[0].seq) +
                              str(seqs[i].seq).replace('-', '') +
                              str(primers[1].seq), generic_dna)
            SeqIO.write(seqs[i], os.path.join(self.__config["work_dir"], '{}.fasta'.format(str(i))), 'fasta')

    def __generate_sam(self):
        fastas = glob.glob(os.path.join(self.__config["work_dir"], '*.fasta'))
        samples = list()
        for f in fastas:
            samples.append(os.path.join(self.__config["work_dir"], f))
        freqs = list()
        sample_len = None
        for sample in samples:
            s = list(SeqIO.parse(sample, 'fasta'))
            if not sample_len:
                sample_len = len(str(s[0].seq))
            freqs.append(float(s[0].id.split('_')[-1].split('%')[0]))
        total_freq = sum(freqs)
        read_number = [int(freqs[i]/total_freq * self.__config["coverage"] * sample_len /
                           (2 * self.__config["read_length"]))
                       for i in range(len(freqs))]

        for i in range(len(fastas)):
            p = os.path.splitext(os.path.basename(fastas[i]))[0]
            prefix = p + '_'
            out = os.path.join(self.__config["work_dir"], p + '.sam')

            command = ["java", "-jar",  "-Xmx2048m", self.__simseq_jar,
                       "-1", str(self.__config["read_length"]), "-2", str(self.__config["read_length"]),
                       "--error", self.__config["simseq_profile"], "--error2", str(self.__config["simseq_profile"]),
                       "--insert_size", str(self.__config["insert_size"]),
                       "--insert_stdev", str(self.__config["insert_stdev"]),
                       "--mate_pair", "--mate_frag", str(self.__config["mate_frag"]),
                       "--mate_frag_stdev", str(self.__config["mate_frag_stdev"]),
                       "--mate_pulldown_error_p", str(self.__config["mate_pulldown_error_p"]),
                       "--read_number", str(read_number[i]), "--read_prefix", prefix,
                       "--reference", samples[i],
                       "--duplicate_probability", str(self.__config["duplicate_probability"]),
                       "--out", out]
            print(" ".join(command))
            subprocess.call(command)

    def __sam2fasta(self):
        sams = glob.glob(os.path.join(self.__config["work_dir"], '*.sam'))
        for f in sams:
            name = os.path.splitext(os.path.basename(f))[0]
            s = os.path.join(self.__config["work_dir"], f)
            b = os.path.join(self.__config["work_dir"], name + '.bam')
            r = os.path.join(self.__config["work_dir"], name + '.fasta')
            bp = name + '_sorted.bam'
            bs = os.path.join(self.__config["work_dir"], bp)
            fasta1 = os.path.join(self.__config["work_dir"], name + '.1.fastq')
            fasta2 = os.path.join(self.__config["work_dir"], name + '.2.fastq')

            command = ["samtools", "view", "-Sb", s, "-T", r]
            with open(b, 'w') as f:
                subprocess.call(command, stdout=f)

            command = ["samtools", "sort", b]
            with open(bp, 'w') as f:
                subprocess.call(command, stdout=f)

            command = ["samtools", "index", bs]
            subprocess.call(command)

            command = ["java", "-jar", "-Xmx2048M", os.path.join(self.__config["picard_home"], "SamToFastq.jar"),
                       "INPUT={}".format(bs), "FASTQ={}".format(fasta1), "SECOND_END_FASTQ={}".format(fasta2),
                       "INCLUDE_NON_PF_READS=true", "VALIDATION_STRINGENCY=SILENT"]
            subprocess.call(command)

            command = ["sed", '-i', "s/\\(^@[0-9][0-9]*_[a-zA-Z0-9]*\\)\\/1/\\1/g", fasta1]
            subprocess.call(command)

            command = ["sed", '-i', "s/\\(^@[0-9][0-9]*_[a-zA-Z0-9]*\\)\\/2/\\1/g", fasta2]
            subprocess.call(command)

    def __join_fastas(self):
        fastas1 = sorted(glob.glob(os.path.join(self.__config["work_dir"], '*.1.fastq')))
        fastas2 = sorted(glob.glob(os.path.join(self.__config["work_dir"], '*.2.fastq')))
        final_fasta1 = os.path.join(self.__config["out_dir"], self.__config["out_prefix"] + '.1.fastq')
        final_fasta2 = os.path.join(self.__config["out_dir"], self.__config["out_prefix"] + '.2.fastq')

        os.makedirs(os.path.dirname(final_fasta1), exist_ok=True)
        os.makedirs(os.path.dirname(final_fasta2), exist_ok=True)

        with open(final_fasta1, 'w'):
            pass
        with open(final_fasta2, 'w'):
            pass

        with open(final_fasta1, 'a') as fo:
            for fi in fastas1:
                subprocess.call(["cat", fi], stdout=fo)
        with open(final_fasta2, 'a') as fo:
            for fi in fastas2:
                subprocess.call(["cat", fi], stdout=fo)

    def __generate_primers(self):
        alph = ["A", "C", "G", "T"]
        seq = "".join(random.choice(alph) for _ in range(self.__config["read_length"]))
        seqs = [SeqRecord(Seq(seq),id=str(i)) for i in range(2)]
        SeqIO.write(seqs, self.__config["primers_fasta"], "fasta")


def parse_args():
    arg_parser = argparse.ArgumentParser("Generator of simulated MiSeq reads for viral populations.")
    arg_parser.add_argument("-c", "--coverage", type=int, help="mean coverage for the sequencing experiment")
    arg_parser.add_argument("-o", "--out-dir", help="output path")
    arg_parser.add_argument("-s", "--simseq-home", help="simseq tool path")
    arg_parser.add_argument("-p", "--picard-home", help="picard tool path")
    arg_parser.add_argument("-i", "--viral-references", help="fasta file with references for viral population")
    return arg_parser.parse_args()


def main():
    config = vars(parse_args())
    config["sample_fasta"] = config.pop("viral_references")
    ReadGenerator(config).generate_reads()


if __name__ == "__main__":
    main()
