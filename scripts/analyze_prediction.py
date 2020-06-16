#!/usr/bin/env python
"""
.. module:: emd_haplo
    ::synopsis: A module which calculate distance between two sets of genomic sequences.
    It uses EMD described at CliequeSNV article (https://www.biorxiv.org/content/early/2018/03/31/264242). 

.. moduleauthor:: Sergey Knyazev <sergey.n.knyazev@gmail.com>
"""

import argparse
from Bio import SeqIO
from emd import emd
"""emd module address: https://github.com/garydoranjr/pyemd.git"""
import numpy as np
import json
import sys


class Haplotypes(object):
    def __init__(self, seqs, freqs):
        self.seqs = seqs
        self.freqs = freqs

    def add_zero_freqs(self, diff):
        self.freqs = np.append(self.freqs, np.zeros(diff))


# Parses the fasta file to get the frequencies. strip is to remove the % and get the frequency out of 100
# Returns a Haplotype object with the attributes: a list of sequences, a list of frequencies out of 1
def load_fasta(fasta_file):
    seqs = list(SeqIO.parse(open(fasta_file), 'fasta'))
    freqs = np.array([float(x.id.split('_')[-1].strip('%')) for x in seqs])
    return Haplotypes([str(x.seq) for x in seqs], freqs/sum(freqs))

# map applies the function lambda to all of the elements in zip(seq1, seq2)
# zip returns iterator of tuples composed of the i-th from each of seq1 and seq2
# the lambda checks if the tuple has unequal elements
# This function returns how many nucleotides are different in the sequences seq1 and seq2
NUCL = ['A', 'C', 'G', 'T']
def h_dist(seq1, seq2):
    return sum(map(lambda x: x[0] in NUCL and x[1] in NUCL and x[0] != x[1], zip(seq1, seq2)))


# Input: 2 lists of sequences
# This returns 2D np array with the number of nucleotide differences for each possible pair in the 2 lists of sequences
# Rows are predictions, Columns are answers
# For example: 1st row of 2D np array = list of h_dist of all key_seqs compared to the 1st prediction_seq
def get_dist(prediction_seqs, key_seqs):
    return np.array([[float(h_dist(prediction_seqs[i], key_seqs[j]))
                      for j in range(len(key_seqs))]
                     for i in range(len(prediction_seqs))])


# Input: an iterable A
# Returns a tuple (smallest element in A, list of indices of the smallest elements in A)
def locate_min(a):
    smallest = min(a)
    return smallest, [index for index, element in enumerate(a)
                      if smallest == element]


# dist is the 2D (m x n) np array from the get_dist function
# locate_min acts on each column (not row!)
# Returns: an array of n tuples in the form (smallest value in column 0<=i<n, indices of smallest value)
def get_prediction_closest_to_answer(dist):
    return [locate_min(dist[:, j]) for j in range(len(dist[0]))]


def get_answer_closest_to_prediction(dist):
    return [locate_min(dist[i, :]) for i in range(len(dist))]


# emd_flows & dist are 2D np arrays with the same shape
# Iterating along the columns again
# sum(emd_flows[:, j] * dist[:, j]) gets dot product of the jth columns of emd_flows & dist
# the sum/sum code essentially uses weights emd_flows with dist and then divides by sum of emd_flows
# this repeats for each column and returns a separate value for each column.
def get_freq_adjusted_mismatches(emd_flows, dist):
    return [sum(emd_flows[:, j] * dist[:, j])/sum(emd_flows[:, j]) for j in range(len(emd_flows[0]))]


def get_closest_freq(test_freqs, closest_to_key_hapls, key_freqs):
    return [min(x)[0]
            for x in [[(test_freqs[j], abs(key_freqs[i] - test_freqs[j]))
                for j in closest_to_key_hapls[i][1]]
                    for i in range(len(key_freqs))]]


def equalize(prediction_hapls, key_hapls, dist):
    diff = len(prediction_hapls.freqs) - len(key_hapls.freqs)
    if diff > 0:
        key_hapls.add_zero_freqs(diff)
        a = np.zeros((len(prediction_hapls.freqs), len(prediction_hapls.freqs)))
        a[:, :-diff] = dist
        return a
    elif diff < 0:
        prediction_hapls.add_zero_freqs(-diff)
        a = np.zeros((len(key_hapls.freqs), len(key_hapls.freqs)), axis=0)
        return a
    else:
        return dist


def get_adc(predictions_closest_to_answer, predictions_closest_to_answer_freqs):
    adc = 0.0
    for i in range(len(predictions_closest_to_answer)):
        adc += predictions_closest_to_answer[i][0]*predictions_closest_to_answer_freqs[i]
    return adc


def get_uadc(predictions_closest_to_answer, seq_len):
    fnr = 0.0
    for i in range(len(predictions_closest_to_answer)):
        fnr += predictions_closest_to_answer[i][0]
    return fnr/len(predictions_closest_to_answer)/seq_len


def report(prediction_hapls, answer_hapls, dist):
    # For the dataset it reports:
    # 1) the count of predicted haplotypes with no errors (TP);
    # 2) the count of predicted haplotypes with at least one error(FP);
    # 3) total count of haplotypes(TP+FP);
    # 4) sensitivity(TP/(TP+FN));
    # 5) Precision(PPV=(TP/(TP+FP));
    # 6) EMD to a consensus.

    # For every true variant it should report:
    # 1) true frequency(TF);
    # 2) editing distance to the closest prediction variant(ECP);
    # 3) frequency of the closest predicted variant(FCP);
    # 4) explanation error for a true variant (EEV).

    #For every predicted variant it should report:
    # 1) editing distance to the closest true variant (ECT).

    emd_res = emd(X=np.ones(len(prediction_hapls.freqs)), Y=np.ones(len(answer_hapls.freqs)),
                  X_weights=prediction_hapls.freqs, Y_weights=answer_hapls.freqs,
                  distance='precomputed', D=dist, return_flows=True)
    pred_freqs_unif = np.array([1./len(prediction_hapls.freqs) for _ in range(len(prediction_hapls.freqs))])
    answer_freqs_unif = np.array([1./len(answer_hapls.freqs) for _ in range(len(answer_hapls.freqs))])

    emd_unif = emd(X=np.ones(len(prediction_hapls.freqs)), Y=np.ones(len(answer_hapls.freqs)),
                   X_weights=pred_freqs_unif, Y_weights=answer_freqs_unif,
                   distance='precomputed', D=dist, return_flows=True)

    ans_hapl_count = len(answer_hapls.seqs)
    pred_hapl_count = len(prediction_hapls.seqs)
    predictions_closest_to_answer = get_prediction_closest_to_answer(dist)
    answer_closest_to_prediction = get_answer_closest_to_prediction(dist)
    predictions_closest_to_answer_freqs = get_closest_freq(prediction_hapls.freqs, predictions_closest_to_answer,
                                                           answer_hapls.freqs)
#    answer_closest_to_prediction_freqs = get_closest_freq(answer_hapls.freqs, answer_closest_to_prediction,
#                                                          prediction_hapls.freqs)
    freq_adjusted_mismatches = get_freq_adjusted_mismatches(emd_res[1], dist)
    report_dict = dict()
    report_dict["TP"] = sum([x[0] == 0 for x in predictions_closest_to_answer[:ans_hapl_count]])
    report_dict["FP"] = len(prediction_hapls.seqs) - report_dict["TP"]
    report_dict["TotalPredicted"] = len(prediction_hapls.seqs)
    report_dict["Sensitivity"] = float(report_dict["TP"])/ans_hapl_count
    report_dict["PPV"] = float(report_dict["TP"])/report_dict["TotalPredicted"]
    report_dict["EMD"] = emd_res[0]
    # Fractional accuracy
    report_dict["UEMD"] = emd_unif[0]
    report_dict["TF"] = [x for x in answer_hapls.freqs[:ans_hapl_count]]
    report_dict["ECP"] = [x[0] for x in predictions_closest_to_answer[:ans_hapl_count]]
    report_dict["ECT"] = [x[0] for x in answer_closest_to_prediction[:pred_hapl_count]]
    report_dict["FCP"] = [x for x in predictions_closest_to_answer_freqs[:ans_hapl_count]]
    report_dict["EEV"] = [x for x in freq_adjusted_mismatches[:ans_hapl_count]]
    report_dict["PCA"] = [x[1][0] for x in predictions_closest_to_answer]
    report_dict["ACP"] = [x[1][0] for x in answer_closest_to_prediction]
    # ADC
    report_dict["ADC"] = get_adc(predictions_closest_to_answer, answer_hapls.freqs)
    # APE
    report_dict["APE"] = get_adc(answer_closest_to_prediction, prediction_hapls.freqs)
    report_dict["UADC"] = get_adc(predictions_closest_to_answer, answer_freqs_unif)
    report_dict["UAPE"] = get_adc(answer_closest_to_prediction, pred_freqs_unif)
    json.dump(report_dict, sys.stdout)


def parse_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("prediction_fasta")
    arg_parser.add_argument("relevant_fasta")
    return arg_parser.parse_args()


def main(args):
    prediction_hapls = load_fasta(args.prediction_fasta)
    relevant_hapls = load_fasta(args.relevant_fasta)
    dist = get_dist(prediction_hapls.seqs, relevant_hapls.seqs)
    report(prediction_hapls, relevant_hapls, dist)


if __name__ == "__main__":
    main(parse_args())
