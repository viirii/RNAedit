#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 14:35:54 2017

@author: TT
"""

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from Bio import AlignIO
import cb21
from sklearn.svm import SVC
from sklearn.externals import joblib
from sklearn.model_selection import cross_val_score


def readfile(filename, format):
    alignment = AlignIO.read(filename, format)
    return alignment


# Transfer a string of DNA sequence into categorical labeled sequence
def nuc_to_label(sequence):
    label = np.zeros(len(sequence))
    for i in range(len(sequence)):
        if sequence[i] == '-':
            label[i] = 0
        elif sequence[i] == 'A':
            label[i] = 1
        elif sequence[i] == 'T':
            label[i] = 2
        elif sequence[i] == 'C':
            label[i] = 3
        elif sequence[i] == 'G':
            label[i] = 4
        elif sequence[i] == 'Y':
            label[i] = 2
    return label


# Calculate the GC content for a given sequence
def gc_content(sequence):
    if isinstance(sequence, str):
        sequence = nuc_to_label(sequence)
    return np.sum(np.logical_or(sequence == 3, sequence == 4)) / sequence.shape[0]


def gen_window(data, window_sz):
    window = np.zeros((data.shape[0] * data.shape[1], window_sz))
    half_window = np.floor(window_sz / 2)
    for row in range(data.shape[0]):
        for col in range(data.shape[1]):
            if half_window <= col < data.shape[1] - window_sz:
                window[row * data.shape[1] + col, :] = data[row, int(col - half_window):int(col + half_window + 1)]
            elif col < half_window:
                window[row * data.shape[1] + col, int(-col + half_window):] = data[row, 0:int(col + half_window + 1)]
            elif col >= data.shape[1] - window_sz:
                window[row * data.shape[1] + col, 0:data.shape[1] - col] = data[row, col:]
    return window


def motif_prob(sequence, motif):
    sequence = nuc_to_label(sequence)
    probability = 0
    if len(sequence) == motif.shape[1]:
        for i in range(len(sequence)):
            probability += np.log(motif[sequence[i], i])
    return probability


def gen_data(raw_data, window_sz, feature_num, motif=''):
    transformed = np.zeros((len(raw_data), len(raw_data[0])))
    for i in range(len(raw_data)):
        transformed[i, :] = nuc_to_label(raw_data[i])
    window_data = gen_window(transformed, window_sz)
    data = np.zeros((window_data.shape[0], window_data.shape[1] + feature_num))
    data[0:window_data.shape[0], 0:window_data.shape[1]] = window_data
    basis = window_data.shape[1]
    if not motif == '':
        for i in range(window_data.shape[0]):
            data[i, basis + 0] = motif_prob(window_data[i, :], motif)
    return data


def convert_label(label):
    converted_label = np.zeros(label.shape[0] * label.shape[1])
    for i in range(label.shape[0]):
        converted_label[i * label.shape[1]:(i + 1) * label.shape[1]] = label[i, :]
    return converted_label


def gen_fold(X, Y, fold):
    X_fold = np.zeros((fold, X.shape[0], X.shape[1]))
    Y_fold = np.zeros((fold,))

    return X_, Y_fold


# def train_with_svm(X, Y):
#     clf = SVC(class_weight='balanced')
#     clf.fit(X, Y)
#     return clf

def score_helper(clf, X, Y):
    pred_Y = clf.predict(X)
    f1 = f1_score(Y, pred_Y, average='weighted')
    cnf_matrix = confusion_matrix(Y, pred_Y)
    print(cnf_matrix)
    return f1


def test_with_svm(X, Y, fold):
    clf = SVC(class_weight='balanced')
    scores = cross_val_score(clf, X, y=Y, scoring='f1_weighted', cv=fold, n_jobs=-1)
    print(scores)


def store_full_classifier(X, Y):
    clf = SVC(class_weight='balanced')
    clf.fit(X, Y)
    joblib.dump(clf, 'trained_model')


def main():
    # Read in *.aln file

    # filename = 'STEP4_TrimmedNT.aln'
    # fmt = 'clustal'

    # alignment = readfile(filename, fmt)
    # sequences = []
    # for record in alignment:
    #     sequences.append(record.seq)

    # Get sequences using Christine's script
    filename = '21.txt'
    sequences, labels = cb21.middleU(filename)

    data = gen_data(sequences, 7, 0)
    labels = convert_label(labels)

    for window_sz in range(3, 15, 2):
        data = gen_data(sequences, window_sz, 0)
        test_with_svm(data, labels, 10)
    store_full_classifier(data, labels)

if __name__ == '__main__':
    main()
