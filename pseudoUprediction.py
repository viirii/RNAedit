##################################################
# Main program provided to users
# Takes in two arguments: filename, file format
#
##################################################

from sklearn.externals import joblib
import sys
import os
from Bio import SeqIO
import classification
import numpy as np


# Pre-process the arguments from terminal
def parse_args(args):
    if not len(args) == 3:
        print('Insufficient parameters.')
        exit(1)
    return args[1], args[2]

if __name__ == '__main__':
    filename, file_format = parse_args(sys.argv)

    classifier = joblib.load('trained_model')
    if not os.path.exists('./result'):
        os.makedirs('./result')

    #sequences = []
    for record in SeqIO.parse(filename, file_format):
        sequences = [str(record.seq)]
        labels = classification.gen_label(sequences)
        data = classification.gen_data(sequences, 19, 0, aligned='no')
        predicted_label = classifier.predict(data)
        sequence = sequences[0]
        #np.savetxt('./result/'+record.id, predicted_label, fmt='%d')
        index = np.arange(0, len(record.seq))
        index = index[predicted_label == 1]
        for i in index:
            sequence = sequence[:i] + 'Y' + sequence[i+1:]
        with open('./result/'+record.id, 'w') as f:
            f.write(sequence)


    #for record in SeqIO.parse(filename, file_format):

