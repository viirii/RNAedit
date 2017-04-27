##################################################
# Main program provided to users
# Takes in two arguments: filename, file format
#
##################################################

from sklearn.externals import joblib
import sys
import os
from Bio import SeqIO
import helper
import numpy as np
import getopt
import shutil

# Pre-process the arguments from terminal
def parse_args(args):
    probability_table = False
    if not len(args) >= 3:
        print('Insufficient parameters.')
        exit(1)
    if len(args) == 4:
        if args[1] == '-p':
            probability_table = True
            return probability_table
        else:
            print('Unidentified parameters.')
            print('To calculate the probability for each base, use "-p inputFile fileFormat" as parameters')
    if len(args) == 5:
        if isinstance(args[2], int):
            print('Wrong parameter type')
            exit(1)
    return probability_table


if __name__ == '__main__':

    probability_trigger = False
    threshold_trigger = False
    site_trigger = False

    threshold = 0.5
    filename = ''
    file_format = ''

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hpst:', [])
    except getopt.GetoptError:
        # TODO add a more elegant error information for input, can directly provide help info
        print('error information')
        sys.exit(2)

    if len(args) != 2:
        # TODO same as above
        print('error information')
    filename, file_format = args[0], args[1]

    for opt, arg in opts:
        if opt == '-h':
            # TODO print help information for users
            # -h help -p output_probability -t set_threshold -s output_position filename file_format
            # remember to inform users that the probability estimation might not be accurate
            print('help information')
            sys.exit()
        elif opt == '-s':
            site_trigger = True
        elif opt == '-p':
            probability_trigger = True
        elif opt == '-t':
            threshold_trigger = True
            threshold = float(arg)

    classifier = joblib.load('trained_model')
    if not os.path.exists('./result'):
        os.makedirs('./result')
    else:
        shutil.rmtree('./result')
        os.makedirs('./result')

    for record in SeqIO.parse(filename, file_format):
        sequences = [str(record.seq).replace('T', 'U')]
        labels = helper.gen_label(sequences)
        data = helper.gen_data(sequences, 19, 0, aligned='no')
        predicted_label = classifier.predict(data)
        sequence = sequences[0]
        # output annotated file
        index = np.arange(0, len(record.seq))
        index = index[predicted_label == 1]
        for i in index:
            sequence = sequence[:i] + 'Y' + sequence[i + 1:]
        with open('./result/predicted.' + file_format, 'a') as f:
            f.write('>' + record.description + '\n')
            f.write(sequence + '\n')
        # output probability for each position
        if probability_trigger:
            probability = classifier.predict_proba(data)
            print(record.description)
            with open('./result/predicted_probability.txt', 'a') as f:
                f.write('>' + record.description + '\n')
                f.write('\tnon-pseudoU\tpseudoU\n')
                for i in range(probability.shape[0]):
                    f.write(str(i + 1) + '\t' + str(probability[i, 0]) + '\t' + str(probability[i, 1]) + '\n')
        if site_trigger and not threshold_trigger:
            with open('./result/predicted_sites', 'a') as f:
                f.write('>' + record.description + '\n')
                f.write('\t'.join([str(i) for i in index]) + '\n')
        if threshold_trigger:
            index = np.arange(0, len(record.seq))
            index = index[probability[:, 1] > threshold]
            with open('./result/predicted_sites', 'a') as f:
                f.write('>' + record.description + '\n')
                f.write('\t'.join([str(i) for i in index]) + '\n')
