##################################################
# Main program provided to users
# Takes in two arguments: filename, file format
#
##################################################

from sklearn.externals import joblib
import sys

# Pre-process the arguments from terminal
def parse_args(args):
    if not len(args) == 2:
        print('Insufficient parameters.')
        exit(1)
    return args[0], args[1]

if __name__ == '__main__':
    filename, file_format = parse_args(sys.argv)
    classifier = joblib.load('trained_model')
    predicted_label = classifier.predict()
    print(predicted_label)
