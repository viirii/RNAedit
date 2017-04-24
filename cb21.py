# http://www.sciencedirect.com/science/article/pii/S0092867414010988 
# this script process above paper's supplemental material, table s3
# into intermediate format, for further processing
# requires 21.txt, which is only the sequences column from mmc3-1.xlsx

import numpy as np

pseudoU = "Y"  # wild-card for 'pseudo U'


def readFile(path):
    with open(path, "rt") as f:
        return f.read()


def middleU(filepath):
    result = []
    labels = []
    read = readFile(filepath)
    for elem in read.split():
        result.append(elem[:10] + pseudoU + elem[11:])

    ulabel = [0 for x in range(9)] + [1] + [0 for x in range(11)]
    for i in range(len(read.split())):
        labels.append(ulabel)

    # Convert to a numpy matrix
    labels = np.matrix(labels)

    return result, labels


result, labels = middleU("21.txt")
