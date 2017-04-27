# Predicting Pseudouridylation of Individual RNA Bases Using Machine Learning
2017 Spring Bioinformatics Practicum @ CMU

Team Members : Christine Baek, Kevin Chon, Deepank Korandla, Tianqi Tang



## Introduction

This program predicts the location of pseduouridines in the given nucleotide sequences using the support vector machine model (SVM). This model was trained on human rRNA and mRNA sequences and uses boundary-based decision-making. The program takes in either single or multiple nucleotide sequences and then outputs a single file containing all the sequences such that all the positions predicted to be pseduouridine are updated accordingly. The user  has the option of asking the program to output a file listing the probabilities that each position in each sequence will either be pseudouridylated or not. Additionally, the user can ask the program to output a file listing the positions of the predicted pseudouridines in each sequence. For the second option, the user can also specify a threshold such that only those positions with a probability greater than the threshold are predicted to be pseudouridines. It should be noted that probability-based prediction is not reliable due to the nature of the SVM.



## How to Use

### Dependencies and Installation

- Python 2 or 3
- Scikit Learn
- bioPython
- numpPy
- scipPy


Below are command for package installation for MacOS

`pip install -U scikit-learn`

`pip install numpy`

`pip install scipy`

`pip install biopython`



### Input & Output

- Input : Standard bioinformatics file format such as `.fasta` or `.mfa` that contains your input sequence (DNA or RNA). For a full list of acceptable input file formats, please visit http://biopython.org/wiki/SeqIO
  * File format command line argument for `.fasta` or `.mfa` files: `fasta`
- Output : The program will always output a text file containing each updated sequence such that the positions predicted to be pseudouridine are marked with the character `Y` (instead of a T or a U). The user can also use command line arguments to ask for two additional text files to be outputted (one or both can be selected). All files will be located in the `results` directory, and separation by sequence occurs within each file.
  1. Positional Probability: For each position in each sequence, the program will list the probability that position is a pseudouridine (first probability column) or not a pseudouridine (second probability column).
    * Command line argument: `-p`
  2. Pseudouridine positions: For each sequence, the program will list the positions of predicted pseudouridines.
    * Command line argument: `-s`
    * The user can also specify the threshold at which to call a particular position a pseudouridine.
        * Command line argument: `-t <floating point number>`
        * The `-s` argument will be ignored if the `-t` argument is present.

### Sample Commands

Below are some sample commands that can be used

- `python pseudoUprediction.py {output commands} {inputfilepath} {fileformat} `
- `python pseudoUprediction.py sample.fasta fasta`
- `python pseudoUprediction.py -p sample.mfa fasta`
- `python pseudoUprediction.py -s sample.mfa fasta`
- `python pseudoUprediction.py -t 0.9 sample.mfa fasta`
- `python pseudoUprediction.py -p -s -t 0.9 sample.mfa fasta`
