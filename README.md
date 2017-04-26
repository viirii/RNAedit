# Predicting Pseudouridylation of Individual RNA Bases Using Machine Learning
2017 Spring Bioinformatics Practicum @ CMU

Team Members : Christine Baek, Kevin Chon, Deepank Korandla, Tianqi Tang



## Introduction

This program predicts the location of pseduouridines in the given nucleotide
sequences using the support vector machine model. This model was trained on
human rRNA and mRNA sequences. The program takes in either single or multiple
nucleotide sequences and then outputs the probability each position is a
pseduouridine or not, the positions of the predicted pseudouridines, and/or the
sequences annotated with pseduouridines. The user can specify which of the three
output types they want, and the output data is printed to separate files based
on output type (instead of a single file for each sequence).



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
--* Command line argument for `.fasta` or `.mfa`: `fasta`
- Output : There are three possible output types. Only one file is used for each output type; separation by sequence occurs within the file. The files will be located in the `results` directory. The user can specify which of the three output types they would like through the command line (more than one can be selected).
--1. Positional Probability: For each position in each sequence, the program will list the probability that position is a pseudouridine (first probability column) or not a pseudouridine (second probability column).
---* Command line argument: `prob`
--2. Pseudouridine positions: For each sequence, the program will list the positions of predicted pseudouridines (one-indexing).
---* Command line argument: `pos`
--3. Annotated sequences: The program will output a text file containing the updated sequences, with predictions for each position, for all sequences. Each base predicted to be pseudouridine will be marked with the character `Y`.
---* Command line argument: `annotate`

### Sample Commands

Below are some sample commands that can be used

- `python pseudoUprediction.py {inputfilepath} {fileformat} {output types}`
- `python pseudoUprediction.py sample.fasta fasta`
- `python pseudoUprediction.py sample.mfa fasta`
