# Predicting Pseudouridylation of Individual RNA Bases Using Machine Learning#
2017 Spring Bioinformatics Practicum @ CMU

Team Members : Christine Baek, Kevin Chon, Deepank Korandla, Tianqi Tang



## How to Use

### Dependencies and Installation

- Python 2 or 3
- Scikit Learn

### Input & Output 

- Input : Standard bioinformatics file format such as`.fasta` or `.mfa` that contains your input sequence (DNA or RNA). For full list of acceptable input file formats, please visit http://biopython.org/wiki/SeqIO
- Output : Text file containing the updated sequence, with predictions for each position. Each base predicted to be pseudouridine will be marked with the character `Y`. This file will be located in the `results/` directory

### Sample Commands

Below are some sample commands that can be used

- `python test.py {inputfilepath} {fileformat}`
- `python test.py seq.fasta fasta`
- `python test.py seq.mfa fasta``



