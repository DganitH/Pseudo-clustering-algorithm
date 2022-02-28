# Pseudo-clustering algorithm

The goal of this project is to design and implement a heuristic algorithm for the clustering problem in DNA storage under high error rate from Nanopore sequencing.
My algorithm has two levels:
1- Division according to the index.
2- Division by closeness to existing clusters according to a certain criterion.

For detailed explanation, see the file 'Pseudo-clustering algorithm- Final Project.pdf'.

## Getting Started

After downloading the folder, install the following:
* regex
* statistics
* numpy
* matplotlib

In order to run the program, run:

```
python main.py
```

## The Project Files

The files in the folder are:

_data_files:_
- 450indices.txt - The indices that I concatenated to the DNA sequences.
- errors_shuffled.txt - A file with all the noisy copies of all of the sequences shuffled.
- evyat.txt - A file with a list of all the original sequences and below each one of them, the erroneous copies of it.

_python files:_
- main.py - The file that contains the clustering algorithm. You should run this file.
- functions.py.
- results_functions.py.

_data files for the graphs:_
- perc1.npy
- perc2.npy


## Author

Dganit Hanania- Buchris 
