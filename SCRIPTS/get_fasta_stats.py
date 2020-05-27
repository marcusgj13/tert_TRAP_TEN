from sys import argv
import numpy as np


def fasta_stats(filename):
    seq_lengths = []
    with open(filename,'r') as fin:
        for line in fin:
            if line[0] != '>':
                seq_lengths.append(len(line))

    print("The longest sequence is %f \n" % np.max(seq_lengths))
    print("The shortest sequence is %f \n" % np.min(seq_lengths))
    print("The average sequence length is %f \n" % np.mean(seq_lengths))


file_in = argv[1]
fasta_stats(file_in)
