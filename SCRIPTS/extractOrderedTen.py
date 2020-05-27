from sys import argv
import matplotlib.pyplot as pl


def find_ordered_region(seq, disorder_values):
    """This function is used to determine the beginning and end of an ordered region
    within the TEN domain of TERT. It uses disopred to first predict disordered regions
    and then looks for *'s at the beginning and end to find the end of N-terminal and
    C-terminal disordered regions."""
    ind = 0
    consecutive_disorder = 0
    for x in disorder_values:
        if x == '*':
            consecutive_disorder += 1
        else:
            consecutive_disorder = 0
        ind += 1
        if consecutive_disorder == 10 and ind > 30:
            break

    partial_seq = ''.join(seq[:ind-9])
    return partial_seq


def parse_diso_file(file_name):
    """This function reads in a .diso file that is the output from disopred and returns
    the disorder profile and the sequence of AAs as lists."""
    full_sequence = []
    disorder_profile = []
    count = 0
    with open(file_name, 'r') as fin:
        for line in fin:
            if count > 2:
                full_sequence.append(line.split()[1])
                disorder_profile.append(line.split()[2])
            count += 1
    return full_sequence, disorder_profile


def write_ten_domain(string_seq, org_name, position):
    """This function writes the ordered region of the TEN domain to file based on the
    result of find_ordered_region"""
    if position == 1:
        with open('trimmedTENDomain.txt', 'w') as out_file:
            out_file.write('> ' + org_name + '\n')
            out_file.write(string_seq + '\n')
    else:
        with open('trimmedTENDomain.txt', 'a') as out_file:
            out_file.write('> ' + org_name + '\n')
            out_file.write(string_seq + '\n')


script, fn = argv

ten_len = []
with open(fn, 'r') as file_list:
    pos = 1
    for f in file_list:
        filename = './RESULTS/' + f.strip()
        organism_name = f.split('.')[0].replace('_', ' ')
        print(organism_name)
        sequence, diso_prof = parse_diso_file(filename)
        ten_seq = find_ordered_region(sequence, diso_prof)
        ten_len.append(len(ten_seq))
        write_ten_domain(ten_seq, organism_name, pos)
        pos += 1

fig, ax = pl.subplots(nrows=1, ncols=1)
ax.grid(axis='y')
ax.hist(ten_len, bins=20, rwidth=0.9, color='#607c8e', alpha=0.85)
ax.set_title('Histogram of ordered TEN domain lengths')
ax.set_xlabel('Number of AAs')
ax.set_ylabel('Frequency')

pl.savefig('Ordered_ten_length_distribution.png')
