from sys import argv
import matplotlib.pyplot as pl

def parse_horiz_file(file_name):
    """Read through the horiz file output by Psipred """
    confidence = []
    prediction = []
    sequence = []

    with open(file_name, 'r') as fin:
        for line in fin:
            if line.strip():
                if line[0] == "C":
                    confidence.append(line.split(':')[1].strip())
                elif line[0] == "P":
                    prediction.append(line.split(':')[1].strip())
                elif line[2] == "A":
                    sequence.append(line.split(':')[1].strip())

        confidence = "".join(confidence)
        prediction = "".join(prediction)
        sequence = "".join(sequence)

    return confidence, prediction, sequence


def locate_helices(ss_prediction, confidence_score):
    """This function will be used to locate the sequence between IFDa and IFDc
    which are both helical structures. These should be well predicted and should
    be a better marker."""
    IFDa_start_found = 0
    IFDa_end = 1
    for i in range(0, len(ss_prediction)):
        if ss_prediction[i] == 'H':
            if int(confidence_score[i]) >= 8:
                if IFDa_start_found == 0:
                    IFDa_start_found = 1
        elif ss_prediction[i] != 'H':
            if IFDa_start_found == 1:
                IFDa_end = i
                break

    IFDc_end_found = 0
    IFDc_start = len(ss_prediction)
    for i in range(len(ss_prediction)-1, 0, -1):
        if ss_prediction[i] == 'H':
            if int(confidence_score[i]) >= 8:
                if IFDc_end_found == 0:
                    print("Found IFDc at residue " + str(i))
                    IFDc_end_found = 1
        elif ss_prediction[i] != 'H':
            if IFDc_end_found == 1:
                IFDc_start = i
                break


    return IFDa_end, IFDc_start


def write_trap(seq, start, end, org_name, position):
    """Write the extracted trap sequence to file"""
    if position == 1:
        with open('TRAPDomain_only.txt', 'w') as out_file:
            out_file.write('> ' + org_name + '\n')
            if seq[start:end].strip():
                out_file.write(seq[start:end] + '\n')
            else:
                out_file.write('No TRAP\n')
    else:
        with open('TRAPDomain_only.txt', 'a') as out_file:
            out_file.write('> ' + org_name + '\n')
            if seq[start:end].strip():
                out_file.write(seq[start:end] + '\n')
            else:
                out_file.write('No TRAP\n')

script, fn = argv

TRAP_len = []
with open(fn, 'r') as file_list:
    pos = 1
    for f in file_list:
        filename = './RESULTS/' + f.strip()
        organism_name = f.split('.')[0].replace('_', ' ')
        print(organism_name)
        psipred_conf, psipred_pred, TRAP_seq = parse_horiz_file(filename)
        seq_start, seq_end = locate_helices(psipred_pred, psipred_conf)
        if (seq_end - seq_start) > 0:
            TRAP_len.append(seq_end - seq_start)
        else:
            TRAP_len.append(0)
        write_trap(TRAP_seq, seq_start, seq_end, organism_name, pos)
        pos += 1

fig, ax = pl.subplots(nrows=1, ncols=1)
ax.grid(axis='y')
ax.hist(TRAP_len, bins=20, rwidth=0.9, color='#607c8e', alpha=0.85)
ax.set_title('Histogram of TRAP domain lengths')
ax.set_xlabel('Number of AAs')
ax.set_ylabel('Frequency')

pl.savefig('TRAP_length_distribution.png')