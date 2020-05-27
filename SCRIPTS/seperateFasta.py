from sys import argv

def writeToFile(name, sequence):

    filename = './INPUT_FILES/' + name.split('>')[1][1:].split(' ')[0][:-1]+ '.fasta'

    with open(filename,'w') as f:
        f.write(name + '\n')
        f.write(sequence)

def readSequences(fileIn):
    count = 1
    for line in open(fileIn,'r'):
        if line[0] == '>':
            protName = line.strip('\n')
        else:
            seq = line.strip('\n')

        if count % 2 == 0:
            writeToFile(protName,seq)

        count += 1

script, fileIn = argv

readSequences(fileIn)
