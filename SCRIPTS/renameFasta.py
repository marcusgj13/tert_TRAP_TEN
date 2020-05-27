
fileIn = 'TRAPDomain_only.txt'
fileOut = 'TRAPDomain_only_renamed.fasta'

orgNames = []
seq = []
with open(fileIn,'r') as fin:
    for line in fin:
        if line[0] == '>':
            orgNames.append(line.split('>')[1].replace(' ','_' ))
        else:
            seq.append(line.strip())

with open(fileOut,'w') as fout:
    for i in range(0, len(seq)):
        fout.write('> ' + orgNames[i][1:])
        fout.write(seq[i]+'\n')
        
