from collections import defaultdict
import matplotlib.pyplot as pl
import numpy as np
import fnmatch

# Script for locating the position in of B prime
# in a list of tert sequences.

# First load the list of motifs to search for
BprimeMotifs = []

for line in open('BprimeMotifs.txt', 'r'):
    if line[0] != '>':
        BprimeMotifs.append(line.strip())

# Next load in the list of file names so that individual
# sequences can be loaded 1 by 1

FileNames = []

for f in open('FileList.txt', 'r'):
    FileNames.append(f.strip())

print(FileNames)
# Next set up parameters for doing a sliding window search
# through the sequences.

windowSize = 5
organismNames = defaultdict()
notFound = []
# Loop through all files, extract name and sequence and position
# of B prime motif, save this in a dictionary.

for f in FileNames:
    seq = []
    for line in open(f, 'r'):
        if line[0] == '>':
            name = line.split('[')[1][:-2].replace(" ","_")
        else:
            seq.append(line.strip())

    fullSeq = ''.join(seq)

    motifsLoc = []
    foundMotif = 0
    for i in range(0, len(fullSeq) - windowSize):
        subSeq = fullSeq[i:i + windowSize]

        for motif in BprimeMotifs:
            if subSeq == motif:
                motifsLoc.append([i, motif, fullSeq])
                foundMotif += 1
                print(motifsLoc[0][0:1])

    if foundMotif > 0:
        organismNames[name] = motifsLoc
    else:
        notFound.append(name+'\n')

# Write out a file containing the organisms name and the
# location of putative Bprime motifs within its sequence.

with open('BprimeMotifLocation.txt', 'w') as Bp:
    for key in organismNames:
        Bp.write(key + '\t')
        Bp.write(str(organismNames[key][0][0:2]).strip('[').strip(']') + '\n')

# Now we have a dictionary that contains both the sequence
# and the location of the Bprime motif we can loop backwards
# until we find the A motif.

AMotifs = []
IFDSeqs = defaultdict()
seqLens = []
IFDLens = []
windowSize = 6
for line in open('AMotifs.txt', 'r'):
    if line[0] != '>':
        AMotifs.append(line.strip())


for name in organismNames:
    seq = organismNames[name][0][2]
    startPoint = organismNames[name][0][0]

    AMotifLoc = []
    foundMotif = 0

    for ind in range(startPoint, windowSize, -1):
        subSeqA = seq[ind - windowSize:ind]
        if name == 'Danio_rerio' or name == 'Cyprinus_carpio':
            motif = 'YFVKVD'
            if foundMotif == 0:
                if fnmatch.fnmatch(subSeqA, motif) and startPoint - ind > 25:
                     AMotifLoc.append([len(seq), ind, startPoint, subSeqA, seq[ind - 1:startPoint + 1]])
                     IFDLens.append(startPoint-ind)
                     seqLens.append(len(seq)) 
                     foundMotif += 1
        elif name == 'Capsaspora_owczarzaki':
             motif = 'YFVSLD'
             if foundMotif == 0:
                if fnmatch.fnmatch(subSeqA, motif) and startPoint - ind > 25:
                     AMotifLoc.append([len(seq), ind, startPoint, subSeqA, seq[ind - 1:startPoint + 1]])
                     IFDLens.append(startPoint-ind)
                     seqLens.append(len(seq))
                     foundMotif += 1
        elif name == 'Opisthorchis_viverrini':
             motif = 'RLAKMD'
             if foundMotif == 0:
                if fnmatch.fnmatch(subSeqA, motif) and startPoint - ind > 25:
                     AMotifLoc.append([len(seq), ind, startPoint, subSeqA, seq[ind - 1:startPoint + 1]])
                     IFDLens.append(startPoint-ind)
                     seqLens.append(len(seq))
                     foundMotif += 1
        elif name == 'Ramazzottius_varieornatus':
             motif = 'YFVTCD'
             if foundMotif == 0:
                if fnmatch.fnmatch(subSeqA, motif) and startPoint - ind > 25:
                     AMotifLoc.append([len(seq), ind, startPoint, subSeqA, seq[ind - 1:startPoint + 1]])
                     IFDLens.append(startPoint-ind)
                     seqLens.append(len(seq))
                     foundMotif += 1
        elif name == 'Tetrahymena_thermophila':
             motif = 'YYVTLD'
             if foundMotif == 0:
                if fnmatch.fnmatch(subSeqA, motif) and startPoint - ind > 25:
                     AMotifLoc.append([len(seq), ind, startPoint, subSeqA, seq[ind - 1:startPoint + 1]])
                     IFDLens.append(startPoint-ind)
                     seqLens.append(len(seq)) 
                     foundMotif += 1 
        elif 'Apis' in name:
             motif = 'YFVCCD'
             if foundMotif == 0:
                if fnmatch.fnmatch(subSeqA, motif) and startPoint - ind > 25:
                     AMotifLoc.append([len(seq), ind, startPoint, subSeqA, seq[ind - 1:startPoint + 1]])
                     IFDLens.append(startPoint-ind)
                     seqLens.append(len(seq)) 
                     foundMotif += 1 
        elif name == 'Crassostrea virginica':
            motif = 'SKGVKD'
            if foundMotif == 0:
                if fnmatch.fnmatch(subSeqA, motif) and startPoint - ind > 25:
                    AMotifLoc.append([len(seq), ind, startPoint, subSeqA, seq[ind - 1:startPoint + 1]])
                    IFDLens.append(startPoint-ind)
                    seqLens.append(len(seq)) 
                    foundMotif += 1 
        else:
            for motif in AMotifs:
                if foundMotif == 0:
                    if fnmatch.fnmatch(subSeqA, motif) and startPoint - ind > 40:
                        AMotifLoc.append([len(seq), ind, startPoint, subSeqA, seq[ind - 1:startPoint + 1]])
                        IFDLens.append(startPoint-ind)
                        seqLens.append(len(seq)) 
                        foundMotif += 1

    print(AMotifLoc)
    if foundMotif > 0:
        IFDSeqs[name] = AMotifLoc
    else:
        notFound.append(name+'\n')

# Create histograms to show length distributions of IFD and TERT
fig, ax = pl.subplots(nrows=1, ncols=1)
n, bins, patches = ax.hist(x=seqLens, bins='auto', color='#0504aa',
        alpha=0.7, rwidth=0.85)
ax.grid(axis='y', alpha=0.75)
pl.xlabel('Sequence Length')
pl.ylabel('Frequency')
maxfreq = n.max()
pl.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
fig.savefig('TertHistogram.png',dpi=300)

fig, ax = pl.subplots(nrows=1, ncols=1)
n, bins, patches = ax.hist(x=IFDLens, bins='auto', color='#0504aa',
        alpha=0.7, rwidth=0.85)
ax.grid(axis='y', alpha=0.75)
pl.xlabel('Sequence Length')
pl.ylabel('Frequency')
maxfreq = n.max()
pl.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
fig.savefig('IFDTrapHistogram.png',dpi=300)

fig, ax = pl.subplots(nrows=1, ncols=1)
ax.scatter(seqLens, IFDLens,marker='o')
pl.xlabel('TERT sequence length')
pl.ylabel('IFD sequence length')
fig.savefig('IFDLengthVsSequenceLength.png',dpi=300)
# Finally print out a file with the Locations of the IFD sequences

with open('IFDTrapLocation.txt', 'w') as Bp:
    for key in IFDSeqs:
        Bp.write(key + '\t')
        Bp.write(str(IFDSeqs[key][0][0:3]).strip('[').strip(']') + '\n')

with open('IFDTrapSequences.fasta', 'w') as output:
    for key in IFDSeqs:
        output.write('> ' + key.replace(' ', '_') + ' ' +
                     str(IFDSeqs[key][0][0:3]).strip('[').strip(']') + '\n')
        output.write(IFDSeqs[key][0][4] + '\n')

with open('TertFirst300Residues.fasta','w') as outFile:
    for name in organismNames:
        outFile.write('>' + name.replace(' ','_') + '\n')
        outFile.write(organismNames[name][0][2][0:301] + '\n')

with open('IFDnotFound.txt','w') as nf:
    for name in notFound:
        nf.write(name)
