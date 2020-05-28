from collections import defaultdict
import matplotlib.pyplot as pl
import numpy as np
import fnmatch

''' Script for locating the position of CP motif
   in a list of tert sequences.'''

# First load the list of motifs to search for
CPMotifs = []

for line in open('./INPUTS/new_CPMotifs.txt', 'r'):
    if line[0] != '>':
        CPMotifs.append(line.strip())

# Next load in the list of file names so that individual
# sequences can be loaded 1 by 1

FileNames = []

for f in open('FileList.txt', 'r'):
    FileNames.append(f.strip())

print(FileNames)
# Next set up parameters for doing a sliding window search
# through the sequences.

windowSize = 10
organismNames = defaultdict()

# Loop through all files, extract name and sequence and position
# of CP motif, save this in a dictionary.
noCP = []
for f in FileNames:
    seq = []
    for line in open(f, 'r'):
        if line[0] == '>':
            name = line.split('[')[1][:-2]
        else:
            seq.append(line.strip())

    fullSeq = ''.join(seq)

    motifsLoc = []
    foundMotif = 0
    for i in range(0, len(fullSeq) - windowSize):
        subSeq = fullSeq[i:i + windowSize]

        for motif in CPMotifs:
            if fnmatch.fnmatch(subSeq, motif):
                if i > 50:
                    motifsLoc.append([i, motif, fullSeq])
                    foundMotif += 1
                    print(motifsLoc[0][0:1])

    if foundMotif > 0:
        organismNames[name] = motifsLoc
    else:
        noCP.append(name)
# Write out a file containing the organisms name and the
# location of putative CP motifs within its sequence.

with open('CPMotifLocation.txt', 'w') as Bp:
    for key in organismNames:
        Bp.write(key + '\t')
        Bp.write(str(organismNames[key][0][0:2]).strip('[').strip(']') + '\n')

# Now we have a dictionary that contains both the sequence
# and the location of the CP motif we can loop backwards
# until we reach the start of the sequence.


UpstreamSeqs = defaultdict()
seqLens = []
UpstreamLens = []

for name in organismNames:
    seq = organismNames[name][0][2]
    startPoint = organismNames[name][0][0]
    UpstreamSeqs[name] = [len(seq), startPoint, seq[:startPoint + 1]]
    UpstreamLens.append(startPoint)
    seqLens.append(len(seq))

print("Longest sequence is " + str(np.max(seqLens)))
print("Shortest sequence is " + str(np.min(seqLens)))
print("Average sequence is " + str(np.mean(seqLens)))
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
n, bins, patches = ax.hist(x=UpstreamLens, bins='auto', color='#0504aa',
        alpha=0.7, rwidth=0.85)
ax.grid(axis='y', alpha=0.75)
pl.xlabel('Sequence Length')
pl.ylabel('Frequency')
maxfreq = n.max()
pl.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
fig.savefig('TenAndLinkerHistogram.png',dpi=300)

fig, ax = pl.subplots(nrows=1, ncols=1)
ax.scatter(seqLens, UpstreamLens,marker='o')
pl.xlabel('TERT sequence length')
pl.ylabel('Ten and linker sequence length')
fig.savefig('TenAndLinkerLengthVsSequenceLength.png',dpi=300)

# Finally print out a file with the sequence upstream of CP

with open('TenAndLinkerSequences.fasta', 'w') as output:
    for key in UpstreamSeqs:
        output.write('> ' + key.replace(' ', '_') + ', ' +
                     str(UpstreamSeqs[key][0:2]).strip('[').strip(']') + '\n')
        output.write(UpstreamSeqs[key][2] + '\n')
with open('noCPFound.txt','w') as outFile:
    for name in noCP:
        outFile.write(name.replace(" ", "_")+'\n')

