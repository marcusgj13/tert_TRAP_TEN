from os import listdir

H_content = []
E_content = []
HtoE_ratio = []
lengths = []
Names = []
for file in listdir('./RESULTS'):
	if file.endswith('.horiz'):
		prediction = []
		Names.append(file.split('.')[0])
		with open('./RESULTS/' + file,'r') as f:
			for line in f:
				if line[0] == 'P':	
					prediction.append(line.split(':')[1].strip())

			prediction = "".join(prediction)
			seq_len = len(prediction)
			print(seq_len)
			lengths.append(str(seq_len))
			H_count = float(prediction.count('H'))
			E_count = float(prediction.count('E')) 

			if seq_len == 0:
				H_content.append(0.0)
				E_content.append(0.0)
			else:
				H_content.append(H_count/float(seq_len))
				E_content.append(E_count/float(seq_len))


			try:
				HtoE_ratio.append(str(H_count/E_count))
			except ZeroDivisionError:
				HtoE_ratio.append(str(H_count))
					



with open('20191102_PsiPred_Results.txt','w') as r:
	r.write('Filename' + '\t')
	r.write('%Helix' + '\t')
	r.write('%Sheet' + '\t')
	r.write('Helix2sheetRatio' + '\t')
	r.write('SequenceLength' + '\n')
	lis = [Names, H_content, E_content, HtoE_ratio, lengths]
	for x in zip(*lis):
		r.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(*x)) 

	
