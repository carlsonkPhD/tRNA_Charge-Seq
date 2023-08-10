import numpy as np
#initialize lists
gene = []
seq = []

count = 0
with open("mm10-tRNAs.fa") as fp:
	for line in fp:
		count += 1
		if count % 3 == 1:
			#adds gene name to gene list
			namefull = line.strip().split(' ', 1)[0]
			name = namefull[1:]
			gene.append(name)
		if count % 3 == 2:
			#add first line of sequence to gene name
			seq.append(line.strip())
		if count % 3 == 0:
			#add 2nd line of sequence to first
			newidx = (count-3)/3
			seq[int(newidx)] = seq[int(newidx)]+line.strip()


#combine gene and seq lists into an array of tuples
myGenes = np.column_stack((gene, seq))
#create empty list to hold unique sequences 
uniqueGenes = []

#identify unique sequences in myGenes array then add the sequence and tRNA ID to uniqueGenes list
for x in myGenes:
    result = next((k for k, v in enumerate(uniqueGenes) if v[1] == x[1]), -1)
    if result == -1:
        uniqueGenes.append(list(x))
    else:
        continue
    
    
f = open("mm10-tRNAsCCA.fa", "a")
for t in uniqueGenes:
    #add "CCA" to each sequence
    t[1] = t[1]+'CCA'
    f.write('>'+ t[0] + '\n' + t[1] + '\n')
f.close()