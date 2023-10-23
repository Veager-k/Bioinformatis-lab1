from Bio.Seq import Seq
import math

    

bac1 = open("data\\bacterial1.fasta").read()
bac2 = open("data\\bacterial2.fasta").read()
bac3 = open("data\\bacterial3.fasta").read()
bac4 = open("data\\bacterial4.fasta").read()

mam1 = open("data\\mamalian1.fasta").read()
mam2 = open("data\\mamalian2.fasta").read()
mam3 = open("data\\mamalian3.fasta").read()
mam4 = open("data\\mamalian4.fasta").read()


sequences = [bac1, bac2, bac3, bac4, mam1, mam2, mam3, mam4]


#remove prefix and '\n'
for i, seq in enumerate(sequences):
    seq = seq[seq.find("\n")+1:]
    sequences[i] = seq.replace("\n", "")


#finding pairs
pairs = []
    
for v, virus in enumerate(sequences):
    pairs.append([[], []])
    pairs[v][0].append([0, 0])
    pairs[v][1].append([0, 0])
    complement = str(Seq(virus).reverse_complement())
    
    for s in range(len(virus)):
        if virus[s:s+3] == "ATG":
            for e in range(s, len(virus), 3):
                if virus[e:e+3] == "TGA" or virus[e:e+3] == "TAG" or virus[e:e+3] == "TAA":
                    if e-s < 100: #check if its not too short
                        break
                    if pairs[v][0][-1][1] == e: #check if end is the same as prev element
                        break
                    pairs[v][0].append([s, e])
                    break
    
    for s in range(len(complement)):
        if complement[s:s+3] == "ATG":
            for e in range(s, len(complement), 3):
                if complement[e:e+3] == "TGA" or complement[e:e+3] == "TAG" or complement[e:e+3] == "TAA":
                    if e-s < 100: #check if its not too short
                        break
                    if pairs[v][1][-1][1] == e: #check if end is the same as prev element
                        break
                    pairs[v][1].append([s, e])
                    break
    
    pairs[v][0].pop(0)
    pairs[v][1].pop(0)
    
    
#convert start-end pairs to sequences
seq = []
    
for v, virus in enumerate(pairs):
    seq.append([])
    
    for pair in virus[0]:
        seq[v].append(sequences[v][pair[0]+3:pair[1]])
    
    for pair in virus[1]:
        complement = str(Seq(sequences[v]).reverse_complement())
        seq[v].append(complement[pair[0]+3:pair[1]])


#translate to acids
for i, virus in enumerate(seq):
    for j, sequence in enumerate(seq[i]):
        while len(sequence) % 3 != 0:
           sequence = sequence + 'N'
        seq[i][j] = str(Seq(sequence).translate())


#form frequency lists
acidList = ['A', 'R', 'N', 'D', 'C', 'Q', 'E',
         'G', 'H', 'I', 'L', 'K', 'M', 'F',
         'P', 'S', 'T', 'W', 'Y', 'V'
        ]

acidListDI = [i+j for i in acidList for j in acidList]


frequencies = []
for i in range(8):
    frequencies.append([])
    for acid in acidList:
        frequencies[i].append([acid, 0])
        
frequenciesDI = []
for i in range(8):
    frequenciesDI.append([])
    for acid in acidListDI:
        frequenciesDI[i].append([acid, 0])


#find frequencies
for i, virus in enumerate(seq):
    for sequence in virus:
        for acid in sequence:
            for ac in frequencies[i]:
                if acid == ac[0]:
                    ac[1] += 1
                    break;

for virus in frequencies:
    sumTotal = 0
    
    for ac in virus:
        sumTotal += ac[1]
        
    for ac in virus:
        ac[1] = (ac[1] / sumTotal)*100
    
                    
for i, virus in enumerate(seq):
    for sequence in virus:
        for a in range(0, len(sequence), 2):
            if a+1 == len(sequence):
                break;
            for acids in frequenciesDI[i]:
                if (sequence[a] + sequence[a+1]) == acids[0]:
                    acids[1] += 1
                    break;

for virus in frequenciesDI:
    sumTotal = 0
    
    for ac in virus:
        sumTotal += ac[1]
        
    for ac in virus:
        ac[1] = (ac[1] / sumTotal)*100


#get average distances between viruses
distances = []
distancesDI = []

for v, virus1 in enumerate(frequencies):
    distances.append([])
    for virus2 in frequencies:
        sumTotal = 0
        
        for i in range(len(virus1)):
            if (virus1[i][1] + virus2[i][1]) == 0:
                continue;
            sumTotal += abs(virus1[i][1] - virus2[i][1]) / (virus1[i][1] + virus2[i][1])
            #sumTotal += abs(math.log(virus1[i][1] / virus2[i][1]))
        
        sumTotal /= len(virus1)
        distances[v].append(sumTotal)


for v, virus1 in enumerate(frequenciesDI):
    distancesDI.append([])
    for virus2 in frequenciesDI:
        sumTotal = 0
        
        for i in range(len(virus1)):
            if (virus1[i][1] + virus2[i][1]) == 0:
                continue;
            sumTotal += abs(virus1[i][1] - virus2[i][1]) / (virus1[i][1] + virus2[i][1])
            #sumTotal += abs(math.log(virus1[i][1] / virus2[i][1]))
        
        sumTotal /= len(virus1)
        distancesDI[v].append(sumTotal)



#find acids with most variation
distancesByAcid = []

for a, acid in enumerate(acidList):
    sumTotal = 0
    for v, virus1 in enumerate(frequencies):
        for virus2 in frequencies[v:]:
            if (virus1[a][1] + virus2[a][1]) == 0:
                continue;
            sumTotal += abs(virus1[a][1] - virus2[a][1]) / (virus1[a][1] + virus2[a][1])
    sumTotal /= 36
    distancesByAcid.append([acid, sumTotal])


distancesByAcidDI = []

for a, acidDI in enumerate(acidListDI):
    sumTotal = 0
    for v, virus1 in enumerate(frequenciesDI):
        for virus2 in frequenciesDI[v:]:
            if (virus1[a][1] + virus2[a][1]) == 0:
                continue;
            sumTotal += abs(virus1[a][1] - virus2[a][1]) / (virus1[a][1] + virus2[a][1])
    sumTotal /= 36
    distancesByAcidDI.append([acidDI, sumTotal])

distancesByAcid = sorted(distancesByAcid, key=lambda x: x[1], reverse=True)
distancesByAcidDI  = sorted(distancesByAcidDI, key=lambda x: x[1], reverse=True)

for i in range(5):
    print(distancesByAcid[i])
print()

for i in range(5):
    print(distancesByAcidDI[i])
print()


#print matrices
viruses = ["bac1", "bac2", "bac3", "bac4",
           "mam1", "mam2", "mam3", "mam4"]

print(8)
for v, virus in enumerate(distances):
    print(viruses[v] + " ", end="")
    for num in virus:
        print(f"{num:4.3f} ", end="")
    print()
    

print(8)
for v, virus in enumerate(distancesDI):
    print(viruses[v] + " ", end="")
    for num in virus:
        print(f"{num:4.3f} ", end="")
    print()
    
    