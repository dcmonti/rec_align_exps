import Levenshtein 
import os
haplos_id = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20']
haplos = []

for i in haplos_id:
    with open ('data/sars-cov-2/haplotypes/sequence_' + i + '.fasta', 'r') as f:
        haplo = ''
        for line in f:
            if line[0] != '>':
                haplo += line.strip()
        haplos.append(haplo)

edits = dict()
for i in range(len(haplos)):
    for j in range(i+1, len(haplos)):
        distance = Levenshtein.distance(haplos[i], haplos[j])
        edits[(i,j)] = distance
sorted_edits = sorted(edits.items(), key=lambda x: x[1])
sorted_edits.reverse()
print(sorted_edits)
haplo1 = haplos[sorted_edits[0][0][0]]
haplo2 = haplos[sorted_edits[0][0][1]]
id1, id2 = sorted_edits[0][0]

haplo3 = haplos[sorted_edits[1][0][0]]
haplo4 = haplos[sorted_edits[1][0][1]]
id3, id4 = sorted_edits[1][0]

haplo5 = haplos[sorted_edits[2][0][0]]
haplo6 = haplos[sorted_edits[2][0][1]]
id5, id6 = sorted_edits[2][0]

rec_haplo1 = haplo1[:int(len(haplo1)/2)] + haplo2[int(len(haplo2)/2):]
rec_haplo2 = haplo3[:int(len(haplo3)/2)] + haplo4[int(len(haplo4)/2):]
rec_haplo3 = haplo5[:int(len(haplo5)/2)] + haplo6[int(len(haplo6)/2):]

with open ('data/sars-cov-2/rec_haplos/sequence_r01.fasta', 'w') as f:
    f.write(f'>rec_{id1}_{id2}\n')
    f.write(rec_haplo1)

with open ('data/sars-cov-2/rec_haplos/sequence_r02.fasta', 'w') as f:
    f.write(f'>rec_{id3}_{id4}\n')
    f.write(rec_haplo2)
    
with open ('data/sars-cov-2/rec_haplos/sequence_r03.fasta', 'w') as f:
    f.write(f'>rec_{id5}_{id6}\n')
    f.write(rec_haplo3)