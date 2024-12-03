import Levenshtein 

haplos_id = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10']
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
print(sorted_edits[0][0])

haplo3 = haplos[sorted_edits[1][0][0]]
haplo4 = haplos[sorted_edits[1][0][1]]
print(sorted_edits[1][0])

haplo5 = haplos[sorted_edits[2][0][0]]
haplo6 = haplos[sorted_edits[2][0][1]]
print(sorted_edits[2][0])

rec_haplo1 = haplo1[:int(len(haplo1)/3)] + haplo2[int(len(haplo2)/3):int(len(haplo2)/3)*2] + haplo1[int(len(haplo1)/3)*2:]
rec_haplo2 = haplo3[:int(len(haplo3)/3)] + haplo4[int(len(haplo4)/3):int(len(haplo4)/3)*2] + haplo3[int(len(haplo3)/3)*2:]
rec_haplo3 = haplo5[:int(len(haplo5)/3)] + haplo6[int(len(haplo6)/3):int(len(haplo6)/3)*2] + haplo5[int(len(haplo5)/3)*2:]

with open ('data/sars-cov-2/rec_haplos/sequence_r01.fasta', 'w') as f:
    f.write('>rec_17\n')
    f.write(rec_haplo1)

with open ('data/sars-cov-2/rec_haplos/sequence_r02.fasta', 'w') as f:
    f.write('>rec_19\n')
    f.write(rec_haplo2)
    
with open ('data/sars-cov-2/rec_haplos/sequence_r03.fasta', 'w') as f:
    f.write('>rec_14\n')
    f.write(rec_haplo3)