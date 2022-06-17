import sys

# This program calculates dinucleotides of nucleosomal DNA.  For single-end reads it extraplolates the nucleosomal DNA by using the leftmost base position
# as a reference.  For paired-end reads it extrapolates a dyad by taking the postion halfway between the two reads.  Once a dyad is calculated, the
# program goes 5' 74 bp and 3' 73 bp.


SAM = open(sys.argv[1],'r') #SAM input file
REF = open(sys.argv[2],'r') #reference genomce sequence
output = str(sys.argv[3]) #name of output file
f = open(output,'w')

#break reference genome into strings by chromosome
chrI = ""
chrII = ""
chrIII = ""
chrIV = ""
chrV = ""
chrX = ""
chrMtDNA = ""

chromoFlag = ""
for line in REF:
    #determine chromoFlag change
    if line.startswith('>CHROMOSOME_IV'):
        chromoFlag = '4'
        continue
    elif line.startswith('>CHROMOSOME_III'):
        chromoFlag = '3'
        continue
    elif line.startswith('>CHROMOSOME_II'):
        chromoFlag = '2'
        continue
    elif line.startswith('>CHROMOSOME_I'):
        chromoFlag = '1'
        continue
    elif line.startswith('>CHROMOSOME_V'):
        chromoFlag = '5'
        continue
    if line.startswith('>CHROMOSOME_X'):
        chromoFlag = 'X'
        continue
    elif line.startswith('>CHROMOSOME_MtDNA'):
        chromoFlag = 'M'
        continue

    if chromoFlag == '4':
        chrIV += str(line)
    elif chromoFlag == '3':
        chrIII += str(line)
    elif chromoFlag == '2':
        chrII += str(line)
    elif chromoFlag == '1':
        chrI += str(line)
    elif chromoFlag == '5':
        chrV += str(line)
    elif chromoFlag == 'X':
        chrX += str(line)
    elif chromoFlag == 'M':
        chrMtDNA += str(line)

#create list with nucleosomal DNA sequences
DNA = []

#take SAM file data and create a list of nucleosomal DNA sequences
for line in SAM:
    if line.startswith('@'):
        continue
    else:
        lineList = line.split('\t')
        if lineList[1] == '16':
            left = int(lineList[3]) - 151 + len(lineList[9])
            right= int(lineList[3]) + len(lineList[9]) - 4
            if int(lineList[3]) >= 150 - len(lineList[9]):
                if lineList[2] == "I":
                    DNA.append(chrI[left:right])
                elif lineList[2] == "II":
                    DNA.append(chrII[left:right])
                elif lineList[2] == "III":
                    DNA.append(chrIII[left:right])
                elif lineList[2] == "IV":
                    DNA.append(chrIV[left:right])
                elif lineList[2] == "V":
                    DNA.append(chrV[left:right])
                elif lineList[2] == "X":
                    DNA.append(chrX[left:right])
                elif lineList[2] == "MtDNA":
                    DNA.append(chrMtDNA[left:right])
                else:
                    continue
            else:
                continue
        elif lineList[1] == "0":
            left = int(lineList[3]) - 1
            right = int(lineList[3]) + 146
            if lineList[2] == "I":
                if int(lineList[3]) <= 15072284:
                    DNA.append(chrI[left:right])
                else:
                    continue
            elif lineList[2] == "II":
                if int(lineList[3]) <= 15279271:
                    DNA.append(chrII[left:right])
                else:
                    continue
            elif lineList[2] == "III":
                if int(lineList[3]) <= 13783651:
                    DNA.append(chrIII[left:right])
                else:
                    continue
            elif lineList[2] == "IV":
                if int(lineList[3]) <= 17493679:
                    DNA.append(chrIV[left:right])
                else:
                    continue
            elif lineList[2] == "V":
                if int(lineList[3]) <= 20924030:
                    DNA.append(chrV[left:right])
                else:
                    continue
            elif lineList[2] == "X":
                if int(lineList[3]) <= 17718792:
                    DNA.append(chrX[left:right])
                else:
                    continue
            elif lineList[2] == "MtDNA":
                if int(lineList[3]) <= 13644:
                    DNA.append(chrMtDNA[left:right])
                else:
                    continue
        elif lineList[1] == '99' or '163':
            dyad = int(int(lineList[8]) // 2) + int(lineList[3]) -1
            left = dyad - 74
            right = dyad + 73
            if left <= 1:
                continue
            elif lineList[2] == "I":
                if right >= 15072434:
                    continue
                else:
                    DNA.append(chrI[left:right])
            elif lineList[2] == "II":
                if right >= 15279421:
                    continue
                else:
                    DNA.append(chrII[left:right])
            elif lineList[2] == "III":
                if right >= 13783801:
                    continue
                else:
                    DNA.append(chrIII[left:right])
            elif lineList[2] == "IV":
                if right >= 17493829:
                    continue
                else:
                    DNA.append(chrIV[left:right])
            elif lineList[2] == "V":
                if right >= 20924180:
                    continue
                else:
                    DNA.append(chrV[left:right])
            elif lineList[2] == "X":
                if right >= 15279421:
                    continue
                else:
                    DNA.append(chrX[left:right])
            elif lineList[2] == "MtDNA":
                if right >= 13794:
                    continue
                else:
                    DNA.append(chrMtDNA[left:right])
            else:
                continue
        else:
            continue

#count the bases for each position in the nucleosomal DNA sequence file
aaCountByPosition = [0] * 147
atCountByPosition = [0] * 147
acCountByPosition = [0] * 147
agCountByPosition = [0] * 147
caCountByPosition = [0] * 147
ctCountByPosition = [0] * 147
ccCountByPosition = [0] * 147
cgCountByPosition = [0] * 147
gaCountByPosition = [0] * 147
gtCountByPosition = [0] * 147
gcCountByPosition = [0] * 147
ggCountByPosition = [0] * 147
taCountByPosition = [0] * 147
ttCountByPosition = [0] * 147
tcCountByPosition = [0] * 147
tgCountByPosition = [0] * 147
for sequence in DNA:
    if len(sequence) == 0:
        continue
    for index in range(0, 146):
        base = sequence[index]
        nextBase = sequence[index + 1]
        if base == 'a' and nextBase == 'a':
            aaCountByPosition[index] += 1
        if base == 'a' and nextBase == 't':
            atCountByPosition[index] += 1
        if base == 'a' and nextBase == 'c':
            acCountByPosition[index] += 1
        if base == 'a' and nextBase == 'g':
            agCountByPosition[index] += 1
        if base == 'c' and nextBase == 'a':
            caCountByPosition[index] += 1
        if base == 'c' and nextBase == 't':
            ctCountByPosition[index] += 1
        if base == 'c' and nextBase == 'c':
            ccCountByPosition[index] += 1
        if base == 'c' and nextBase == 'g':
            cgCountByPosition[index] += 1
        if base == 'g' and nextBase == 'a':
            gaCountByPosition[index] += 1
        if base == 'g' and nextBase == 't':
            gtCountByPosition[index] += 1
        if base == 'g' and nextBase == 'c':
            gcCountByPosition[index] += 1
        if base == 'g' and nextBase == 'g':
            ggCountByPosition[index] += 1
        if base == 't' and nextBase == 'a':
            taCountByPosition[index] += 1
        if base == 't' and nextBase == 't':
            ttCountByPosition[index] += 1
        if base == 't' and nextBase == 'c':
            tcCountByPosition[index] += 1
        if base == 't' and nextBase == 'g':
            tgCountByPosition[index] += 1

for position in range(0,147):
    f.write('Position: ' + str(position + 1) + '\n' +
    'AA: ,' + str(aaCountByPosition[position]) + ',' + '\n' +
    'AC: ,' + str(acCountByPosition[position]) + ',' + '\n' +
    'AG: ,' + str(agCountByPosition[position]) + ',' + '\n' +
    'AT: ,' + str(atCountByPosition[position]) + ',' + '\n' +
    'CA: ,' + str(caCountByPosition[position]) + ',' + '\n' +
    'CC: ,' + str(ccCountByPosition[position]) + ',' + '\n' +
    'CG: ,' + str(cgCountByPosition[position]) + ',' + '\n' +
    'CT: ,' + str(ctCountByPosition[position]) + ',' + '\n' +
    'GA: ,' + str(gaCountByPosition[position]) + ',' + '\n' +
    'GC: ,' + str(gcCountByPosition[position]) + ',' + '\n' +
    'GG: ,' + str(ggCountByPosition[position]) + ',' + '\n' +
    'GT: ,' + str(gtCountByPosition[position]) + ',' + '\n' +
    'TA: ,' + str(taCountByPosition[position]) + ',' + '\n' +
    'TC: ,' + str(tcCountByPosition[position]) + ',' + '\n' +
    'TG: ,' + str(tgCountByPosition[position]) + ',' + '\n' +
    'TT: ,' + str(ttCountByPosition[position]) + ',' + '\n')
f.close()
