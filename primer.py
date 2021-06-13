from Bio.Seq import Seq #sudo pip3 install biopython
from Bio.SeqUtils import GC
from Bio.SeqUtils import MeltingTemp as mt

#parameters 
name = "A. vaginae" #for file name
gc = 50 #GC > X
temp = 60; tempvar = 5 #temp = X +/- Y 
minl = 18; maxl = 22 #range of primer length 
maxdist = 90 #distance away from 16S
#can change these accordingly if not desired results 

#sequence to find primers for 
seq = "GATGAACGCTGGCGGCGCGCCTAACACATGCAAGTCGAACGGTTAAAGCATCTTCGGATGTGTATAAAGTGGCGAACGGCTGAGTAACACGTGGGCAACCTGCCCTTTGCACTGGGATAGCCTCGGGAAACCGAGGTTAATACCGGATACTCCATATTTGTCGCATGGCGAATATGGGAAAGCTCCGGCGGCAAAGGATGGGCCCGCGGCCTGTTAGCTAGTTGGTGGGGTAGTGGCCTACCAAGGCAATGATGGGTAGCCGGGTTGAGAGACCGACCGGCCAGATTGGGACTGAGACACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATCTTGCACAATGGGCGAAAGCCTGATGCAGCGACGCCGCGTGCGGGATGAAGGCCTTCGGGTTGTAAACCGCTTTCAGCAGGGACGAGGCCGCAAGGTGACGGTACCTGCAGAAGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCAAGCGTTATCCGGATTCATTGGGCGTAAAGCGCGCGTAGGCGGTCTGTTAGGTCAGGAGTTAAATCTGGGGGCTCAACCCCTATCCGCTCCTGATACCGGCAGGCTTGAGTCTGGTAGGGGAAGATGGAATTCCAAGTGTAGCGGTGAAATGCGCAGATATTTGGAAGAACACCGGTGGCGAAGGCGGTCTTCTGGGCCATGACTGACGCTGAGGCGCGAAAGCTAGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCTAGCTGTAAACGATGGACACTAGGTGTGGGGAGATTATACTTTCCGTGCCGCAGCTAACGCATTAAGTGTCCCGCCTGGGGAGTACGGTCGCAAGACTAAAACTCAAAGGAATTGACGGGGGCCCGCACAAGCAGCGGAGCATGTGGCTTAATTCGAAGCAACGCGAAGAACCTTACCAGGGCTTGACATTTAGGTGAAGCAGTGGAAACACTGTGGCCGAAAGGAGCCTAAACAGGTGGTGCATGGCTGTCGTCAGCTCGTGTCGTGAGATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTGTCGCATGTTGCCAGCGGTTCGGCCGGGCACCCATGCGAGACCGCCGGCGTTAAGCCGGAGGAAGGTGGGGACGACGTCAAGTCATCATGCCCCTTATGTCCTGGGCTGCACACGTGCTACAATGGCCGGCACAGAGGGCTGCTACTGCGCGAGCAGAAGCGAATCCCTAAAGCCGGTCCCAGTTCGGATTGGAGGCTGCAACTCGCCTCCATGAAGTCGGAGTTGCTAGTAATCGCGGATCAGCACGCCGCGGTGAATGCGTTCCCGGGCCTTGTACACACCGCCCGTCACACCACCCGAGTCGTCTGCACCCGAAGTCGTCGGCCTAACCCGCAAGGGAGGGAGGCGCCGAAGGTGTGGAGGGTAAGGGGGGT"
F16s = "AAACTCAAAGGAATTGACGG" #starting seq for variable region
R16s = "GTCAGCTCGTGTCGTGAG" #ending seq for variable region

def findTemp(x):
	y = mt.Tm_NN(x, dnac1=300, dnac2=300) # nearst neighbor with Primer concentration (nm)
	return round(y, 2)

def GCr(x):
	return round(GC(x),2)

F16si = seq.index(F16s) 
# print (F16si+1, F16si+len(F16s))
# exit()

R16si = seq.index(R16s) 

forward = Seq(seq[0:F16si])
reverse = Seq(seq[R16si+len(R16s):])

# print(reverse.reverse_complement())


# test= Seq("GGCAACATGCGACAAGGGTT")
# print (GC(test))
# print('%0.2f' % findTemp(test))

#starting with forward 
forwardT = forward[:int(len(forward)-30)] #30 bp gap bt 16s and this primer 
# print(str(forwardT))

Flist = [] #"start","end", "seq", "gc", "temp"

partsF = [forwardT[i:i+j] for i in range(len(forwardT)-minl) for j in range(minl,maxl+1)]
partsF = list( dict.fromkeys(partsF) )
#print(len(partsF))

Fgc = []
for part in partsF:
	partsFtr = str(part)
	if GC(part) > gc and findTemp(part) <= temp + tempvar and findTemp(part) >= temp - tempvar and seq.index(partsFtr)+1 > F16si - maxdist:
		# print (seq.index(partsFtr)+1, seq.index(partsFtr)+len(partsFtr), partsFtr, GCr(partsFtr), findTemp(partsFtr))
		Fgc.append((seq.index(partsFtr)+1, seq.index(partsFtr)+len(partsFtr), partsFtr, GCr(partsFtr), findTemp(partsFtr))) 


# for x in Fgc:
# 	print (x)
# print (len(Fgc))
# print (seq.index(part)+1, seq.index(part)+len(part), part, GC(part), findTemp(part))


#doing reverse
reverseT = reverse[30:]

partsR = [reverseT[i:i+j] for i in range(len(reverseT)-minl) for j in range(minl,maxl+1)]
partsR = list( dict.fromkeys(partsR) )

# print(len(partsR))
Rgc = []
for part in partsR:
	partsRstr = str(part)
	part_rc = part.reverse_complement()
	part_rcstr = str(part_rc)
	if GC(part_rc) > gc and findTemp(part_rc) <= temp + tempvar and findTemp(part_rc) >= temp - tempvar and seq.index(partsRstr)+1 < R16si + maxdist:
		# print (seq.index(partsRstr)+1, seq.index(partsRstr)+len(partsRstr), partsRstr, GCr(partsRstr), findTemp(partsRstr))
		Rgc.append((seq.index(partsRstr)+1, seq.index(partsRstr)+len(partsRstr), part_rcstr, GCr(part_rcstr), findTemp(part_rcstr))) 

# for x in Rgc:
# 	print (x)
# print (len(Rgc))

i = 1
f = open(name+" F primer.txt", "w")
for x in Fgc[::-1]:
	f.write(str(">" + str(i) + ". F"+ str(x[0])+"-"+str(x[1])+" GC = "+ str(x[3])+ " tm = "+ str(x[4])+"\n"))
	f.write(x[2]+"\n")
	i +=1 
f.close()

i = 1
f = open(name+" R primer.txt", "w")
for x in Rgc:
	f.write(str(">" + str(i) + ". R"+ str(x[0])+"-"+str(x[1])+" GC = "+ str(x[3])+ " tm = "+ str(x[4])+"\n"))
	f.write(x[2]+"\n")
	i +=1 
f.close()






