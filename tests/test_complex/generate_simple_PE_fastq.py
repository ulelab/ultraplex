import math
import random
import gzip

def rev_c(seq):
	# simple function that reverse complements a given sequence
	
	# first reverse the sequence
	seq = seq[::-1]
	# and then complement
	d = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N"}
	seq2 = ""
	for element in seq:
		seq2 = seq2 + d[element]

	return seq2

def make_rand_seq(l):
	d = {1:'A', 2:'C', 3:'G', 4:'T'}
	seq = ""
	for _ in range(l):
		seq = seq + d[random.randint(1,4)]
	return seq

def quality_string(l, chance):
	# make list
	qualities = "#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMN"
	index = len(qualities)-1

	q = ""
	for _ in range(l):
		q = q + qualities[index]
		if random.randint(0,chance) == 0 and index > 0:
			index = index-1
	return q

def write_fastq_gz(reads, file_name):
	with gzip.open(file_name, "wb") as file:
		for i, read in enumerate(reads):
			file.write(("@"+str(i)+'\n').encode())
			file.write((read+'\n').encode())
			file.write(b'+\n')
			file.write((quality_string(len(read), chance = 6)+'\n').encode())


# Main script

n_seqs = 5000
min_length = 20
max_length = 70
five_p_bcs = ['AAAATGCC', 'AAACCTCC', 'AAAATACC']
three_p_bcs = ['CC', 'AT', 'TG']
#five_p_bcs = ['AAAATGCC', ]
#three_p_bcs = ["NNAT"]
adapter1 = "AGATCGGAAGAGCGGTTCAG"
adapter2 = "AGATCGGAAGAGCGTCGTG"
max_read_length = 100

reads1 = []
reads2 = []

for i in range(n_seqs):
	# first, generate a random sequence
	seq = make_rand_seq(random.randint(min_length, max_length))

	# now add barcodes
	seq = five_p_bcs[random.randint(0,len(five_p_bcs)-1)] + seq + three_p_bcs[random.randint(0,len(three_p_bcs)-1)]

	read1 = seq + adapter1
	read2 = rev_c(seq) + adapter2

	if len(read1) > max_read_length:
		read1 = read1[0:max_read_length]
	if len(read2) > max_read_length:
		read2 = read2[0:max_read_length]

	reads1.append(read1)
	reads2.append(read2)

write_fastq_gz(reads1, "reads1.fastq.gz")
write_fastq_gz(reads2, "reads2.fastq.gz")