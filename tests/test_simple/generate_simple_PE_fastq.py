import math
import random
import gzip


def rev_c(seq):
    # simple function that reverse complements a given sequence

    # first reverse the sequence
    seq = seq[::-1]
    # and then complement
    d = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    seq2 = ""
    for element in seq:
        seq2 = seq2 + d[element]

    return seq2


def make_rand_seq(l):
    d = {1: "A", 2: "C", 3: "G", 4: "T"}
    seq = ""
    for _ in range(l):
        seq = seq + d[random.randint(1, 4)]
    return seq


def quality_string(l, chance):
    # make list
    qualities = "#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMN"
    index = len(qualities) - 1

    q = ""
    for _ in range(l):
        q = q + qualities[index]
        if random.randint(0, chance) == 0 and index > 0:
            index = index - 1
    return q


def write_fastq_gz(reads, seqnames, file_name):
    with gzip.open(file_name, "wb") as file:
        for read, seqname in zip(reads, seqnames):
            file.write(("@" + seqname + "\n").encode())
            file.write((read + "\n").encode())
            file.write(b"+\n")
            file.write((quality_string(len(read), chance=6) + "\n").encode())


# Main script

n_seqs = 300
five_p_bcs = ["TTATGTT", "TTCCATT", "TTGCGTT"]  # form NNBBBNN
three_p_bcs = ["GGCC", "GGAT", "GGTG"]  # form NNBB
adapter1 = "AGATCGGAAGAGCGGTTCAG"
adapter2 = "AGATCGGAAGAGCGTCGTG"
max_read_length = 100

reads1 = []
reads2 = []
seqnames = []

for i in range(n_seqs):
    for five_p_bc in five_p_bcs:
        for three_p_bc in three_p_bcs:
            # first, generate a random sequence
            seq = "A" * i  # insert has form AAAAA....
            seq = five_p_bc + seq + three_p_bc

            read1 = seq + adapter1
            read2 = rev_c(seq) + adapter2

            if len(read1) > max_read_length:
                read1 = read1[0:max_read_length]
            if len(read2) > max_read_length:
                read2 = read2[0:max_read_length]

            seqname = str(i) + "_" + five_p_bc + "_" + three_p_bc

            reads1.append(read1)
            reads2.append(read2)
            seqnames.append(seqname)

write_fastq_gz(reads1, seqnames, "tests/test_simple/reads1.fastq.gz")
write_fastq_gz(reads2, seqnames, "tests/test_simple/reads2.fastq.gz")
