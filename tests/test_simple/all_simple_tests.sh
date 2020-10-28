python3 tests/test_simple/generate_simple_PE_fastq.py

# verify that all reads are assigned correctly, except long reads which have no match. UMI should be same in each case
# also verify that reads are correct length
python3 ultraplex.py -i tests/test_simple/reads1.fastq.gz -b tests/test_simple/barcodes_5_and_3.csv -d tests/test_simple/SE -o single_end

# verify that reads are renamed correctly
python3 ultraplex.py -i tests/test_simple/reads1.fastq.gz -b tests/test_simple/barcodes_5_and_3_named.csv -d tests/test_simple/SE_named -o single_end_named

# verify that all reads are the correct length, that reverse reads are TTT..., that UMI is the same for both reads, that ALL reads are assigned
python3 ultraplex.py -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3.csv -d tests/test_simple/PE -o paired_end

# verify that all reads are the correct length, that reverse reads are TTT..., that UMI is the same for both reads, that ALL reads are assigned
python3 ultraplex.py -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3_named.csv -d tests/test_simple/PE_named -o paired_end_named
