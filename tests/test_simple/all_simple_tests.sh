python3 tests/test_simple/generate_simple_PE_fastq.py

# verify that all reads are assigned correctly, except long reads which have no match. UMI should be same in each case
python3 ultraplex.py -i tests/test_simple/reads1.fastq.gz -b tests/test_simple/barcodes_5_and_3.csv -d tests/test_simple/SE -o single_end

# verify that reads are renamed correctly
python3 ultraplex.py -i tests/test_simple/reads1.fastq.gz -b tests/test_simple/barcodes_5_and_3_named.csv -d tests/test_simple/SE_named -o single_end_named
