python3 tests/test_simple/generate_simple_PE_fastq.py

python3 ultraplex.py -i tests/test_simple/reads1.fastq.gz -b tests/test_simple/barcodes_5_and_3.csv -d tests/test_simple/SE -o single_end
