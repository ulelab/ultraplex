##### RUN THIS SCRIPT FROM THE MAIN ULTRAPLEX DIRECTORY #####
# eg bash tests/test_simple/all_simple_tests.sh

# generate fastq files
python3 tests/test_simple/generate_simple_PE_fastq.py

# verify that all reads are assigned correctly, except long reads which have no match. UMI should be same in each case
# also verify that reads are correct length
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -b tests/test_simple/barcodes_5_and_3.csv -d tests/test_simple/SE -o single_end

# verify that reads are renamed correctly
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -b tests/test_simple/barcodes_5_and_3_named.csv -d tests/test_simple/SE_named -o single_end_named

# verify that all reads are the correct length, that reverse reads are TTT..., that UMI is the same for both reads, that ALL reads are assigned
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3.csv -d tests/test_simple/PE -o paired_end

# verify that all reads are the correct length, that reverse reads are TTT..., that UMI is the same for both reads, that ALL reads are assigned
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3_named.csv -d tests/test_simple/PE_named -o paired_end_named

# verify that only some are renamed
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3_some_named.csv -d tests/test_simple/PE_some_named -o paired_end_some_named

# verify that samples 1-3 are NOT assigned
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_sample123_wrong_5.csv -d tests/test_simple/wrong_1_2_3 -o wrong_1_2_3 -m5 0

# verify that samples 1-3 are NOT assigned and the no matches are not written out
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_sample123_wrong_5.csv -d tests/test_simple/wrong_1_2_3_inm -o wrong_1_2_3_inm -m5 0 -inm

# verify that samples 1-3 are assigned
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_sample123_wrong_5.csv -d tests/test_simple/wrong_1_2_3_mm_allowed -o wrong_1_2_3_mm_allowed -m5 1

# verify that sample 3 is NOT assigned
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_sample3_wrong_3.csv -d tests/test_simple/wrong_3 -o wrong_3

# verify that sample 3 IS assigned
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_sample3_wrong_3.csv -d tests/test_simple/wrong_3_mm_allowed -o wrong_3_mm_allowed -m3 1

# verify that samples 7 8 and 9 are combined
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_789_combined.csv -d tests/test_simple/789_combined -o 789_combined

# verify that samples 1-3, 4-6 and 7-9 have 4-6, 5-7 and 6-8 umis respectively. Verify that sample 6 has UMI of TTTAGG, and fwd read TAAAA.....AAA, rev read TTTT...TTA
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_different_umis.csv -d tests/test_simple/umi_lengths -o umi_lengths

# verify that samples 1-3, 4-6 and 7-9 have 4-6, 5-7 and 6-8 umis respectively. Verify that sample 9 has UMI of TTTTAAG, and fwd read AAAA.....AAA, rev read TTTT...TTT
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_3nt_bc_789.csv -d tests/test_simple/3nt_789 -o 3nt_789

# verify that mate adapter is correctly removed
python3 ultraplex -i tests/test_simple/example_read_1.fastq.gz -i2 tests/test_simple/example_read_2.fastq.gz -b tests/test_simple/barcodes_3nt_bc_789.csv -d tests/test_simple/mate -o mate

# verify that length filter works correctly
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -b tests/test_simple/barcodes_5_and_3.csv -d tests/test_simple/SE_length60 -o single_end_length60 -l 60

# verify that length filter works correctly
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3.csv -d tests/test_simple/PE_length60 -o paired_end_length60 -l 60

# verify that short reads do not cause crash
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3_very_long.csv -d tests/test_simple/PE_long -o paired_end_long

# verify that short reads do not cause crash
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3_very_long_barcode.csv -d tests/test_simple/PE_long_barcode -o paired_end_long_barcode

# check it works without UMIs
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_no_umi_3.csv -d tests/test_simple/no_umi -o paired_no_umi

# Check that it works with v v long UMIs
python3 ultraplex -i tests/test_simple/reads1.fastq.gz -i2 tests/test_simple/reads2.fastq.gz -b tests/test_simple/barcodes_5_and_3_very_very_long_umi.csv -d tests/test_simple/v_v_long_umi -m5 0

