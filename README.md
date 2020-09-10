# ultraplex
Ultra-fast 5' and 3' demultiplexer. Able to quality-trim, adapter-trim and demultiplex a whole lane in 30 minutes.

# Usage
Create a csv of barcodes. The first column should be the five prime barcodes. The second column (optional) should be the 3' barcodes associated with that 5' barcode. Each 3' barcode should be separated by a semi-colon. There should be no other characters except: "A", "C", "T", "G", "N" (for UMIs), ",", and ";"

For example:

NNNNCCTGCNNNNNN
NNNNGGAGCNNNNNN,NNATT;NNAGG;NNCTG
NNNNAGTCCNNNNNN,NNATT;NNAGG;NNCTG
NNNNTCGCCNNNNNN,NNATT;NNAGG;NNCTG
NNNNCATACNNNNNN,NNATT;NNAGG
NNNNGTAACNNNNNN,NNATT;NNAGG;NNCTG
