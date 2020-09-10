![Ultraplex Logo](ultraplex-full.png)
Ultra-fast 5' and 3' demultiplexer. Able to quality-trim, adapter-trim and demultiplex a whole lane in 30 minutes.

# Key Features

* Quality trimming and adaptor removal at 3' end by cutadapt
* Can move UMI to header for downstream analysis (eg with UMI-TOOLS)
* Can handle different length UMIs (eg for mixed demultiplexing of irCLIP and iiCLIP)


# Usage
Create a csv of barcodes. The first column should be the five prime barcodes. The second column (optional) should be the 3' barcodes associated with that 5' barcode. Each 3' barcode should be separated by a semi-colon. There should be no other characters except: "A", "C", "T", "G", "N" (for UMIs), ",", and ";"


