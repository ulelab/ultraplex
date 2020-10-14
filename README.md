![Ultraplex Logo](gif_logo.gif)

Ultra-fast 5' and 3' demultiplexer. Able to quality-trim, adapter-trim and demultiplex a whole lane in 30 minutes.

# Key Features

* Quality trimming and adaptor removal at 3' end by cutadapt
* Can move UMI to header for downstream analysis (eg with UMI-TOOLS)
* Can handle different length UMIs (eg for mixed demultiplexing of irCLIP and iiCLIP)


# Usage
Create a csv of barcodes. The first column should be the five prime barcodes. The second column (optional) should be the 3' barcodes associated with that 5' barcode. Each 3' barcode should be separated by a semi-colon. There should be no other characters except: "A", "C", "T", "G", "N" (for UMIs), ",", and ";"

```
usage: v6.py [-h] -i INPUTFASTQ -b BARCODES [-m5 [FIVEPRIMEMISMATCHES]]
             [-m3 [THREEPRIMEMISMATCHES]] [-q [PHREDQUALITY]] [-t [THREADS]]
             [-a [ADAPTER]] [-o [OUTPUTPREFIX]] [-sb]

Ultra-fast demultiplexing of fastq files.

required arguments:
  -i INPUTFASTQ, --inputfastq INPUTFASTQ
                        fastq file to be demultiplexed
  -b BARCODES, --barcodes BARCODES
                        barcodes for demultiplexing in csv format

optional arguments:
  -h, --help            show this help message and exit
  -m5 [FIVEPRIMEMISMATCHES], --fiveprimemismatches [FIVEPRIMEMISMATCHES]
                        number of mismatches allowed for 5prime barcode
                        [DEFAULT 1]
  -m3 [THREEPRIMEMISMATCHES], --threeprimemismatches [THREEPRIMEMISMATCHES]
                        number of mismatches allowed for 3prime barcode
                        [DEFAULT 0]
  -q [PHREDQUALITY], --phredquality [PHREDQUALITY]
                        phred quality score for 3prime end trimming
  -t [THREADS], --threads [THREADS]
                        threads [DEFAULT 8]
  -a [ADAPTER], --adapter [ADAPTER]
                        sequencing adapter to trim [DEFAULT Illumina
                        AGATCGGAAGAGCGGTTCAG]
  -o [OUTPUTPREFIX], --outputprefix [OUTPUTPREFIX]
                        prefix for output sequences [DEFAULT demux]
  -sb, --sbatchcompression
                        whether to compress output fastq using SLURM sbatch
```

