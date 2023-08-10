# tRNA_Charge-Seq
This repository contains script files used for analysis of tRNA Charge-Seq data.


Fasta files for spliced tRNA transcripts for human and mouse were obtained from gtrnadb.ucsc.edu and saved as file named hg38-mature-tRNAs.fa (human) and mm10-tRNAs.fa (mouse)

The python scripts appendCCA.py and appendCCAhuman.py was used on the respective human or mouse raw tRNA fasta files to convert RNA to DNA sequences, combine duplicate genes and append the suffix 5'-CCA-3' to the 3' end of each sequence.
The produced files hg38tRNAsCCA.fa and mm10-tRNAsPheCCa.fa were used as the reference genomes for alignment using bowtie2 read aligner.

Following alignment the python script parseSamCharge.py was used to identify the charging status of each read contained within the alinged read files and counts of charged and uncharged reads for each isodecoder were saved as text files.

The text files containing counts of charged and uncharged tRNA were then utilized by condenseCountsMEF.R and condenseCountsNTERT.R to calculate statistics and create graphs.


These scripts were used in the following publications:

Misra, J., Carlson, K., Spandau, D., Wek, R. Multiple mechanisms activate GCN2 eIF2 kinase in response to diverse stress conditions. Nucleic Acids Research. (Under Revision)
