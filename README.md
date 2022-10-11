# tRNA_Charge-Seq
This repository contains script files used for analysis of tRNA Charge-Seq data.


Fasta files for spliced tRNA transcripts were obtained from gtrnadb.ucsc.edu and saved as file named hg38-mature-tRNAs.fa

The python script appendCCAhuman.py was used to convert RNA to DNA sequences, combine duplicate genes and append the suffix 5'-CCA-3' to the 3' end of each sequence.
The produced file hg38tRNAsCCA.fa was used as the reference genome for alignment using bowtie read aligner.

Following alignment the python script parseSamCharge.py was used to identify the charging status of each read contained within the alinged read files and counts of charged and uncharged reads for each isodecoder were saved as text files.

The text files containing counts of charged and uncharged tRNA were then utilized by condenseCountsRicardo.R to calculate statistics and create graphs.


These scripts were used in the following publications:
Cordova, R. A., Misra, J., Amin, P. H., Klunk, A. J., Damayanti, N. P., Carlson, K. R., Elmendorf, A. J., Kim, H.-G., Mirek, E. T., Elzey, B. D., Miller, M. J., Dong, X. C., Cheng, L., Anthony, T. G., Pili, R., Wek, R. C., and Staschke, K. A. (2022) GCN2 eIF2 kinase promotes prostate cancer by maintaining amino acid homeostasis. eLife 11, e81083. [https://doi.org/10.7554/eLife.81083](https://doi.org/10.7554/eLife.81083)
