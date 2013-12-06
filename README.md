HMMSplicerBEDToSAMParser
========================

This is a perl script to parse the non-colpased BED file generated from HMMSplicer into SAM alignment. It will first extract the successfully aligned (unspliced) sequences from the initial SAM file generated from HMMSplicer (have to turn on the -S option in bowtie command of the HMMSPLicer python script), and it reads the BED file for spliced reads, and takes the sequence from the SAM file (unaligned ones). User have to specify paired-end or single-end. In case of paired-end, the read name has to be ended with 1 and 2. The aim of this script is to generate a SAM file that is comparable to that of Tophat so that it can be loaded to Cufflinks.
