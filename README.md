## variantPipe
A pipeline for calling variants with freebayes.

[Reference publication for using BBMap + freebayes as a preferred variant calling method](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0174446)

#### Pipeline overview
1. Produce bam files from paired fastq files aligned against a reference fasta file with BBMap or bowtie2
2. Tag all sorted bam files with bamaddrg, pass them into a single large bam
3. Sort bam file with samtools
4. Remove duplicate reads with SAMtools
5. Index data with SAMtools
6. Run tagged bam file through freebayes to produce a VCF
7. Filter VCF to a minimum read depth of 3, QUAL > 20

#### Requirements
*The following must be available in your PATH*
- Python >= 3.5
- freebayes
- vcflib
- bbmap
- bamaddrg
- samtools
- bowtie2

#### Usage
This is a command line tool that requires at least three inputs.
1. Set of paths to folders that contain a pair of .fastq.bz2 files in each (_R1, _R2). Note that samples will be labelled by the parent folder name.
2. Path to the reference genome `-r, --refgenome`
3. Path to your desired output directory `-o, --outputdir`

*Example command:*

`python3 variantpipe.py /mnt/nas/bio_requests/9707/fastq/2017-SEQ-0766 /mnt/nas/bio_requests/9707/fastq/2017-SEQ-0767 /mnt/nas/bio_requests/9707/fastq/2017-SEQ-0768
-r /mnt/nas/bio_requests/9707/ref/Pseudomonas_simple_name.fna
-o /mnt/nas/bio_requests/9707/CLI_testing`

BBMap can be swapped out for bowtie2 as the aligner with the flag `-bt, --bowtie2`

