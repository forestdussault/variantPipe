"""

Step 1: Produce BAM files from fastq's aligned against reference (Pseudomonas_aeruginosa_PAO1) with BBMap

Step 2: Tag all sorted bam files with bamaddrg, pass them into a single large bam

Step 3: Sort BAM file with samtools

Step 4: Remove duplicate reads with SAMtools

Step 5: Index data with SAMtools

Step 6: Run tagged bam file through freebayes to produce a VCF

Step 7: Filter VCF to a minimum read depth of 3

NOTES: After building this the first iteration of this pipeline I found a publication released in 2017 that recommends
using BBMap + FreeBayes for variant calling. This is good news.

Here's the link: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0174446
Supplementary File S5 in the publication offers reference commands for running the programs.

Requirements:
    - Python 3
    - FreeBayes
    - vcflib
    - bbmap
    - bamaddrg
    - samtools

"""

import subprocess
import glob
import os


# STEP 1:
def run_bbmap(read1, read2, base_id, ref_genome):
    workdir = '/mnt/nas/bio_requests/9707/pipeline/'
    bam_filename = workdir + base_id + '/' + base_id + '.bam'

    p = subprocess.Popen('bbmap.sh in={0} in2={1} '
                         'out={2} '
                         'covstats={2}_covstats.txt '
                         'covhist={2}_covhist.txt '
                         'basecov={2}_basecov.txt '
                         'bincov={2}_bincov.txt '
                         'sam=1.3 '  # Required for freebayes to work
                         'ref={3} '
                         'ambiguous=toss '  # Recommended by Smith et al. 2017 (see link above for details)
                         'nodisk '
                         'overwrite=true '.format(read1,
                                                  read2,
                                                  bam_filename,
                                                  ref_genome),
                         shell=True,
                         executable="/bin/bash")
    p.wait()

    return bam_filename


# STEP 2:
def run_bamaddrg(*args):
    print('\nRunning bamaddrg...')
    root_cmd = './bamaddrg '
    for arg in args:
        arg_base = os.path.basename(arg)
        root_cmd += ('-b {} -s {} '.format(arg, arg_base))
    root_cmd += ' > /mnt/nas/bio_requests/9707/pipeline/combined_tagged.bam'

    p = subprocess.Popen(root_cmd,
                         shell=True,
                         executable='/bin/bash',
                         cwd='/home/dussaultf/Applications/bamaddrg/bamaddrg')
    p.wait()

    combined_bam = glob.glob('/mnt/nas/bio_requests/9707/pipeline/*_tagged.bam')
    return combined_bam[0]


# STEP 3:
def sort_bamfile(combined_bam):
    print("\nSorting bam file: {}".format(combined_bam))
    # samtools sort -O bam -o SEQUENCE_DATA_SORTED.BAM -T TEMP -@ 4 SEQUENCE_DATA.BAM
    cmd = 'samtools sort -T TEMP -@ 8 -o {0} {1}'.format(combined_bam.replace('.bam', '_sorted.bam'),
                                                         combined_bam)

    p = subprocess.Popen(cmd,
                         shell=True,
                         executable='/bin/bash')
    p.wait()

    sorted_bam = glob.glob('/mnt/nas/bio_requests/9707/pipeline/*_sorted.bam')
    return sorted_bam[0]


# STEP 4:
def dedupe_bam(sorted_bam):
    print("\nDeduping bam file: {}".format(sorted_bam))
    # samtools rmdup -s SEQUENCE_DATA_SORTED.BAM SEQUENCE_DATA_DEDUP.BAM
    cmd = 'samtools rmdup -s {0} {1}'.format(sorted_bam, sorted_bam.replace('_sorted.bam', '_deduped.bam'))

    p = subprocess.Popen(cmd,
                         shell=True,
                         executable='/bin/bash')
    p.wait()

    deduped_bam = glob.glob('/mnt/nas/bio_requests/9707/pipeline/*_deduped.bam')
    return deduped_bam[0]


# STEP 5:
def index_bam(deduped_bam):
    print("\nIndexing bam file: {}".format(deduped_bam))

    cmd = 'samtools index {0}'.format(deduped_bam)
    p = subprocess.Popen(cmd,
                         shell=True,
                         executable='/bin/bash')
    p.wait()


# STEP 6:
def run_freebayes(*args):
    # TODO: use freebayes-parallel instead. Must be run from the freebayes installation directory.
    print("\nRunning freebayes...")
    bams = ''
    for arg in args:
        bams += arg

    bams += ' > /mnt/nas/bio_requests/9707/pipeline/var.vcf'

    # Set ploidy to 1 with -p 1. If this isn't specified the program will assume diploid.
    cmd = 'freebayes -p 1 -f {0} {1}'.format(
        '/mnt/nas/bio_requests/9707/ref/Pseudomonas_simple_name.fna', bams)
    p = subprocess.Popen(cmd,
                         shell=True,
                         executable='/bin/bash')
    p.wait()

    filtered_vcf = glob.glob('/mnt/nas/bio_requests/9707/pipeline/*.vcf')
    return filtered_vcf[0]


# STEP 7:
def run_vcffilter(vcf_file):
    """
    See some vcffilter examples at https://www.biostars.org/p/51439/
    """
    print("\nRunning vcffilter on {}".format(vcf_file))

    cmd = 'vcffilter -f "DP > 2 & QUAL > 20" {0} > {1}'.format(vcf_file, vcf_file.replace('var', 'var_filtered'))

    p = subprocess.Popen(cmd,
                         shell=True,
                         executable='/bin/bash')
    p.wait()


########################
# RUNNING THE PIPELINE #
########################
def main():
    # 1. Producing sorted BAM files for 0766, 0767, 0768
    ref_genome_fasta = '/mnt/nas/bio_requests/9707/ref/Pseudomonas_simple_name.fna'

    bam1 = run_bbmap(read1='/mnt/nas/bio_requests/9707/fastq/2017-SEQ-0766/2017-SEQ-0766_R1_trimmed.fastq.bz2',
                     read2='/mnt/nas/bio_requests/9707/fastq/2017-SEQ-0766/2017-SEQ-0766_R2_trimmed.fastq.bz2',
                     base_id='2017-SEQ-0766',
                     ref_genome=ref_genome_fasta)

    bam2 = run_bbmap(read1='/mnt/nas/bio_requests/9707/fastq/2017-SEQ-0767/2017-SEQ-0767_R1_trimmed.fastq.bz2',
                     read2='/mnt/nas/bio_requests/9707/fastq/2017-SEQ-0767/2017-SEQ-0767_R2_trimmed.fastq.bz2',
                     base_id='2017-SEQ-0767',
                     ref_genome=ref_genome_fasta)

    bam3 = run_bbmap(read1='/mnt/nas/bio_requests/9707/fastq/2017-SEQ-0768/2017-SEQ-0768_R1_trimmed.fastq.bz2',
                     read2='/mnt/nas/bio_requests/9707/fastq/2017-SEQ-0768/2017-SEQ-0768_R2_trimmed.fastq.bz2',
                     base_id='2017-SEQ-0768',
                     ref_genome=ref_genome_fasta)

    # 2. Tagging/merging all BAM files with bamaddrg
    combined_bam = run_bamaddrg(bam1, bam2, bam3)

    # 3. Sort BAM
    sorted_bamfile = sort_bamfile(combined_bam)

    # 4. Dedupe BAM
    deduped_bam = dedupe_bam(sorted_bamfile)

    # 5. Index BAM
    index_bam(deduped_bam)

    # 6. Run combined/tagged bam files through freebayes to determine SNVs/indels
    filtered_vcf = run_freebayes(deduped_bam)

    # 7. Run vcffilter with minimum read depth of 3
    run_vcffilter(filtered_vcf)


if __name__ == '__main__':
    main()
