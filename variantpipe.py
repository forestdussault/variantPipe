import subprocess
import shutil
import time
import argparse
import glob
import os


class VariantCaller(object):
    def sample_filer(self):
        sample_dict = {}
        for item in self.sample_targets:
            base_id = os.path.basename(item)
            fastq_list = glob.glob(item+'/*.fastq.bz2')

            if len(fastq_list) > 2:
                print('Invalid number of .fastq.bz2 files detected. Max. 2 in folder.')
                quit()

            sample_dict[base_id] = fastq_list
        return sample_dict

    # STEP 1A:
    def run_bbmap(self, sample_dict):
        workdir = self.outputdir
        mapped_bam_list = []
        unmapped_read_list = []

        for key in sample_dict:
            print('\nRunning BBMap on {}...'.format(key))
            mapped_bam = workdir + key + '/' + key + '.bam'
            unmapped_reads = workdir + key + '/' + key + '_unmapped.fq'

            mapped_bam_list.append(mapped_bam)
            unmapped_read_list.append(unmapped_reads)

            # Make result folders for BBmap
            try:
                os.mkdir(workdir + key)
            except:
                shutil.rmtree(workdir + key)
                os.mkdir(workdir + key)

            print('\nFile 1: {0}\nFile 2: {1}\nOutput dir: {2}\n'.format(sample_dict[key][0],
                                                                         sample_dict[key][1],
                                                                         os.path.dirname(mapped_bam)))

            p = subprocess.Popen('bbmap.sh in={0} in2={1} '
                                 'outm={2} '
                                 'outu={4} '
                                 'covstats={2}_covstats.txt '
                                 'covhist={2}_covhist.txt '
                                 'basecov={2}_basecov.txt '
                                 'bincov={2}_bincov.txt '
                                 'sam=1.3 '  # Required for freebayes to work
                                 'ref={3} '

                                 # # High sensitivity
                                 # 'slow '
                                 # 'k=11 '
                                 # 'minratio=0.1 '

                                 'ambiguous=toss '  # Recommended by Smith et al. 2017 (see link above for details)
                                 'nodisk '
                                 'overwrite=true '.format(sample_dict[key][0],
                                                          sample_dict[key][1],
                                                          mapped_bam,
                                                          self.refgenome,
                                                          unmapped_reads),
                                 shell=True,
                                 executable="/bin/bash")
            p.wait()

        return mapped_bam_list

    # Alternate STEP 1B:
    # TODO: Create another function to create a bowtie2 index with bowtie2-build from the ref genome
    # def run_bowtie2(self, sample_dict):
    #     mapped_bam = self.outputdir + base_id + '/' + base_id + '.bam'
    #
    #     try:
    #         os.mkdir(workdir + base_id)
    #     except OSError:
    #         shutil.rmtree(workdir + base_id)
    #         os.mkdir(workdir + base_id)
    #
    #     p = subprocess.Popen('bowtie2  \
    #                          -p 8 \
    #                          -x {3}  \
    #                          -1 {0}  \
    #                          -2 {1}  \
    #                          --sensitive-local \
    #                          | samtools view -bS > {2}'.format(read1,
    #                                                            read2,
    #                                                            mapped_bam,
    #                                                            bt2_index),
    #                          shell=True,
    #                          executable="/bin/bash")
    #     p.wait()
    #
    #     return mapped_bam

    # STEP 2:
    def run_bamaddrg(self, mapped_bam_list):
        print('\nRunning bamaddrg on the following:')
        for bam in mapped_bam_list:
            print(os.path.basename(bam))

        root_cmd = 'bamaddrg '
        for bam in mapped_bam_list:
            bam_base = os.path.basename(bam)
            root_cmd += ('-b {} -s {} '.format(bam, bam_base))
        root_cmd += ' > {}combined_tagged.bam'.format(self.outputdir)

        p = subprocess.Popen(root_cmd,
                             shell=True,
                             executable='/bin/bash')
        p.wait()

        combined_bam = glob.glob('{}*_tagged.bam'.format(self.outputdir))
        return combined_bam[0]

    # STEP 3:
    def sort_bamfile(self, combined_bam):
        print("\nSorting bam file: {}".format(combined_bam))
        cmd = 'samtools sort -T TEMP -@ 8 -o {0} {1}'.format(combined_bam.replace('.bam', '_sorted.bam'),
                                                             combined_bam)

        p = subprocess.Popen(cmd,
                             shell=True,
                             executable='/bin/bash')
        p.wait()

        sorted_bam = glob.glob('{}*_sorted.bam'.format(self.outputdir))
        return sorted_bam[0]

    # STEP 4:
    def dedupe_bam(self, sorted_bam):
        print("\nDeduping bam file: {}".format(sorted_bam))
        # samtools rmdup -s SEQUENCE_DATA_SORTED.BAM SEQUENCE_DATA_DEDUP.BAM
        cmd = 'samtools rmdup -s {0} {1}'.format(sorted_bam, sorted_bam.replace('_sorted.bam', '_deduped.bam'))

        p = subprocess.Popen(cmd,
                             shell=True,
                             executable='/bin/bash')
        p.wait()

        deduped_bam = glob.glob('{}*_deduped.bam'.format(self.outputdir))
        return deduped_bam[0]

    # STEP 5:
    @staticmethod
    def index_bam(deduped_bam):
        print("\nIndexing bam file: {}".format(deduped_bam))

        cmd = 'samtools index {0}'.format(deduped_bam)
        p = subprocess.Popen(cmd,
                             shell=True,
                             executable='/bin/bash')
        p.wait()

    # STEP 6:
    def run_freebayes(self, *args):
        """
        Can take several BAM files, though currently I'm combining all BAM files with bamaddrg prior to this step.
        """
        # TODO: use freebayes-parallel instead
        print("\nRunning freebayes on the following:")
        for arg in args:
            print(os.path.basename(arg))

        bams = ''
        for arg in args:
            bams += arg

        bams += ' > {}var.vcf'.format(self.outputdir)

        # Set ploidy to 1 with -p 1. If this isn't specified the program will assume diploid.
        cmd = 'freebayes -p 1 -f {0} {1}'.format(self.refgenome, bams)
        p = subprocess.Popen(cmd,
                             shell=True,
                             executable='/bin/bash')
        p.wait()

        filtered_vcf = glob.glob('{}*.vcf'.format(self.outputdir))
        return filtered_vcf[0]

    # STEP 7:
    @staticmethod
    def run_vcffilter(vcf_file):
        """
        See some vcffilter examples at https://www.biostars.org/p/51439/
        """
        print("\nRunning vcffilter on {}".format(vcf_file))

        # vcffilter set to minimum read depth of 3 and quality > 20
        cmd = 'vcffilter -f "DP > 2 & QUAL > 20" {0} > {1}'.format(vcf_file, vcf_file.replace('var', 'var_filtered'))

        p = subprocess.Popen(cmd,
                             shell=True,
                             executable='/bin/bash')
        p.wait()

    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nVariantCaller' + '\033[0m')

        # Get args
        self.args = args

        # Required input
        self.sample_targets = args.sample_targets
        self.refgenome = args.refgenome
        self.outputdir = args.outputdir

        # Optional input
        self.bowtie2flag = args.bowtie2

        # Validate output dir path
        if self.outputdir.endswith('/'):
            pass
        else:
            self.outputdir += '/'

        # Get samples to work withb
        sample_dict = self.sample_filer()

        # Run BBMap on given samples
        mapped_bam_list = self.run_bbmap(sample_dict=sample_dict)

        # Pass mapped_bam_list to bamaddrg
        combined_bam = self.run_bamaddrg(mapped_bam_list)

        # 3. Sort BAM
        sorted_bamfile = self.sort_bamfile(combined_bam)

        # 4. Dedupe BAM
        deduped_bam = self.dedupe_bam(sorted_bamfile)

        # 5. Index BAM
        self.index_bam(deduped_bam)

        # 6. Run combined/tagged bam files through freebayes to determine SNVs/indels
        filtered_vcf = self.run_freebayes(deduped_bam)

        # 7. Run vcffilter with minimum read depth of 3, Q > 20
        self.run_vcffilter(filtered_vcf)

if __name__ == '__main__':
    start = time.time()
    parser = argparse.ArgumentParser()

    parser.add_argument('sample_targets',
                        help='Path(s) to the folders that contain paired .fastq.bz2 files.',
                        nargs='*')
    parser.add_argument('-r', '--refgenome',
                        help='Specify the path to the reference genome.',
                        required=True)
    parser.add_argument('-o', '--outputdir',
                        help='Specify the path to your desired output directory.',
                        required=True)
    # TODO: Make the bowtie2 flag actually run bowtie2 instead of BBmap
    parser.add_argument('-bt', '--bowtie2',
                        default=False,
                        action='store_true',
                        help='Uses bowtie2 as the aligner instead of the default BBmap')

    arguments = parser.parse_args()

    x = VariantCaller(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    print('\033[92m' + '\033[1m' + '\nFinished VariantCaller functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')
