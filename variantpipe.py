"""
To do:
- Cleanup intermediary .bam files
- Actually implement support for bowtie2
- Allow support for more than just .fastq.bz2 files
- Dockerize package
"""

import subprocess
import shutil
import time
import argparse
import glob
import os


class VariantPipe(object):
    def sample_filer(self):
        sample_dict = {}
        for item in self.sample_targets:
            base_id = os.path.basename(item)
            fastq_list = glob.glob(item + '/*.fastq.bz2')

            if len(fastq_list) > 2:
                print('Invalid number of .fastq.bz2 files detected. Max. 2 in folder.')
                quit()

            sample_dict[base_id] = fastq_list
        return sample_dict

    @staticmethod
    def exec_command(cmd):
        p = subprocess.Popen(cmd,
                             shell=True,
                             executable="/bin/bash")
        p.wait()

    @staticmethod
    def delete_file(file):
        try:
            os.remove(file)
        except OSError:
            print('\nCould not remove {}'.format(file))

    def run_qualimap(self, bam):
        # Must sort prior to running qualimap
        sorted_bam = self.sort_bamfile(bam)
        cmd = 'qualimap bamqc -bam {} -nw 400 -hm 3'.format(sorted_bam)
        self.exec_command(cmd)

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
            except OSError:
                shutil.rmtree(workdir + key)
                os.mkdir(workdir + key)

            print('\nFile 1: {0}\nFile 2: {1}\nOutput dir: {2}\n'.format(sample_dict[key][0],
                                                                         sample_dict[key][1],
                                                                         os.path.dirname(mapped_bam)))

            cmd = 'bbmap.sh ' \
                  'in={0} ' \
                  'in2={1} ' \
                  'outm={2} ' \
                  'outu={4} ' \
                  'covstats={2}_covstats.txt ' \
                  'covhist={2}_covhist.txt ' \
                  'basecov={2}_basecov.txt ' \
                  'bincov={2}_bincov.txt ' \
                  'sam=1.3 ' \
                  'ref={3} ' \
                  'ambiguous=toss ' \
                  'nodisk ' \
                  'slow ' \
                  'overwrite=true'.format(sample_dict[key][0],
                                          sample_dict[key][1],
                                          mapped_bam,
                                          self.refgenome,
                                          unmapped_reads),
            self.exec_command(cmd)

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

        self.exec_command(root_cmd)

        combined_bam = glob.glob('{}*_tagged.bam'.format(self.outputdir))

        # Cleanup
        for bam in mapped_bam_list:
            self.delete_file(bam)

        return combined_bam[0]

    # STEP 3:
    def sort_bamfile(self, bam):
        print("\nSorting bam file: {}".format(bam))

        sorted_bam = bam.replace('.bam', '_sorted.bam')
        cmd = 'samtools sort -T TEMP -@ 8 -o {0} {1}'.format(sorted_bam, bam)

        self.exec_command(cmd)

        return sorted_bam

    # STEP 4:
    def dedupe_bam(self, sorted_bam):
        print("\nDeduping bam file: {}".format(sorted_bam))
        deduped_bam = sorted_bam.replace('_sorted.bam', '_deduped.bam')
        cmd = 'samtools rmdup -s {0} {1}'.format(sorted_bam, deduped_bam)

        self.exec_command(cmd)

        # Cleanup
        self.delete_file(sorted_bam)

        return deduped_bam

    # STEP 5:
    def index_bam(self, deduped_bam):
        print("\nIndexing bam file: {}".format(deduped_bam))

        cmd = 'samtools index {0}'.format(deduped_bam)
        self.exec_command(cmd)

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

        bams += ' > {}variants.vcf'.format(self.outputdir)

        # Set ploidy to 1 with -p 1. If this isn't specified the program will assume diploid.
        cmd = 'freebayes -p 1 -f {0} {1}'.format(self.refgenome, bams)
        self.exec_command(cmd)

        filtered_vcf = glob.glob('{}*.vcf'.format(self.outputdir))
        return filtered_vcf[0]

    # STEP 7:
    def run_vcffilter(self, vcf_file):
        """
        See some vcffilter examples at https://www.biostars.org/p/51439/
        """
        print("\nRunning vcffilter on {}".format(vcf_file))

        # vcffilter set to minimum read depth of 3 and quality > 20
        cmd = 'vcffilter -f "DP > 2 & QUAL > 20" {0} > {1}'.format(vcf_file, vcf_file.replace('variants',
                                                                                              'variants_filtered'))

        self.exec_command(cmd)

    def __init__(self, args):
        print('\033[92m' + '\033[1m' + '\nVariantPipe' + '\033[0m')

        # Get args
        self.args = args

        # Required input
        self.sample_targets = args.sample_targets
        self.refgenome = args.refgenome
        self.outputdir = args.outputdir
        self.num_samples = len(self.sample_targets)

        # Optional input
        self.bowtie2flag = args.bowtie2
        self.qualimapflag = args.qualimap

        # Validate output dir path
        if self.outputdir.endswith('/'):
            pass
        else:
            self.outputdir += '/'

        # Create output dir if it doesn't exist. Overwrite otherwise.
        try:
            os.mkdir(self.outputdir)
        except OSError:
            shutil.rmtree(self.outputdir)
            os.mkdir(self.outputdir)

        # Get samples
        sample_dict = self.sample_filer()

        mapped_bam_list = []

        # Run BBMap on given samples
        if self.bowtie2flag:
            # Run bowtie2
            pass
        else:
            mapped_bam_list = self.run_bbmap(sample_dict=sample_dict)

        # Optionally run qualimap
        if self.qualimapflag:
            for bam in mapped_bam_list:
                self.run_qualimap(bam)

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
    parser.add_argument('-qm', '--qualimap',
                        default=False,
                        action='store_true',
                        help='Set this flag to run qualimap on the individual output BAM files')

    arguments = parser.parse_args()

    x = VariantPipe(arguments)

    end = time.time()
    m, s = divmod(end - start, 60)
    h, m = divmod(m, 60)

    print('\033[92m' + '\033[1m' + '\nFinished VariantPipe functions in %d:%02d:%02d ' % (h, m, s) + '\033[0m')
