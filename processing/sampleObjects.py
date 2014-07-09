import sys
sys.path.append('lib/python2.7/site-packages/')
import os
import logging
import subprocess
from glob import glob
from tools import *
from processing.mappedReadsTools import *
from multiprocessing import cpu_count
from math import floor

def getMaxCpu(max_cpu):
    '''
    # set the CPU number for threading
    '''
    cpu_num = cpu_count() - 1
    if cpu_num > max_cpu:
        cpu_num = max_cpu

    return max_cpu

class SampleContainer(dict):
    def __init__(self, sample_dir):
        self.sample_dir = sample_dir
        self._loadSampleMetadata()
        self._initLogger()

    def _loadSampleMetadata(self):
        metadata_dict = self.parseSampleMetaData()
        for key in metadata_dict:
            setattr(self, key, metadata_dict[key])

    def parseSampleMetaData(self,):
        metadata_fn = self.sample_dir + '/SampleID'
        metadata_fh = open(metadata_fn)

        sample_metadata_dict = {}

        for line in metadata_fh:
            if ';' not in line:
                continue
            line = line.split(';')
            try:
                sample_metadata_dict[line[0]] = int(line[1].strip())
            except ValueError:
                sample_metadata_dict[line[0]] = str(line[1].strip())

        return sample_metadata_dict

    def _initLogger(self):
        log_name = '%s processing: %s' % (self.experiment_type, self.sample_name)
        self.logger = openLog(log_name, 'log/' + log_name + '.log', console_log=True)

    def _reportLoadedParams():
        self.logger.info("Processing sample with files in %s..." % self.sample_dir)
        self.logger.info("\tread type: %s" % self.read_type)
        self.logger.info("\tgenome version: %s" % self.genome_version)
        self.logger.info("\tsingle-end adapter sequence: %s" % self.adapter_se)
        self.logger.info("\tpaired-end adapter sequence: %s" % self.adapter_pe)
        self.logger.info("\tinsert size: %s" % self.insert_size)

    def gunzipFastqs(self, fastq_dir_prefix='FASTQ', single_substr='R1', paired_substr='R2', combined_fastq_fn_suffix='.all.fastq'):
        '''
        Unzips and concatenates the FASTQ files from an Illumina experiment
        for downstream analysis.

        Returns a dict of the names of the single and paired end file names for other functions to use.

        The concatenated file is stored in a subdirectory of the
        sample directory.

        The substr identifies the read file as single or paired end.
        '''
        def _getFastqgzFilenameList():
            fastqgz_fn_list = glob(self.sample_dir + '/' + '*fastq.gz')

            if fastqgz_fn_list == []:
                raise ValueError("Error finding the FASTQ files in the sample directory %s... exiting." % (self.sample_dir))
            else:
                self.logger.info("Found gzipped fastq files in %s..." % (self.sample_dir))
                return fastqgz_fn_list

        def _checkForCombinedFastq():
            '''
            Check for existing combine fastq files.

            If they exist, return them inthe combined_fastq_fn_dict,
            if not,
            '''
            self.logger.info("Checking for previously unzipped fastq files...")

            combined_fastq_fn_dict = {}
            self.combined_fastq_fn_dict = combined_fastq_fn_dict

            single_combined_fastq_fn = self.fastq_dir + '/' + self.sample_name + '_' + single_substr + combined_fastq_fn_suffix
            paired_combined_fastq_fn = self.fastq_dir + '/' + self.sample_name + '_' + paired_substr + combined_fastq_fn_suffix

            if os.path.isfile(single_combined_fastq_fn):
                self.logger.info("...found the single-end combined fastq file.")
                self.combined_fastq_fn_dict['single'] = single_combined_fastq_fn
            else:
                return False

            if os.path.isfile(paired_combined_fastq_fn):
                self.logger.info("...found the paired-end combined fastq file.")
                self.combined_fastq_fn_dict['paired'] = paired_combined_fastq_fn

            return True

        def _combineFastqgzFiles():
            '''
            Uncompresses and concatenates gzipped fastq files.
            '''
            # get the list of fastq.gz files. This catches ValueError and exits if files don't exist.
            try:
                fastqgz_fn_list = _getFastqgzFilenameList()
                self.logger.info("Found fastq.gz files to combine in %s..." % (self.sample_dir))
            except ValueError as error:
                self.logger.error(error)
                sys.exit()

            self.logger.info("...creating a new FASTQ directory for this analysis at %s." % (self.fastq_dir))
            try:
                os.mkdir(self.fastq_dir)
            except:
                self.logger.warning("...oh, wait, it was already there, but without combined fastq files.")


            single_fastqgz_fn_list, paired_fastqgz_fn_list = [], []
            for fastqgz_fn in fastqgz_fn_list:
                if single_substr in fastqgz_fn:
                    single_fastqgz_fn_list.append(fastqgz_fn)
                elif paired_substr in fastqgz_fn:
                    paired_fastqgz_fn_list.append(fastqgz_fn)

            # this sort ensures that paired reads stay in the correct order, but
            # it also assumes that there are order-identifying characters in the
            # read filenames, e.g. the R001, R002 in :
            # LL1_index1_ATCACG_L007_R2_001.fastq.gz, LL1_index1_ATCACG_L007_R2_002.fastq.gz
            single_fastqgz_fn_list.sort()
            paired_fastqgz_fn_list.sort()

            for fastqgz_fn in fastqgz_fn_list:
                if (fastqgz_fn not in single_fastqgz_fn_list) and (fastqgz_fn not in paired_fastqgz_fn_list):
                    self.logger.warning("Skipping %s because it doesn't contain a single or paired end FASTQ substring!" % (fastqgz_fn))

            combined_fastq_fn_dict = {}

            if len(single_fastqgz_fn_list) > 0:
                single_combined_fastq_fn = self.fastq_dir + '/' + self.sample_name + '_' + single_substr + combined_fastq_fn_suffix
                os.system('zcat ' + ' '.join(single_fastqgz_fn_list) + '> ' + single_combined_fastq_fn)
                combined_fastq_fn_dict['single'] = single_combined_fastq_fn
            if len(paired_fastqgz_fn_list) > 0:
                paired_combined_fastq_fn = self.fastq_dir + '/' + self.sample_name + '_' + paired_substr + combined_fastq_fn_suffix
                os.system('zcat ' + ' '.join(paired_fastqgz_fn_list) + '> ' + paired_combined_fastq_fn)
                combined_fastq_fn_dict['paired'] = paired_combined_fastq_fn

            self.combined_fastq_fn_dict = combined_fastq_fn_dict

        #####
        ## function body beings here
        #####
        # set the output directory
        self.fastq_dir = self.sample_dir + '/' + '_'.join([fastq_dir_prefix, self.read_type, self.genome_version])

        # check for existing unzipped files. If the they exist, the
        # combined_fastq_fn_dict is updated and we return
        if _checkForCombinedFastq():
            return

        _combineFastqgzFiles()

    def trimFastqFiles(self, trimmed_fastq_fn_suffix='.adaptercleaned.fastq'):
        '''
        Trims given adapter sequences out of the fastq files.
        '''
        def _checkForTrimmedFastqFn():
            '''
            '''
            if os.path.isfile(trimmed_fastq_fn):
                self.logger.warning("It looks like the adapter sequences were already trimmed from fastq file %s !" % (trimmed_fastq_fn))
                return True

            return False

        def trimFastqFile(combined_fastq_fn_key):
            '''
            '''
            def _finishWhileLoop():
                trim_counts_fh.write('Total Reads = %s\n' % (read_count))
                trim_counts_fh.write('Trimmed Reads = %s\n' % (trim_count))
                percent = float(trim_count) / read_count * 100
                trim_counts_fh.write('Percent trimmed = %s\n' % (percent))

                self.logger.info("\tFrom a total of %s, %s were trimmed (%s percent)." % (read_count, trim_count, percent))

                combined_fastq_fh.close()
                trimmed_fastq_fh.close()
                trim_counts_fh.close()

            def _removeAdapter(read, adapter_seq):
                '''
                Removes an adapter_sequence from a read by taking a 3-nt match and extending it.
                '''
                idx = read.rfind(adapter_seq)
                if idx != -1:
                    self.logger.debug('Found the adapter at pos %s in read %s' % (idx, read))
                    return read[0:idx]
                else:
                    idx = read.rfind(adapter_seq[0:3])
                    if idx != -1:
                        pad = len(read) - idx
                        idx = read.rfind(adapter_seq[0:pad])
                        if idx != -1:
                            self.logger.debug('Found the incomplete adapter at pos %s in read %s' % (idx, read))
                            return read[0:idx]
                        else:
                            self.logger.debug('Not even incomplete adapter found in read %s' % (read))
                            return read
                    else:
                        self.logger.debug('No adapter found in read %s' % (read))
                        return read
            #####
            ## trimming function body starts here
            #####
            combined_fastq_fh = open(self.combined_fastq_fn_dict[combined_fastq_fn_key])
            trimmed_fastq_fn = self.trimmed_fastq_fn_dict[combined_fastq_fn_key]
            trimmed_fastq_fh =  open(trimmed_fastq_fn, 'w')
            trim_counts_fn = trimmed_fastq_fn.split('.')[0] + '.trimcounts'
            trim_counts_fh = open(trim_counts_fn, 'w')

            if combined_fastq_fn_key == 'single':
                adapter_seq = self.adapter_se
            elif combined_fastq_fn_key == 'paired':
                adapter_seq = self.adapter_pe

            self.logger.info("\nRemoving adapter sequences from %s..." % (trimmed_fastq_fn))

            # we want to count and report the number of reads and the number trimmed
            read_count = 0
            trim_count = 0
            while (True):
                for i in range(4):
                    read_name = combined_fastq_fh.readline().strip()
                    if read_name == '':
                        _finishWhileLoop()
                        return
                    read = combined_fastq_fh.readline().strip()
                    untrimmed_read_len = len(read)
                    read = _removeAdapter(read, adapter_seq)
                    if len(read) != untrimmed_read_len:
                        trim_count += 1
                    # if the read is entirely composed of adapter, then leave 'Ns' to
                    # preserve paired-end read order.  At least Two Ns are required for bowtie to handle them.
                    if read == '':
                        read = 'N' * untrimmed_read_len
                    comment = combined_fastq_fh.readline().strip()
                    qscore = combined_fastq_fh.readline().strip()
                    trimmed_fastq_fh.write(read_name + '\n')
                    trimmed_fastq_fh.write(read + '\n')
                    trimmed_fastq_fh.write(comment + '\n')
                    trimmed_fastq_fh.write(qscore[0:len(read)] + '\n')
                    read_count += 1

                if read_count % 1e6 == 0:
                    self.logger.info("%s\t%s" % (read_count, read) )
        #####
        ## main function starts here
        #####
        self.trimmed_fastq_fn_dict = {}

        for combined_fastq_fn_key in self.combined_fastq_fn_dict:
            trimmed_fastq_fn = self.combined_fastq_fn_dict[combined_fastq_fn_key].split('.')[0] + trimmed_fastq_fn_suffix
            if _checkForTrimmedFastqFn():
                self.trimmed_fastq_fn_dict[combined_fastq_fn_key] = trimmed_fastq_fn
                continue
            self.trimmed_fastq_fn_dict[combined_fastq_fn_key] = trimmed_fastq_fn
            trimFastqFile(combined_fastq_fn_key)

    def mapBWA(self, max_cpu=1, min_mapq_score=37, remove_duplicate_reads_flag=False):
        '''
        Wrapper fucntion for mapping reads with BWA (as opposed to something else).
        '''
        def _selectFastqType():
            '''
            The Snyder RNA-Seq reads don't get trimmed, and so we need
            to cheat that step without breaking the rest of the code.
            '''
            # used untrimmed reads if rna-seq
            if self.experiment_type == 'rna-seq':
                self.trimmed_fastq_fn_dict = self.combined_fastq_fn_dict
            elif self.experiment_type == 'mnase-seq':
                self.trimmed_fastq_fn_dict = self.trimmed_fastq_fn_dict
            elif self.experiment_type == 'chip-seq':
                self.trimmed_fastq_fn_dict = self.trimmed_fastq_fn_dict
            else:
                self.logger.error("experiment_type invalid: %s" % self.experiment_type )

        def _checkForExistingBAM():
            test_bam_fn = os.path.dirname(self.trimmed_fastq_fn_dict[self.read_type]) + '/' + 'bwa' + '/' +\
                os.path.basename(self.trimmed_fastq_fn_dict[self.read_type]).split('.fastq')[0] + '.bam'

            if os.path.isfile(test_bam_fn):
                self.logger.warning("BAM file already exists!  %s" % (test_bam_fn) )
                self.bam_fn = test_bam_fn
                return True
            else:
                return False

        _selectFastqType()
        if _checkForExistingBAM():
            return

        ## if no BAM file found, we make it
        # first make the index
        self.indexReadsBWA(max_cpu=max_cpu)
        # then align to it
        self.alignBWA(max_cpu=max_cpu)

        # then turn the SAM file into a BAM file
        self.bam_fn = makeBAMfromSAM(self.sam_fn, min_mapq_score, remove_duplicate_reads_flag=remove_duplicate_reads_flag)

    def indexReadsBWA(self, single_substr='R1', paired_substr='R2', max_cpu=1):
        '''
        Index the reads with BWA.

        the file name of the BWA index directory is hardcoded in bwa_idx_dirname
        '''
        if self.genome_version == 'sacCer3':
            bwa_idx_prefix = 'reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BWAIndex/genome.fa'
        else:
            self.logger.error("Invalid genome_version name for %s" % (self.genome_version))
            sys.exit()

        sai_fn_dict = {}
        for trimmed_fastq_fn_key in self.trimmed_fastq_fn_dict:
            # make a directory to store things in
            path = os.path.dirname(self.trimmed_fastq_fn_dict[trimmed_fastq_fn_key])
            bwa_out_dir=path + '/bwa'
            if not os.path.isdir(bwa_out_dir):
                os.mkdir(bwa_out_dir)

            sai_fn = bwa_out_dir + '/' + os.path.basename(self.trimmed_fastq_fn_dict[trimmed_fastq_fn_key]).split('.fastq')[0] + '.sai'
            if sai_fn.find(single_substr) != -1:
                sai_fn_dict['single'] = sai_fn
            elif sai_fn.find(paired_substr) != -1:
                sai_fn_dict['paired'] = sai_fn
            sai_out_fn = bwa_out_dir + '/' + os.path.basename(self.trimmed_fastq_fn_dict[trimmed_fastq_fn_key]).split('.fastq')[0] + '.saiout'
            sai_out_fh = open(sai_out_fn, 'w')
            self.logger.info("Making BWA idx file %s." % (sai_fn))

            max_cpu = getMaxCpu(max_cpu)

            bwa_bin = 'bin/bwa'
            bwa_idx_fn = bwa_idx_prefix + 'bwaIndex'
            bwa_args = [bwa_bin, 'aln', '-t', str(max_cpu), bwa_idx_prefix, self.trimmed_fastq_fn_dict[trimmed_fastq_fn_key], '-f', sai_fn]
            proc = subprocess.Popen(bwa_args, stdout=sai_out_fh, stderr=sai_out_fh)
            proc.communicate()
            sai_out_fh.close()

        self.sai_fn_dict = sai_fn_dict


    def alignBWA(self, single_substr='R1', paired_substr='R2', max_cpu=1):
        '''
        Align reads with BWA.
        '''
        if self.genome_version == 'sacCer3':
            bwa_idx_prefix = 'reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BWAIndex/genome.fa'
        else:
            self.logger.error("Invalid genome_version name for %s" % (self.genome_version))
            sys.exit()

        max_cpu = getMaxCpu(max_cpu)

        if self.read_type == 'single':
            self.logger.info('Aligning single-end reads...')
            bwa_out_dir = os.path.dirname(self.sai_fn_dict['single'])
            sam_fn = self.sai_fn_dict[self.read_type].split('.sai')[0] + '.sam'
            sam_out_fn = self.sai_fn_dict['single'].split('.sai')[0] + '.samout'
            sam_out_fh = open(sam_out_fn, 'w')
            self.logger.info("\tAligning single-end reads with BWA to %s." % (self.sai_fn_dict[self.read_type]))
            self.logger.info("\tsingle-end reads: %s" % (self.trimmed_fastq_fn_dict['single']) )
            self.logger.info("\tsingle-end index: %s" % (self.sai_fn_dict['single']))
            self.logger.info("\tsam file: %s" % ( sam_fn ))
            bwa_bin = 'bin/bwa'
            bwa_args = [ bwa_bin, 'samse', bwa_idx_prefix, self.sai_fn_dict['single'], self.trimmed_fastq_fn_dict['single'], '-f', sam_fn]
            proc = subprocess.Popen(bwa_args, stdout=sam_out_fh, stderr=sam_out_fh)
            proc.communicate()

            self.sam_fn = sam_fn

        elif self.read_type == 'paired':
            bwa_out_dir = os.path.dirname(self.sai_fn_dict[self.read_type])
            sam_fn = self.sai_fn_dict['paired'].split('.sai')[0] + '.sam'
            sam_out_fn = self.sai_fn_dict['paired'].split('.sai')[0] + '.samout'
            sam_out_fh = open(sam_out_fn, 'w')

            self.logger.info("Aligning paired-end reads with BWA to %s." % (self.sai_fn_dict[self.read_type]))
            self.logger.info("\tsingle-end reads: %s" % ( self.trimmed_fastq_fn_dict['single'] ))
            self.logger.info("\tsingle-end index: %s" % ( self.sai_fn_dict['single'] ))
            self.logger.info("\tpaired-end reads: %s" % ( self.trimmed_fastq_fn_dict['paired'] ))
            self.logger.info("\tpaired-end index: %s" % ( self.sai_fn_dict['paired'] ))
            self.logger.info("\tsam file: %s" % ( sam_fn ))

            bwa_bin = 'bin/bwa'
            bwa_args = [ bwa_bin, 'sampe', bwa_idx_prefix, self.sai_fn_dict['single'], self.sai_fn_dict['paired'],\
                             self.trimmed_fastq_fn_dict['single'], self.trimmed_fastq_fn_dict['paired'], '-f', sam_fn]
            proc = subprocess.Popen(bwa_args, stdout=sam_out_fh, stderr=sam_out_fh)
            proc.communicate()

            self.sam_fn = sam_fn

    def makeNormalizedBedGraph(self, extension_size=200, scaled_library_size=1e7):
        '''
        For many comparisons, it is useful to scale the total read count to an arbitrarily large
        number for all samples.

        This function converts a BAM file into a scaled BedGraph file, which can then be used
        for visualization with the UCSC Genome Browser, or used to represent scaled read count
        in figures.
        '''
        def _checkForBedgraph():
            '''
            Check to see if there is already a BedGraph, in which case we can skip making it.
            '''
            self.logger.info("Checking to see if a BedGraph already exists for this sample...")
            bg_fn = self.bam_fn.replace('.bam', '.scaled.bedGraph')
            if os.path.isfile(bg_fn):
                self.logger.warning("Looks like the BAM file has already been scaled and processed to a BedGraph file...")
                self.bg_fn = bg_fn
                return True
            else:
                self.logger.info("\t...no BedGraph file found.")

        def _setChromSizeFileName():
            self.logger.info("Getting the chromosome sizes for the %s genome..." % (self.genome_version))
            if self.genome_version == 'sacCer3':
                self.chromosome_size_fn = 'reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/ChromInfo.txt'
            else:
                raise ValueError("Error:  invalid genome_version name for %s" % (self.genone_version))

        def _setChromosomeSizes():
            '''
            Retun the chromosome sizes for the given genome_version.
            '''
            self.chromosome_size_dict = {}
            chromosome_size_fh = open(self.chromosome_size_fn)
            for line in chromosome_size_fh:
                line = line.strip().split('\t')
                if line == ['']:
                    continue
                else:
                    self.logger.info("...%s, %s" % (line[0], line[1]) )
                    self.chromosome_size_dict[line[0]] = int(line[1])

            self.logger.info("...done getting chromosome sizes.")

        def _writeMappedLibrarySize():
            self.logger.info("There are %s reads in the mapped library." % (mapped_library_size))

            with open( self.bam_fn.replace('.bam', '.mappedLibrarySize'), 'w') as mapped_library_size_fh:
                mapped_library_size_fh.write(str(mapped_library_size))

        def _makeBEDfromBAM():
            '''
            uses bamToBed program to convert the BAM file to BED format
            line-by-line. We process the stream to extend the coverage
            of each read position by the given extension_size, and then
            write the file to a BED.

            returns the name of the bed file
            '''
            self.bed_fn = self.bam_fn.replace('.bam', '.bed')
            bed_fh = open(self.bed_fn, 'w')
            cmd_args = ['bin/bamToBed', '-i', self.bam_fn]
            proc = subprocess.Popen(cmd_args, stdout=bed_fh, stderr=sys.stderr)
            proc.wait()


        def _makePaddedBEDfromBAM():
            '''
            uses bamToBed program to convert the BAM file to BED format
            line-by-line. We process the stream to extend the coverage
            of each read position by the given extension_size, and then
            write the file to a BED.

            returns the name of the bed file
            '''
            self.bed_fn = self.bam_fn.replace('.bam', '.bed')
            bed_fh = open(self.bed_fn, 'w')
            cmd_args = ['bin/bamToBed', '-i', self.bam_fn]
            proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            bed_PIPE = proc.communicate()[0].split('\n')

            count = 0

            for line in bed_PIPE:
                count += 1
                line = line.split()
                # something the bamToBed program introduces a blank line at the end of the file, this handles it
                try:
                    chrom = line[0]
                except IndexError as error:
                    #print line, count
                    #print error
                    continue
                start = int(line[1])
                stop = int(line[2])
                read_length = stop - start
                name = str(line[3])
                score = str(line[4])
                strand = line[5]
                if strand == '+':
                    stop = start + extension_size
                    if stop > self.chromosome_size_dict[chrom]:
                        self.logger.warning("Chromosome end %s overran while making bed file for sample %s" % (chrom, self.sample_name) )
                        stop = self.chromosome_size_dict[chrom]
                elif strand == '-':
                    stop = stop + read_length
                    start = stop - extension_size
                    if start < 0:
                        self.logger.warning("Chromosome start %s underran while making bed file for sample %s" % (chrom, self.sample_name) )
                        start = 0
                if start > stop:
                    import pdb
                    pdb.set_trace()
                bed_fh.write('\t'.join([chrom,str(start),str(stop),name, score, strand]) + '\n')
            bed_fh.close()

        ####
        # function body starts here for makeNormalizedBedGraph
        ####
        if _checkForBedgraph():
            return

        _setChromSizeFileName()

        try:
            _setChromosomeSizes()
        except ValueError as error:
            self.logger.error(error.value)
            return

        self.logger.info("Making Normalized BedGraph for %s..." % (self.bam_fn))

        mapped_library_size = calculateMappedLibrarySize(self.bam_fn)
        _writeMappedLibrarySize()

        #if extension_size == 0:
        #    self.bed_fn = ''
        if self.experiment_type == 'mnase-seq':
            _makeBEDfromBAM()
        elif self.experiment_type == 'chip-seq':
            _makePaddedBEDfromBAM()
        elif self.experiment_type == 'rna-seq':
            _makeBEDfromBAM()

        sortBED(self.bed_fn)
        self.wig_fn = makeWigFromBED(self.bed_fn, self.chromosome_size_fn)
        self.scaled_wig_fn = scaleWigFile(self.wig_fn, scaled_library_size, mapped_library_size)
        self.scaled_bw_fn = makeBigWigFromWIG(self.scaled_wig_fn, self.chromosome_size_fn)
        self.bg_fn = makeBedGraphFromBigWig(self.scaled_bw_fn)

    def _cleanupSeqFiles(self, fn_suffix_list=['sam', 'sai', 'bw', 'bed', 'wig', 'withduplicates.bam']):
        '''
        Deletes unneeded files to save disk space.

        fn_suffix_list is a list of the file extensions for files that will be removed.

        By default, removes everything but BAM and BedGraph files.
        '''
        fastq_fn_list = glob(self.fastq_dir + '/*.fastq')

        for fastq_fn in fastq_fn_list:
            os.remove(fastq_fn)

        bwa_dir = self.fastq_dir + '/bwa/'

        deleted_fn_list = []
        for fn_suffix in fn_suffix_list:
            deleted_fn_list.extend( glob(bwa_dir + '*.' + fn_suffix) )

        for deleted_fn in deleted_fn_list:
            os.remove(deleted_fn)

class MnaseSampleContainer(SampleContainer):
    def _generateNPSCalls(self,):
        '''
        '''
        self._runNPS()
        self._writeNPSBed()


    def _runNPS(self,):
        '''
        '''
        nps_bin = 'bin/SeqTag.py'
        nps_par_fn = 'processing/NPS.par'

        self.logger.info("Calculating NPS peaks for sample %s" % (self.sample_name) )

        self.nps_out_fn = self.bed_fn.replace('.bed', '.nps')

        nps_log_fn = self.bed_fn.replace('.bed', '.nps.log')
        nps_log_fh = open(nps_log_fn, 'w')

        cmd_args = [nps_bin, nps_par_fn, self.bed_fn, self.nps_out_fn]

        proc = subprocess.Popen(cmd_args, stderr=nps_log_fh, stdout=nps_log_fh)
        proc.wait()

    def _writeNPSBed(self,):
        '''
        Write the NPS file in proper bed format.  Add the center as thick
        '''
        self.nps_bed_fn = self.nps_out_fn.replace('.nps', '.processed.nps')

        nps_fh = open(self.nps_out_fn)
        nps_bed_fh = open(self.nps_bed_fn, 'w')

        # burn the header
        nps_fh.readline()

        self.logger.info("Writing full BED-format file for NPS data to %s..." % (self.nps_bed_fn) )
        for line in nps_fh:
            line = line.strip().split('\t')
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            name = line[3]
            score = line[4]
            strand = '.'
            thickStart = int(start + floor((stop - start) / 2))
            thickEnd = thickStart + 1
            newline = '\t'.join([chrom, str(start), str(stop), name, score, strand, str(thickStart), str(thickEnd)])
            nps_bed_fh.write(newline + '\n')

        nps_fh.close()
        nps_bed_fh.close()

class ChipSampleContainer(SampleContainer):

    def _runMACS(self, macs_version='macs1.4', genome_size='12000000', bam_fn_suffix='.adaptercleaned.bam'):
        '''
        '''
        def _selectMACSbinary():
            if macs_version == 'macs2.0':
                macs_bin = 'bin/macs2'
            if macs_version == 'macs1.4':
                macs_bin = 'bin/macs14'
            if macs_version == 'macs1.3.7':
                macs_bin = 'bin/macs'
            return macs_bin

        def _constructCommandArgs():
            if macs_version == 'macs2.0':
                macs_bin = 'bin/macs2'
                if control_fn:
                    cmd_args = [macs_bin, 'callpeak', '-t', treat_bam_fn, '-c', control_fn,
                                '--name', sample_out_fn_prefix, '--mfold', '1', '50',
                                #'--broad', '--broad-cutoff', '0.01',
                                #'--slocal', '500',
                                #'--llocal', '5000',
                                #'--qvalue', self.macs_qvalue,
                                '--pvalue', '0.001',
                                '--bw='+ str(self.insert_size), '--format=BAM', '--gsize='+genome_size]
                else:
                    cmd_args = [macs_bin, 'callpeak', '-t', treat_bam_fn, '--name',
                                sample_out_fn_prefix, '--format=BAM', '--mfold',
                                '2', '50', '--bw='+ str(self.insert_size), '--gsize='+genome_size]

            if macs_version == 'macs1.4':
                macs_bin = 'bin/macs14'
                if control_fn:
                    cmd_args = [macs_bin, '-t', treat_bam_fn, '-c', control_fn,
                                '--name', sample_out_fn_prefix,
                                #'--mfold=1,50',
                                #'--slocal', '1000',
                                #'--llocal', '10000',
                                '--bw='+ str(self.insert_size),
                                '--pvalue', '0.001',
                                #'--bw=150',
                                '--format=BAM', '--gsize='+genome_size]
                else:
                    cmd_args = [macs_bin, '-t', treat_bam_fn, '--name', sample_out_fn_prefix,\
                                    '--format=BAM', '--mfold=2,50', '--bw='+ str(self.insert_size),\
                                    '--gsize='+genome_size]

            if macs_version == 'macs1.3.7':
                macs_bin = 'bin/macs'
                if control_fn:
                    cmd_args = [macs_bin, '-t', treat_bam_fn, '-c', control_fn, '--name',
                                    sample_out_fn_prefix, '--mfold=2', '--bw='+
                                    str(self.insert_size), '--format=BAM', '--gsize='+genome_size]
                else:
                    cmd_args = [macs_bin, '-t', treat_bam_fn, '--name', sample_out_fn_prefix,\
                                    '--format=BAM', '--mfold=2', '--bw='+ str(self.insert_size),\
                                    '--gsize='+genome_size]

            return cmd_args

        def _makeMacsDir():
            macs_dir = os.path.dirname(self.bam_fn) + '/macs/'
            if not os.path.isdir(macs_dir):
                self.logger.info("Creating directory for MACS output: %s"  % (macs_dir))
                os.mkdir(macs_dir)

            self.macs_dir = macs_dir

        def _getBAMFilename(base_dir):
            dir_str = base_dir + '/FASTQ_' + self.read_type + '_' + self.genome_version + '/bwa/*' + bam_fn_suffix
            tmp_glob = glob(dir_str)
            if len(tmp_glob) > 1:
                print "Error: multiple BAM files found for sample in %s" % (base_dir)
                sys.exit()
            elif len(tmp_glob) == 0:
                self.logger.error("Could not find BAM files for sample %s with %s." % (self.sample_name, dir_str) )
            else:
                return tmp_glob[0]

        def _addLocalPythonPath():
            env = os.environ.copy()
            if 'PYTHONPATH' in env:
                env['PYTHONPATH'] += ':lib/python2.7/site-packages/'
            else:
                env['PYTHONPATH'] = 'lib/python2.7/site-packages/'
            return env

        #####
        ## runMacs function body starts here
        ####
        self.macs_version = macs_version

        treat_bam_fn = _getBAMFilename(self.sample_dir)

        _makeMacsDir()

        if self.control_data_dir != 'none':
            control_fn = _getBAMFilename(self.control_data_dir)
        else:
            control_fn = None

        macs_log_fn = self.macs_dir + '/' + self.sample_name + '.macs.log'
        macs_log_fh = open(macs_log_fn, 'w')

        sample_out_fn_prefix = self.macs_dir + '/' + self.sample_name

        if control_fn:
            self.logger.info("Running MACS for sample: %s" % (self.sample_name))
            self.logger.info("\ttreatment BAM file: %s" % (treat_bam_fn))
            self.logger.info("\tcontrol BAM file: %s" % (control_fn))
            #cmd_args = [macs_bin, '-t', treat_bam_fn, '-c', control_fn, '--name', sample_out_fn_prefix, '--bw='+ self.insert_size, '--format=BAM', '--gsize='+genome_size]
            cmd_args = _constructCommandArgs()
        else:
            self.logger.info("Running MACS for sample: %s" % (self.sample_name))
            self.logger.info("\ttreatment BAM file: %s" % (treat_bam_fn))
            self.logger.info("\tcontrol BAM file: %s" % ('none'))
            #cmd_args = [macs_bin, '-t', treat_bam_fn, '--name', sample_out_fn_prefix, '--format=BAM', '--bw='+ self.insert_size, '--gsize='+genome_size]
            cmd_args = _constructCommandArgs()

        macs_env = _addLocalPythonPath()

        proc = subprocess.Popen(cmd_args, stdout=macs_log_fh, stderr=macs_log_fh, env=macs_env)
        proc.wait()

        if macs_version == 'macs2.0':
            self.macs_peaks_bed_fn = sample_out_fn_prefix + '_peaks.narrowPeak'
        else:
            self.macs_peaks_bed_fn = sample_out_fn_prefix + '_peaks.bed'
        self.macs_summits_bed_fn = sample_out_fn_prefix + '_summits.bed'
        self.logger.info('...done.')

    def _runPeakSplitter(self,):
        self.logger.info("Splitting peaks for sample in %s" % (self.sample_dir))

        # skip files that have no peaks
        macs_peaks_bed_fh = open(self.macs_peaks_bed_fn)
        if macs_peaks_bed_fh.readline() == "":
            return

        # if no track line, add it to the bedGraph file
        bg_fh = open(self.bg_fn)
        first_line = bg_fh.readline()
        if not first_line.startswith('track'):
            new_bg_fn = self.bg_fn + '.tmp'
            new_bg_fh = open(new_bg_fn, 'w')
            new_bg_fh.write('track type=bedGraph\n')
            while True:
                line = bg_fh.readline()
                if line == '':
                    break
                else:
                    new_bg_fh.write(line)
            new_bg_fh.close()
            bg_fh.close()
            os.rename(new_bg_fn, self.bg_fn)

        print "Splitting peaks for %s..." % (self.macs_peaks_bed_fn)
        peaksplitter_bin='bin/PeakSplitter'
        cmd_args = [peaksplitter_bin, '-p', self.macs_peaks_bed_fn, '-w', self.bg_fn, '-o', self.macs_dir, '-f']

        proc = subprocess.Popen(cmd_args, stdout=sys.stdout, stderr=sys.stderr)
        proc.wait()

        self.subpeaks_bed_fn = self.macs_peaks_bed_fn.replace('.bed', '.subpeaks.bed')

    def _postProcessPeakSplitter(self,):
        self.logger.info("Post-processing subpeaks file %s:" % (self.subpeaks_bed_fn) )

        # this has the subpeak summits too
        subpeaks_bed_fh = open(self.subpeaks_bed_fn)

        subpeaks_summits_bed_fn = self.macs_summits_bed_fn.replace('summits.bed', 'subpeaks.summits.bed')
        subpeaks_summits_bed_fh = open(subpeaks_summits_bed_fn, 'w')

        subpeak_counter = 0
        for line in subpeaks_bed_fh:
            # skip header
            if line.startswith('Chromosome'):
                continue
            subpeak_counter += 1
            line = line.strip().split('\t')
            chrom = line[0]
            start = line[1]
            stop = line[2]
            if self.macs_version == 'macs2.0':
                score = line[4]
            else:
                score = line[3]
            summit_start = line[4]
            summit_stop = str(int(summit_start) + 1)
            name = 'Subpeak_' + str(subpeak_counter)

            out_line = "%s\t%s\t%s\t%s\t%s\n" % (chrom, summit_start, summit_stop, name, score)
            #out_line = "%s\t%s\t%s\t%s\n" % (chrom, summit_start, summit_stop, name)
            subpeaks_summits_bed_fh.write(out_line)

        subpeaks_bed_fh.close()
        subpeaks_summits_bed_fh.close()

