#!/usr/bin/python
import os
import subprocess
import tarfile
import urllib
from multiprocessing import cpu_count
from glob import glob
import gzip
from Bio import SeqIO
from mappedReadsTools import *

def getMaxCpu(max_cpu):
    '''
    # set the CPU number for threading
    '''
    cpu_num = cpu_count() - 1
    if cpu_num > max_cpu:
        cpu_num = max_cpu

    return max_cpu

def selectFastqType(sample_metadata_dict):
    '''
    The Snyder RNA-Seq reads don't get trimmed, and so we need
    to cheat that step without breaking the rest of the code.
    '''
    # used untrimmed reads if rna-seq
    if sample_metadata_dict['experiment_type'] == 'rna-seq':
        trimmed_fastq_fn_dict = sample_metadata_dict['combined_fastq_fn_dict']
    elif sample_metadata_dict['experiment_type'] == 'mnase-seq':
        trimmed_fastq_fn_dict = sample_metadata_dict['trimmed_fastq_fn_dict']
    elif sample_metadata_dict['experiment_type'] == 'chip-seq':
        trimmed_fastq_fn_dict = sample_metadata_dict['trimmed_fastq_fn_dict']
    else:
        print "#####"
        print "## ERROR: experiment_type invalid: %s" % (sample_metadata_dict['experiment_type'] )
        print "#####"

    return trimmed_fastq_fn_dict

def mapBWA(sample_metadata_dict, max_cpu=1, remove_duplicate_reads=True):
    '''
    Wrapper fucntion for mapping reads with BWA (as opposed to something else).
    '''
    read_type = sample_metadata_dict['read_type']
    genome_version = sample_metadata_dict['genome_version']

    trimmed_fastq_fn_dict = selectFastqType(sample_metadata_dict)

    # check for existing BAM file
    test_bam_fn = os.path.dirname(trimmed_fastq_fn_dict[read_type]) + '/' + 'bwa' + '/' + os.path.basename(trimmed_fastq_fn_dict[read_type]).split('.fastq')[0] + '.bam'
    if os.path.isfile(test_bam_fn):
        print "######"
        print "## WARNING:  BAM file already exists!  %s" % (test_bam_fn)
        print "######"
        sample_metadata_dict['bam_fn'] = test_bam_fn

        return

    ## if no BAM file found, we make it
    # first make the index
    indexReadsBWA(sample_metadata_dict)  # updates sample_metadata_dict['sai_fn_dict']
    # then align to it
    alignBWA(sample_metadata_dict)

    # then turn the SAM file into a BAM file
    sample_metadata_dict['bam_fn'] = makeBAMfromSAM(sample_metadata_dict['sam_fn'], remove_duplicate_reads=remove_duplicate_reads)

def prepBWA(genome_version, bwa_index_name='bwaIndex'):
    '''
    Ensures index for genome_version is built, or builds it.
    '''

    def _checkForBWAIdx():
        if os.path.isfile(bwa_idx_prefix + '.bwt'):
            print "Found an existing BWA Index at %s..." % (bwa_idx_prefix)
            return True
        else:
            print "No existing BWA Index found, making it..."
            return False

    def _prepReferenceGenome():
        '''
        download genome from hard-coded location and untar into reference/ location

        note: this function has a lot of hard-coded values for file names. Recreation
              using a different reference genome will require recoding these variables.
        '''

        def _makeChromSizesFile():
            chrom_sizes_fn = genome_version_dict[genome_version]['untar_path'] + 'chrom.sizes'
            chrom_sizes_fh = open(chrom_sizes_fn, 'w')
            chrom_seq_dict = SeqIO.to_dict(SeqIO.parse(genome_version_dict[genome_version]['genome_fasta_fn'], format='fasta'))

            for chrom in chrom_seq_dict:
                chrom_sizes_fh.write('%s\t%s\n' % ( chrom, len(chrom_seq_dict[chrom]) ))
            chrom_sizes_fh.close()

        # get the genome tarball and untar it
        if not os.path.isdir(genome_version_dict[genome_version]['untar_path']):
            os.mkdir(genome_version_dict[genome_version]['untar_path'])

        urllib.urlretrieve(genome_version_dict[genome_version]['remote_genome_location'], filename=genome_version_dict[genome_version]['local_genome_tar_fn'])

        try:
            genome_tarfile = tarfile.open(genome_version_dict[genome_version]['local_genome_tar_fn'], 'r|gz')
            genome_tarfile.extractall(path=genome_version_dict[genome_version]['untar_path'])
        except tarfile.ReadError:
            # if it's just a gz file and not tarred, we should end up here
            in_fh = gzip.open(genome_version_dict[genome_version]['local_genome_tar_fn'])
            out_fh = open(genome_version_dict[genome_version]['genome_fasta_fn'], 'w')
            out_fh.write( in_fh.read() )
            in_fh.close()
            out_fh.close()

        # stupid UCSC doesn't come with a combined genome fasta file, so we make one
        if genome_version == 'sacCer3':
            whole_genome_fasta_fn = genome_version_dict[genome_version]['untar_path'] + genome_version + '.fa'
            # if it exists, delete it and make it again
            if os.path.isfile(whole_genome_fasta_fn):
                os.remove(whole_genome_fasta_fn)

            chrom_fasta_fn_list = glob(genome_version_dict[genome_version]['untar_path'] + '*fa')

            whole_genome_fasta_fh = open(whole_genome_fasta_fn, 'w')

            for chrom_fasta_fn in chrom_fasta_fn_list:
                chrom_fasta_fh = open(chrom_fasta_fn)
                whole_genome_fasta_fh.writelines( chrom_fasta_fh.readlines() )
                chrom_fasta_fh.close()
            whole_genome_fasta_fh.close()

            _makeChromSizesFile()

            return  genome_version_dict[genome_version]['genome_fasta_fn']

        elif genome_version == 'sacCer3masked':
            whole_genome_fasta_fn = genome_version_dict[genome_version]['untar_path'] + genome_version + '.fa'
            # if it exists, delete it and make it again
            if os.path.isfile(whole_genome_fasta_fn):
                os.remove(whole_genome_fasta_fn)

            chrom_fasta_fn_list = glob(genome_version_dict[genome_version]['untar_path'] + '*fa.masked')

            whole_genome_fasta_fh = open(whole_genome_fasta_fn, 'w')

            for chrom_fasta_fn in chrom_fasta_fn_list:
                chrom_fasta_fh = open(chrom_fasta_fn)
                whole_genome_fasta_fh.writelines( chrom_fasta_fh.readlines() )
                chrom_fasta_fh.close()
            whole_genome_fasta_fh.close()

            _makeChromSizesFile()

            return  genome_version_dict[genome_version]['genome_fasta_fn']
        elif genome_version == 'S288':

            _makeChromSizesFile()

            return  genome_version_dict[genome_version]['genome_fasta_fn']

    def _buildIndexBWA():
        '''
        Creates the BWA index files to align with.
        '''
        # create the bwaIndex dir if it doesn't exist
        bwa_bin = 'bin/bwa'
        bwa_args = [bwa_bin, 'index', '-a', 'is', '-p', bwa_idx_prefix, genome_version_dict[genome_version]['genome_fasta_fn']]

        bwa_proc = subprocess.Popen(bwa_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        bwa_proc.communicate()

    #####
    ## Function body starts here
    #####
    genome_version_dict = { 'S288C' : { 'remote_genome_location' : 'http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz',
                                        'local_genome_tar_fn' : 'reference/genome/S288C_reference_genome_R64-1-1_20110203.tgz',
                                        'untar_path' : 'reference/genome/',
                                        'genome_fasta_fn' : 'reference/genome/S288C_reference_genome_R64-1-1_20110203/S288C_reference_sequence_R64-1-1_20110203.fsa',
                                        },
                            'sacCer3' : { 'remote_genome_location' : 'http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz',
                                          'local_genome_tar_fn' : 'reference/genome/sacCer3/chromFa.tar.gz',
                                          'untar_path' : 'reference/genome/sacCer3/',
                                          'genome_fasta_fn' : 'reference/genome/sacCer3/sacCer3.fa',
                                          },
                            'sacCer3masked' : { 'remote_genome_location' : 'http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/chromFaMasked.tar.gz',
                                                'local_genome_tar_fn' : 'reference/genome/sacCer3masked/chromFaMasked.tar.gz',
                                                'untar_path' : 'reference/genome/sacCer3masked/',
                                                'genome_fasta_fn' : 'reference/genome/sacCer3masked/sacCer3masked.fa',
                                                },
                            'sacCer3rna' : { 'remote_genome_location' : 'http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/mrna.fa.gz',
                                             'local_genome_tar_fn' : 'reference/genome/sacCer3rna/mrna.fa.gz',
                                             'untar_path' : 'reference/genome/sacCer3rna/',
                                             'genome_fasta_fn' : 'reference/genome/sacCer3rna/mrna.fa',
                                             },
                            }


    bwa_idx_prefix = genome_version_dict[genome_version]['untar_path'] + '/' + bwa_index_name
    if _checkForBWAIdx():
        return

    genome_fasta_fn = _prepReferenceGenome()
    _buildIndexBWA()

def indexReadsBWA(sample_metadata_dict, single_substr='R1', paired_substr='R2', max_cpu=1):
    '''
    Index the reads with BWA.

    the file name of the BWA index directory is hardcoded in bwa_idx_dirname
    '''
    from glob import glob
    from multiprocessing import cpu_count
    genome_version = sample_metadata_dict['genome_version']
    trimmed_fastq_fn_dict = selectFastqType(sample_metadata_dict)

    if genome_version == 'sacCer3':
        bwa_idx_prefix = 'reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BWAIndex/genome.fa'
    else:
        print "Error:  invalid genome_version name for %s" % (genome_version)
        sys.exit()

    sai_fn_dict = {}
    for trimmed_fastq_fn_key in trimmed_fastq_fn_dict:
        # make a directory to store things in
        path = os.path.dirname(trimmed_fastq_fn_dict[trimmed_fastq_fn_key])
        bwa_out_dir=path + '/bwa'

        if not os.path.isdir(bwa_out_dir):
            os.mkdir(bwa_out_dir)

        sai_fn = bwa_out_dir + '/' + os.path.basename(trimmed_fastq_fn_dict[trimmed_fastq_fn_key]).split('.fastq')[0] + '.sai'
        if sai_fn.find(single_substr) != -1:
            sai_fn_dict['single'] = sai_fn
        elif sai_fn.find(paired_substr) != -1:
            sai_fn_dict['paired'] = sai_fn
        sai_out_fn = bwa_out_dir + '/' + os.path.basename(trimmed_fastq_fn_dict[trimmed_fastq_fn_key]).split('.fastq')[0] + '.saiout'
        sai_out_fh = open(sai_out_fn, 'w')
        print "Making BWA idx file %s." % (sai_fn)

        max_cpu = getMaxCpu(max_cpu)

        bwa_bin = 'bin/bwa'
        bwa_idx_fn = bwa_idx_prefix + 'bwaIndex'
        bwa_args = [bwa_bin, 'aln', '-t', str(max_cpu), bwa_idx_prefix, trimmed_fastq_fn_dict[trimmed_fastq_fn_key], '-f', sai_fn]
        proc = subprocess.Popen(bwa_args, stdout=sai_out_fh, stderr=sai_out_fh)
        proc.communicate()
        sai_out_fh.close()

    sample_metadata_dict['sai_fn_dict'] = sai_fn_dict

def alignBWA(sample_metadata_dict, single_substr='R1', paired_substr='R2', max_cpu=1):
    '''
    Align reads with BWA.
    '''
    read_type = sample_metadata_dict['read_type']
    sai_fn_dict = sample_metadata_dict['sai_fn_dict']
    trimmed_fastq_fn_dict = selectFastqType(sample_metadata_dict)
    genome_version = sample_metadata_dict['genome_version']

    if genome_version == 'sacCer3':
        bwa_idx_prefix = 'reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BWAIndex/genome.fa'
    else:
        print "Error:  invalid genome_version name for %s" % (genome_version)
        sys.exit()
    max_cpu = getMaxCpu(max_cpu)

    if read_type == 'single':
        print 'Aligning single-end reads...'
        bwa_out_dir = os.path.dirname(sai_fn_dict['single'])
        sam_fn = sai_fn_dict[read_type].split('.sai')[0] + '.sam'
        sam_out_fn = sai_fn_dict['single'].split('.sai')[0] + '.samout'
        sam_out_fh = open(sam_out_fn, 'w')
        print "\tAligning single-end reads with BWA to %s." % (sai_fn_dict[read_type])
        print "\t\tReads from %s \n\t\twith %s for the index..." % (trimmed_fastq_fn_dict['single'], sai_fn_dict['single'])
        bwa_bin = 'bin/bwa'
        bwa_args = [ bwa_bin, 'samse', bwa_idx_prefix, sai_fn_dict['single'], trimmed_fastq_fn_dict['single'], '-f', sam_fn]
        proc = subprocess.Popen(bwa_args, stdout=sam_out_fh, stderr=sam_out_fh)
        proc.communicate()

        sample_metadata_dict['sam_fn'] = sam_fn

    elif read_type == 'paired':
        bwa_out_dir = os.path.dirname(sai_fn_dict[read_type])
        sam_fn = sai_fn_dict['paired'].split('.sai')[0] + '.sam'
        sam_out_fn = sai_fn_dict['paired'].split('.sai')[0] + '.samout'
        sam_out_fh = open(sam_out_fn, 'w')

        print "Aligning paired-end reads with BWA to %s." % (sai_fn_dict[read_type])
        print "\tsingle-end reads: %s" % ( trimmed_fastq_fn_dict['single'] )
        print "\tsingle-end index: %s" % ( sai_fn_dict['single'] )
        print "\tpaired-end reads: %s" % ( trimmed_fastq_fn_dict['paired'] )
        print "\tpaired-end index: %s" % ( sai_fn_dict['paired'] )
        print "\tsam file: %s" % ( sam_fn )

        bwa_bin = 'bin/bwa'
        bwa_args = [ bwa_bin, 'sampe', bwa_idx_prefix, sai_fn_dict['single'], sai_fn_dict['paired'],\
                         trimmed_fastq_fn_dict['single'], trimmed_fastq_fn_dict['paired'], '-f', sam_fn]
        proc = subprocess.Popen(bwa_args, stdout=sam_out_fh, stderr=sam_out_fh)
        proc.communicate()

        sample_metadata_dict['sam_fn'] = sam_fn
        
