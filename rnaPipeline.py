#!/usr/bin/python
import sys
sys.path.append('./lib/python2.7/site-packages/')
import os
from processing import gunzipFastqs
from processing.bwaMapper import *
from processing.runTophat import *
from tools import *

def printUsage():
    '''
    Print the usage message for the program.
    '''
    print "Usage:  python ./pipeline.py <sample_dir> <read_type> <genome_version>"
    sys.exit()

def cleanupSeqFiles(sample_metadata_dict, fn_suffix_list=['sam', 'sai', 'bw', 'bed', 'wig', 'withduplicates.bam']):
    '''
    Deletes unneeded files to save disk space.

    fn_suffix_list is a list of the file extensions for files that will be removed.

    By default, removes everything but BAM and BedGraph files.
    '''
    fastq_dir = os.path.dirname(sample_metadata_dict['combined_fastq_fn_dict']['single'])
    fastq_fn_list = glob(fastq_dir + '/*.fastq')

    for fastq_fn in fastq_fn_list: 
        os.remove(fastq_fn)

    bwa_dir = fastq_dir + '/bwa/'
    
    deleted_fn_list = []
    for fn_suffix in fn_suffix_list:
        deleted_fn_list.extend( glob(bwa_dir + '*.' + fn_suffix) )

    for deleted_fn in deleted_fn_list:
        os.remove(deleted_fn)

def processSample(sample_dir):
    '''
    '''    
    sample_metadata_dict = parseSampleMetaData(sample_dir)
    sample_metadata_dict['sample_dir'] = sample_dir

    print "Processing sample with files in %s..." % sample_metadata_dict['sample_dir']
    print "\tread type: %s" % sample_metadata_dict['read_type']
    print "\tgenome version: %s" % sample_metadata_dict['genome_version']
    #print "\tsingle-end adapter sequence: %s" % sample_metadata_dict['adapter_se']
    #print "\tpaired-end adapter sequence: %s" % sample_metadata_dict['adapter_pe']
    #print "\tinsert size: %s" % sample_metadata_dict['insert_size']
    
    # gunzip and concatenate the fq.gz files
    # fqs is a dictionary with the concatenated filename for 
    # single and paired files
    gunzipFastqs(sample_metadata_dict)  # updates sample_metadata_dict['combined_fastq_fn_dict']

    # trim the adapter sequences from the reads
    # capture the filenames in the return to pass
    # to the mapping function of your choice
    trimFastqFiles(sample_metadata_dict)  # updates sample_metadata_dict['trimmed_fastq_fn_dict'] 

    runTophat(sample_metadata_dict, segment_length=18, max_cpu=1)                                                                                                                     
    runCufflinks(sample_metadata_dict, max_cpu=1)  

    cleanupSeqFiles(sample_metadata_dict, fn_suffix_list=['sam', 'sai', 'bw', 'bed', 'wig', 'withduplicates.bam'])

def runSamples(sample_dir_list, max_cpu=1):
    '''
    '''
    from multiprocessing import Pool
    from glob import glob

    def _initMultiprocPool():
        '''
        Figure out how many CPUs to use in the multiproc pool,
        initialize the pool, and then return the pool.
        '''
        sample_num = len('sample_dir_list')
        if sample_num > max_cpu:
            pool = Pool(processes=max_cpu)
        else:
            pool = Pool(processes=sample_num)
        
        return pool

    pool = _initMultiprocPool()
    
    pool.map(processSample, sample_dir_list)

if __name__ == '__main__':

    runSamples(sample_dir_list=glob('data/mnaseData/Sample_LL1_index*'))
    processSample(sample_dir='data/snyderRNASeq/original_dt/')


