#!/usr/bin/python
import sys
sys.path.append('./lib/python2.7/site-packages/')
import os
from processing import gunzipFastqs
from processing.removeAdapters import *
from processing.bwaMapper import *
from processing.runNPS import *
from tools import *
from processing.sampleObjects import *

def processSample(sample_dir):
    '''
    '''    
    sampleContainer = MnaseSampleContainer(sample_dir)

    # gunzip and concatenate the fq.gz files
    # fqs is a dictionary with the concatenated filename for 
    # single and paired files
    sampleContainer.gunzipFastqs()

    # trim the adapter sequences from the reads
    # capture the filenames in the return to pass
    # to the mapping function of your choice
    sampleContainer.trimFastqFiles()  # updates sample_metadata_dict['trimmed_fastq_fn_dict'] 

    sampleContainer.mapBWA(max_cpu=2, remove_duplicate_reads_flag=False)

    sampleContainer.makeNormalizedBedGraph(extension_size=0)

    sampleContainer._generateNPSCalls()

    sampleContainer._cleanupSeqFiles(fn_suffix_list=['sam', 'sai', 'bw', 'bed', 'wig', 'withduplicates.bam'])

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

    #runSamples(sample_dir_list=glob('data/mnaseData/Sample_LL1_index*'), max_cpu=4)
    processSample(sample_dir='data/mnaseTest/Sample_LL1_index99/')



