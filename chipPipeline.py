#!/usr/bin/python
import sys
sys.path.append('./lib/python2.7/site-packages/')
import os
from tools import *
from processing.sampleObjects import *

def printUsage():
    '''
    Print the usage message for the program.
    '''
    logger.error("Usage:  python ./pipeline.py <sample_dir> <read_type> <genome_version>")
    sys.exit()

def processSample(sample_dir):
    '''
    '''
    sampleContainer = ChipSampleContainer(sample_dir)

    # gunzip and concatenate the fq.gz files
    # fqs is a dictionary with the concatenated filename for
    # single and paired files
    sampleContainer.gunzipFastqs()

    # trim the adapter sequences from the reads
    # capture the filenames in the return to pass
    # to the mapping function of your choice
    sampleContainer.trimFastqFiles()  # updates sample_metadata_dict['trimmed_fastq_fn_dict']

    sampleContainer.mapBWA(max_cpu=2, min_mapq_score=37, remove_duplicate_reads_flag=True)

    sampleContainer.makeNormalizedBedGraph(extension_size=200)

    sampleContainer._runMACS(macs_version='macs1.4')
    sampleContainer._runPeakSplitter()
    sampleContainer._postProcessPeakSplitter()

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


def runAllSamples():
    all_chip_sample_dir_list = glob('data/chipData/Sample_LL2_index*')
    control_chip_sample_dir_list = []
    treatment_chip_sample_dir_list = []

    for chip_sample_dir in all_chip_sample_dir_list:
        sample_metadata_dict = parseSampleMetaData(chip_sample_dir)
        if sample_metadata_dict['control_data_name'] == 'none':
            control_chip_sample_dir_list.append(chip_sample_dir)
        else:
            treatment_chip_sample_dir_list.append(chip_sample_dir)

    sample_dir_list = []
    sample_dir_list.extend(control_chip_sample_dir_list)
    sample_dir_list.extend(treatment_chip_sample_dir_list)

    runSamples(sample_dir_list=sample_dir_list, max_cpu=4)



if __name__ == '__main__':

    runSamples(sample_dir_list=glob('data/chipData/Sample_LL2_index*'), max_cpu=4)
    #processSample(sample_dir='data/chipData/Sample_LL2_index3/')



