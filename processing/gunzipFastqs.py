#!/usr/bin/python
import os
from glob import glob
import gzip
import logging

def gunzipFastqs(sample_metadata_dict, fastq_dir_prefix='FASTQ', single_substr='R1', paired_substr='R2', combined_fastq_fn_suffix='.all.fastq'):
    '''
    Unzips and concatenates the FASTQ files from an Illumina experiment
    for downstream analysis.

    Returns a dict of the names of the single and paired end file names for other functions to use.

    The concatenated file is stored in a subdirectory of the
    sample directory.

    The substr identifies the read file as single or paired end.
    '''
    def _getFastqgzFilenameList():
        fastqgz_fn_list = glob(sample_dir + '/' + '*fastq.gz')

        if fastqgz_fn_list == []:
            raise ValueError("Error finding the FASTQ files in the sample directory %s... exiting." % (sample_dir))
        else:
            sample_metadata_dict['logger'].info("Found gzipped fastq files in %s..." % (sample_dir))
            return fastqgz_fn_list

    def _checkForCombinedFastq():
        '''
        Check for existing combine fastq files.

        If they exist, return them inthe combined_fastq_fn_dict,
        if not, 
        '''
        print "Checking for previously unzipped fastq files..."

        #combined_fastq_fn_dict = { 'single' : None,
        #                           'paired' : None,
        #                           }
        combined_fastq_fn_dict = {}
        sample_metadata_dict['combined_fastq_fn_dict'] = combined_fastq_fn_dict

        single_combined_fastq_fn = fastq_dir + '/' + sample_name + '_' + single_substr + combined_fastq_fn_suffix
        paired_combined_fastq_fn = fastq_dir + '/' + sample_name + '_' + paired_substr + combined_fastq_fn_suffix

        if os.path.isfile(single_combined_fastq_fn):
            sample_metadata_dict['logger'].info("...found the single-end combined fastq file.")
            combined_fastq_fn_dict['single'] = single_combined_fastq_fn
        else:
            return False

        if os.path.isfile(paired_combined_fastq_fn):
            sample_metadata_dict['logger'].info("...found the paired-end combined fastq file.")
            combined_fastq_fn_dict['paired'] = paired_combined_fastq_fn

        return True

    def _combineFastqgzFiles():
        '''
        Uncompresses and concatenates gzipped fastq files.
        '''
        # get the list of fastq.gz files. This catches ValueError and exits if files don't exist.
        try:
            fastqgz_fn_list = _getFastqgzFilenameList()
            print "Found fastq.gz files to combine in %s..." % (sample_dir)
        except ValueError as error:
            print error
            sys.exit()

        print "...creating a new FASTQ directory for this analysis at %s." % (fastq_dir)
        try:
            os.mkdir(fastq_dir)
        except:
            print "...oh, wait, it was already there, but without combined fastq files."
            

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
                print "#####"
                print "## WARNING:  Skipping %s because it doesn't contain a single or paired end FASTQ substring!" % (fastqgz_fn)
                print "#####"

        combined_fastq_fn_dict = {}

        if len(single_fastqgz_fn_list) > 0:
            single_combined_fastq_fn = fastq_dir + '/' + sample_name + '_' + single_substr + combined_fastq_fn_suffix
            os.system('zcat ' + ' '.join(single_fastqgz_fn_list) + '> ' + single_combined_fastq_fn)
            combined_fastq_fn_dict['single'] = single_combined_fastq_fn
        if len(paired_fastqgz_fn_list) > 0:
            paired_combined_fastq_fn = fastq_dir + '/' + sample_name + '_' + paired_substr + combined_fastq_fn_suffix
            os.system('zcat ' + ' '.join(paired_fastqgz_fn_list) + '> ' + paired_combined_fastq_fn)
            combined_fastq_fn_dict['paired'] = paired_combined_fastq_fn

        sample_metadata_dict['combined_fastq_fn_dict'] = combined_fastq_fn_dict

    #####
    ## function body beings here
    #####

    # make my typing easier
    sample_dir = sample_metadata_dict['sample_dir']
    read_type = sample_metadata_dict['read_type']
    genome_version = sample_metadata_dict['genome_version']
    sample_name = sample_metadata_dict['sample_name']

    # set the output directory
    fastq_dir = sample_dir + '/' + '_'.join([fastq_dir_prefix, read_type, genome_version])

    # check for existing unzipped files. If the they exist, the 
    # combined_fastq_fn_dict is updated and we return
    if _checkForCombinedFastq():
        return

    _combineFastqgzFiles()

if __name__ == "__main__":
    print "Module not intended to by called directly."

