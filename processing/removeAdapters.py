#!/usr/bin/python
import sys
import os

def trimFastqFiles(sample_metadata_dict, trimmed_fastq_fn_suffix='.adaptercleaned.fastq'):
    '''
    Trims given adapter sequences out of the fastq files.
    '''
    def _checkForTrimmedFastqFn():
        '''
        '''
        if os.path.isfile(trimmed_fastq_fn):
            print "#####"
            print "## WARNING:  It looks like the adapter sequences were already trimmed from fastq file %s !" % (trimmed_fastq_fn)
            print "#####"
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

            print "\tFrom a total of %s, %s were trimmed (%s percent)." % (read_count, trim_count, percent)

            combined_fastq_fh.close()
            trimmed_fastq_fh.close()
            trim_counts_fh.close()

        def _removeAdapter(read, adapter_seq):
            '''
            Removes an adapter_sequence from a read by taking a 3-nt match and extending it.
            '''
            #print read
            idx = read.rfind(adapter_seq)
            if idx != -1:
                #print 'Found the adapter at pos %s in read %s' % (idx, read)
                return read[0:idx]
            else:
                idx = read.rfind(adapter_seq[0:3])
                if idx != -1:
                    pad = len(read) - idx
                    idx = read.rfind(adapter_seq[0:pad])
                    if idx != -1:
                        #print 'Found the incomplete adapter at pos %s in read %s' % (idx, read)
                        return read[0:idx]
                    else:
                        #print 'Not even incomplete adapter found in read %s' % (read)
                        return read
                else:
                    #print 'No adapter found in read %s' % (read)
                    return read

        #####
        ## trimming function body starts here
        #####
        combined_fastq_fh = open(combined_fastq_fn_dict[combined_fastq_fn_key])
        trimmed_fastq_fn = sample_metadata_dict['trimmed_fastq_fn_dict'][combined_fastq_fn_key]
        trimmed_fastq_fh =  open(trimmed_fastq_fn, 'w')
        trim_counts_fn = trimmed_fastq_fn.split('.')[0] + '.trimcounts'
        trim_counts_fh = open(trim_counts_fn, 'w')
        
        if combined_fastq_fn_key == 'single':
            adapter_seq = sample_metadata_dict['adapter_se']
        elif combined_fastq_fn_key == 'paired':
            adapter_seq = sample_metadata_dict['adapter_pe']

        print "\nRemoving adapter sequences from %s..." % (trimmed_fastq_fn)

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
                print read_count, read

    #####
    ## main function starts here
    #####
    # make my typing easier
    combined_fastq_fn_dict = sample_metadata_dict['combined_fastq_fn_dict']
    #trimmed_fastq_fn_dict = { 'single' : None,
    #                          'paired' : None,
    #                          }
    
    trimmed_fastq_fn_dict = {}

    sample_metadata_dict['trimmed_fastq_fn_dict'] = trimmed_fastq_fn_dict

    for combined_fastq_fn_key in combined_fastq_fn_dict:
        trimmed_fastq_fn = combined_fastq_fn_dict[combined_fastq_fn_key].split('.')[0] + trimmed_fastq_fn_suffix
        if _checkForTrimmedFastqFn():
            trimmed_fastq_fn_dict[combined_fastq_fn_key] = trimmed_fastq_fn
            continue
        trimmed_fastq_fn_dict[combined_fastq_fn_key] = trimmed_fastq_fn
        trimFastqFile(combined_fastq_fn_key)



