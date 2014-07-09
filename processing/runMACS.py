#!/usr/bin/python
import sys
import os
import subprocess
from glob import glob

def runMacs(sample_metadata_dict, macs_version='macs1.4', genome_size='12e6', bam_fn_suffix='.adaptercleaned.bam'):
    '''
    '''
    def _selectMACSbinary():
        if macs_version == 'macs2.0':
            pass
        if macs_version == 'macs1.4':
            macs_bin = 'bin/macs14'
        return macs_bin

    def _makeMacsDir():
        macs_dir = os.path.dirname(sample_metadata_dict['bam_fn']) + '/macs/'
        if not os.path.isdir(macs_dir):
            print "Creating directory for MACS output: %s."  % (macs_dir)
            os.mkdir(macs_dir)

        return macs_dir

    def _getBAMFilename(base_dir):
        dir_str = base_dir + '/FASTQ_' + sample_metadata_dict['read_type'] + '_' + sample_metadata_dict['genome_version'] + '/bwa/*' + bam_fn_suffix
        tmp_glob = glob(dir_str)
        if len(tmp_glob) > 1:
            print "Error: multiple BAM files found for sample in %s" % (base_dir)
            sys.exit()
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
    sample_dir = sample_metadata_dict['sample_dir']

    treat_bam_fn = _getBAMFilename(sample_dir)

    if sample_metadata_dict['control_data_dir'] != 'none':
        control_fn = _getBAMFilename(sample_metadata_dict['control_data_dir'])
    else:
        control_fn = None

    macs_dir = _makeMacsDir()
    macs_log_fn = macs_dir + '/' + sample_metadata_dict['sample_name'] + '.macs.log'
    macs_log_fh = open(macs_log_fn, 'w')

    sample_out_fn = macs_dir + '/' + sample_metadata_dict['sample_name']

    macs_bin = _selectMACSbinary()
    if control_fn: 
        print "Running MACS for sample: %s" % (sample_metadata_dict['sample_name'])
        print "\ttreatment BAM file: %s" % (treat_bam_fn)
        print "\tcontrol BAM file: %s" % (control_fn)
        cmd_args = [macs_bin, '-t', treat_bam_fn, '-c', control_fn, '-n', sample_out_fn, '-fBAM', '-g', genome_size]
    else:
        print "Running MACS for sample: %s" % (sample_metadata_dict['sample_name'])
        print "\ttreatment BAM file: %s" % (treat_bam_fn)
        print "\tcontrol BAM file: %s" % ('none')
        cmd_args = [macs_bin, '-t', treat_bam_fn, '-n', sample_out_fn, '-fBAM', '-g', genome_size]

        macs_env = _addLocalPythonPath()

    proc = subprocess.Popen(cmd_args, stdout=macs_log_fh, stderr=macs_log_fh, env=macs_env)
    proc.wait()
    print '...done.'


if __name__ == '__main__':
    '''
    sample_dirs = glob(sys.argv[1] + '/*')
    for sample_dir in sample_dirs:
        runMacs(sample_dir, macs_bin='/g/software/bin/macs14', genome_size='12e6')
    '''
    pass
