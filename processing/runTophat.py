#!/usr/bin/python
import sys
import os
import subprocess
from glob import glob
from multiprocessing import cpu_count


def addLocalBinPath():
    env = os.environ.copy()
    env['PATH'] += ':bin/'
    return env

def getMaxCpu(max_cpu):
    system_cpu_count = cpu_count()
    if system_cpu_count > max_cpu:
        return str(max_cpu)
    elif system_cpu_count <= max_cpu:
        return str(system_cpu_count - 1)

def runTophat(sample_metadata_dir, segment_length='18', bowtie_version='bowtie2', max_cpu=1):

    def _makeTophatOutdir(tophat_out_dirname):
        if not os.path.isdir(tophat_out_dirname):
            os.mkdir(tophat_out_dirname)
        return tophat_out_dirname

    fastq_dirname = '%(sample_dir)s/FASTQ_%(read_type)s_%(genome_version)s' % (sample_metadata_dir)
    tophat_out_dirname = _makeTophatOutdir(fastq_dirname + '/tophat/')

    fastq_fn_list = glob(fastq_dirname + '/*fastq')

    print "Mapping reads from %s to yeast transcriptome using Tophat." % (fastq_dirname)

    tophat_bin = 'bin/tophat'

    if bowtie_version == 'bowtie2':
        bowtie_idx_prefix='reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome'
    elif bowtie_version == 'bowtie':
        bowtie_idx_prefix='reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/BowtieIndex/genome'
    else:
        print "#####"
        print "## Error: bowtie_version should be either 'bowtie' or 'bowtie2'"
        print "#####"
        
        return 
        
    gtf_fn='reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/sgdGene.gtf'
    cpu_num = getMaxCpu(max_cpu)

    tophat_log_fn = tophat_out_dirname + '/tophat.log'
    tophat_log_fh = open(tophat_log_fn, 'w')

    local_env = addLocalBinPath()
    
    cmd_args = [tophat_bin, '-p', cpu_num, '--no-novel-juncs', '--segment-length', str(segment_length), '-o',\
                    tophat_out_dirname, '-G', gtf_fn, bowtie_idx_prefix, ' '.join(fastq_fn_list)]
    proc = subprocess.Popen(cmd_args, stderr=tophat_log_fh, stdout=tophat_log_fh, env=local_env)
    proc.wait()

    if not os.path.isfile(tophat_out_dirname + '/accepted_hits.bam'):
        print "#####"
        print "## Error Tophat failed. See log file %s for details." % (tophat_log_fn)
        return

    sample_metadata_dir['accepted_hits_fn'] = tophat_out_dirname + '/accepted_hits.bam'

def runCufflinks(sample_metadata_dict, max_cpu=1):
    '''
    '''
    input_bam_fn = sample_metadata_dict['accepted_hits_fn']

    cufflinks_bin = 'bin/cufflinks'
    cufflinks_out_dir = sample_metadata_dict['sample_dir'] + '/cufflinks/'
    gtf_fn = 'reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/sgdGene.gtf'

    cufflinks_log_fn = cufflinks_out_dir + '/cufflinks.log'
    cufflinks_log_fh = open(cufflinks_log_fn, 'w')

    local_env = addLocalBinPath()

    print "Running Cufflinks to assign FPKM for sample %s." % (sample_metadata_dict['sample_name'])

    cpu_num = getMaxCpu(max_cpu)
    cmd_args = [ cufflinks_bin, '-p', cpu_num, '-G', gtf_fn, '-o', cufflinks_out_dir, input_bam_fn]

    proc = subprocess.Popen(cmd_args, stderr=cufflinks_log_fh, stdout=cufflinks_log_fh, env=local_env)
    proc.wait()


