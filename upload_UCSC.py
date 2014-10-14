#!/usr/bin/python
import sys
import os
import glob
import shutil

def add_header_to_bed(fn, dest, sample_attrs_dict):
    '''
    .bed files from MACS need the UCSC header line added in order to track
    properly.  This add a minimal header using the file name.
    '''
    print 'Adding header to the %s...' % (fn)

    assert 'sample_name' in sample_attrs_dict
    assert 'sample_desc' in sample_attrs_dict

    if not 'color' in sample_attrs_dict:
        sample_attrs_dict['color'] = '255,0,0'
    if not 'visibility' in sample_attrs_dict:
        sample_attrs_dict['visibility'] = 'full'

    header = 'track type=%(track_type)s name="%(sample_name)s" description="%(sample_desc)s" visibility="%(visibility)s" color=%(color)s autoScale=on maxHeightPixels=100:24:21 graphType=bar alwaysZero=on yLineOnOff=on gridDefault=on windowingFunction=mean\n' % sample_attrs_dict
    print '\t' + header

    in_fh = open(fn, 'r')
    lines = in_fh.readlines()
    #if lines[0].startswith('track'):
    #    print "This file already has a header...not modifying it."
    #    return fn
    in_fh.close()
    out_fn = dest + '/' + os.path.basename(fn)
    print 'Writing updated file to %s' % (out_fn)
    out_fh = open(out_fn, 'w')
    out_fh.write(header)
    out_fh.writelines(lines[1:])
    print '...done.'

    return out_fn

def upload_file(hgsid, track_type, sample_attrs_dict, dest='/home/matt/public_html/ucsc/', local_server_addr='http://dounce.icmb.utexas.edu/~matt/ucsc/'):
    '''
    Copies the given filename to a tmpfile in the http directory, uploads it 
    with wget, and deletes the tmpfile when done.
    '''
    relative_fn = os.path.relpath(sample_attrs_dict['fn'])

    new_fn = add_header_to_bed(relative_fn, dest, sample_attrs_dict)
    local_url = local_server_addr + os.path.basename(new_fn)

    remote_server_addr = 'http://genome.ucsc.edu/cgi-bin/hgCustom?hgsid=' 
    remote_url = remote_server_addr + hgsid + '&hgct_customText=' + local_url

    #cmd = "wget --tries=1 --quiet -O /dev/null '" + remote_url + "'"
    cmd = "curl -I '" + remote_url + "'"
    print '#######', remote_url, cmd
    print "Pushing to \n\t%s \nfrom \n\t%s with command \n\t%s" % (remote_url, local_url, cmd)
    os.system(cmd)
    print "\nCompleted push."

def usage():
    if (len(sys.argv) < 3):
        print "Usage:  python upload_UCSC.py [filename] [hgsid];"
        print "optional: trackname, samplename, sampledesc, color"
        print "can be assigned with a csv list (i.e. samplename,'Sample A')."
        sys.exit()    

if __name__ == '__main__':
    '''
    usage()

    fn = sys.argv[1]
    hgsid = sys.argv[2]
    sample_attrs_dict = dict([arg.split(',', 1) for arg in sys.argv[3:]])
    sample_attrs_dict['ftype'] = fn.split('.')[-1]
    upload_file(fn, hgsid, sample_attrs_dict)

    '''

    ucsc_id = sys.argv[1]
    track_type = sys.argv[2]
    data_type = sys.argv[3]

    if data_type == 'chip':
        top_data_dir = 'data/chipData/'
    elif data_type == 'mnase':
        top_data_dir = 'data/mnaseData/'

    data_dir_list = glob.glob(top_data_dir + '/Sample*index1/')

    for data_dir in data_dir_list:
        sample_attrs_dict = {}
        sample_id_fh = open(data_dir + '/SampleID')
        for line in sample_id_fh:
            if track_type not in line:
                continue
            line = line.split(';')
            assert len(line) == 2  # only one attr per line
            key = line[0].strip()
            attr = line[1].strip()
            attr.replace(' ', '_')
            key = key.replace('ucsc_' + track_type + '_', '')
            sample_attrs_dict[key] = attr

        try:

            if track_type == 'coverage' and data_type == 'chip':
                sample_attrs_dict['fn'] = glob.glob(data_dir + '/FASTQ_single_sacCer3/bwa/*bedGraph')[0]
            elif track_type == 'peaks' and data_type == 'chip':
                sample_attrs_dict['fn'] = glob.glob(data_dir + '/FASTQ_single_sacCer3/bwa/macs/*_peaks.bed')[0]
            elif track_type == 'peaks' and data_type == 'mnase':
                sample_attrs_dict['fn'] = glob.glob(data_dir + '/FASTQ_paired_sacCer3/bwa/*processed.nps')[0]
            elif track_type == 'coverage' and data_type == 'mnase':
                sample_attrs_dict['fn'] = glob.glob(data_dir + '/FASTQ_paired_sacCer3/bwa/*scaled.bedGraph')[0]
            elif track_type == 'centers' and data_type == 'mnase':
                sample_attrs_dict['fn'] = glob.glob(data_dir + '/FASTQ_paired_sacCer3/bwa/*pe_centers.bedGraph')[0]
            elif track_type == 'density' and data_type == 'mnase':
                sample_attrs_dict['fn'] = glob.glob(data_dir + '/FASTQ_paired_sacCer3/bwa/*pe_density.bedGraph')[0]

            if track_type == 'coverage':
                sample_attrs_dict['track_type'] = 'bedGraph'
            elif track_type == 'centers':
                sample_attrs_dict['track_type'] = 'bedGraph'
            elif track_type == 'density':
                sample_attrs_dict['track_type'] = 'bedGraph'
            elif track_type == 'peaks':
                sample_attrs_dict['track_type'] = 'bed'

        except IndexError:
            continue

        upload_file(ucsc_id, track_type, sample_attrs_dict)

