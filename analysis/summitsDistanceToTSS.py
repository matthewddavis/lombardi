#!/g/software/bin/python2.7
import sys
import numpy as np
from glob import glob

class ChromVectorDict(dict):
    def __init__(self, chrom_info_fn='reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/ChromInfo.txt'):
        self.chrom_info_fn = chrom_info_fn
        self._makeChromVectorDict()
        self._loadTSS()

    def _makeChromVectorDict(self,):
        chrom_info_fh = open(self.chrom_info_fn)
        for line in chrom_info_fh:
            line = line.strip().split('\t')
            self[line[0]] = np.zeros(int(line[1]), dtype=int)

    def _loadTSS(self, genes_fn='reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/sgdGene.txt'):
        genes_fh = open(genes_fn)
        lines = genes_fh.readlines()

        for line in lines[1:]:
            line = line.strip().split('\t')
            chrom = line[2]
            tx_start = int(line[4])
            self[chrom][tx_start] = True

    def _shuffleChromVectors(self,):
        for chrom in self:
            np.random.shuffle(self[chrom])


def loadSummitsList():
    summits_fn_list = []
    summits_fn_list.extend(glob('../data/chipData/Sample_LL2_index2/FASTQ_single_sacCer3/bwa/macs/*_summits.bed'))
    summits_fn_list.extend(glob('../data/chipData/Sample_LL2_index3/FASTQ_single_sacCer3/bwa/macs/*_summits.bed'))
    summits_fn_list.extend(glob('../data/chipData/Sample_LL2_index22/FASTQ_single_sacCer3/bwa/macs/*_summits.bed'))
    summits_fn_list.extend(glob('../data/chipData/Sample_LL2_index21/FASTQ_single_sacCer3/bwa/macs/*_summits.bed'))
    summits_fn_list.extend(glob('../data/chipData/Sample_LL2_index4/FASTQ_single_sacCer3/bwa/macs/*_summits.bed'))
    summits_fn_list.extend(glob('../data/chipData/Sample_LL2_index8/FASTQ_single_sacCer3/bwa/macs/*_summits.bed'))

                       
    summits_list = []
    
    for summit_fn in summits_fn_list:
        print "Loading summits from %s..." % summit_fn
        summit_fh = open(summit_fn)
        for line in summit_fh:
            line = line.strip().split('\t')
            chrom = line[0]
            summit_pos = int(line[1])
            summits_list.append((chrom, summit_pos))

    return summits_list

def calcDistArray(summits_list, chrom_vector_dict):
    summit_dists = []
    summit_counter = 0
    for summit in summits_list:
        summit_counter += 1
        if summit_counter % 100 == 0:
            print summit_counter
        chrom = summit[0]
        pos = summit[1]
        dist = min(abs(np.where(chrom_vector_dict[chrom])[0] - [pos]))
        summit_dists.append(dist)
    
    dists_a = np.array(summit_dists, dtype=int)
        
    return dists_a

def reportPercentage(dists_a, cutoff=450, out_fn='output/percent_summits_from_TSS.txt'):
    out_fh = open(out_fn, 'w')
    percent = len(np.where(dists_a <= cutoff)[0]) / float(len(dists_a)) * 100.0
    
    out_fh.write('%2.2f percent of summits are within %s' % (percent, cutoff) )
    print '%2.2f percent of summits are within %s' % (percent, cutoff) 
    

if __name__ == '__main__':
    chrom_vector_dict = ChromVectorDict()
    summits_list = loadSummitsList()
    dists_array = calcDistArray(summits_list, chrom_vector_dict)
    reportPercentage(dists_array, out_fn='output/percent_summits_from_TSS.txt')

    chrom_vector_dict._shuffleChromVectors()
    dists_array = calcDistArray(summits_list, chrom_vector_dict)
    reportPercentage(dists_array, out_fn='output/percent_summits_from_TSS_shuffled.txt')

