#!/g/software/bin/python2.7
'''
Written by Matt Davis.

Library for analysis of MNase data containing classes and methods
employed in the analyzeNucs file.
'''
import sys
sys.path.append('../lib/python2.7/site-packages/')
import numpy as np
from glob import glob
from tools import *
from matplotlib import use
use('Agg')
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

class MnaseDataSamplesContainer(dict):
    '''
    Container for nucleosome annotations by gene.

    The dictionary is keyed by the name of the sample, and the value of each
        key is a nucDataSampele object.
    '''
    def __init__(self):
        pass

    def annotateSample(self, mnaseDataSample):
        '''
        Annotates a sample for nucleosome positions near each gene.
        '''
        sampleID_fn = mnaseDataSample.sample_dir + '/SampleID'

        print "Annotating nucleosomes for %s sample..." % (mnaseDataSample.sample_name)
        nps_fn_list = glob(mnaseDataSample.sample_dir + '/FASTQ*/bwa/*processed.nps')
        #nps_fn_list = glob(mnaseDataSample.sample_dir + '/bwa*/*processed.bed')
        assert len(nps_fn_list) == 1
        mnaseDataSample.nps_fn = nps_fn_list[0]
        mnaseDataSample.annotateAllNucsTSS()
        self[mnaseDataSample.sample_name] = mnaseDataSample
        print "...done."

    def annotateChipBindingAllSamples(self, chipDataSamplesContainer, feature='TSS', dist5p=500, dist3p=300):
        '''
        Annotates all the mnaseDataSamples with binding information from the chipDataSamplesContainer
        provided.
        '''
        for mnaseDataSample in self:
            self[mnaseDataSample].annotateChipBinding(chipDataSamplesContainer, feature=feature, dist5p=dist5p, dist3p=dist3p)


class MnaseSampleDict(dict):
    '''
    A modified dictionary with additional attributes for the nucleosome
    data.

    Keyed by genes, which are incorporated from a reference file such as
    the Pugh 2009 file.
    '''
    def __init__(self, pugh_genes_reference_fn):
        self.pugh_genes_reference_fn = pugh_genes_reference_fn
        self.real_nuc_count = 0
        self.fake_nuc_count = 0
        self.dataDir = ''
        self.nps_fn = ''
        self.nps_data_dict = {}
        self.annotation_model = 0
        self.nucs_vector_dict = {}
        self.annotated_nuc_dists_dict = {}
        self.chip_masked_nuc_dists_dict = {}
        self.annotateNoneChipBinding()
        self.bedGraph_fn = ''
        self.restrict_to_verified = False
        self.chromosome_size_fn = 'reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/ChromInfo.txt'

    def loadMetadata(self, sample_dir):
        '''
        Loads metadata from SampleID file to object attributes.

        Also adds the sample_dir to the object attributes.
        '''
        sample_metadata_dict = parseSampleMetaData(sample_dir)
        for key, value in sample_metadata_dict.items():
            setattr(self, key, value)
        self.sample_dir = sample_dir

    def loadChromScoreVectorDict(self):
        '''
        Accepts a bedGraph filename and loads all chromosomes into float vectors
        of the scores across the bedGraph regions.
        '''

        def _setChromosomeSizes():
            '''
            '''
            self.chromosome_size_dict = {}
            chromosome_size_fh = open(self.chromosome_size_fn)
            for line in chromosome_size_fh:
                line = line.strip().split('\t')
                if line == ['']:
                    continue
                else:
                    print "...%s, %s" % (line[0], line[1]) 
                    self.chromosome_size_dict[line[0]] = int(line[1])

            print "...done getting chromosome sizes."

        try:
            bg_fh = open(self.bedGraph_fn, 'r')
        except:
            print "Error opening bedGraph file %s." % (self.bedGraph_fn)

        print "Loading bedGraph scores from %s..." % (self.bedGraph_fn)

        _setChromosomeSizes()
        chrom_score_vector_dict = {}

        for chrom in self.chromosome_size_dict:
            chrom_score_vector_dict[chrom] = np.zeros(self.chromosome_size_dict[chrom], dtype=float)

        for line in bg_fh:
            if line.startswith('track'):
                continue
            line = line.strip('\n').split('\t')
            chrom = line[0]
            start = int(line[1]) - 1  # bed format is 1-indexed
            stop = int(line[2]) - 1
            score = float(line[3])
            if chrom in chrom_score_vector_dict:
                chrom_score_vector_dict[chrom][start:stop] = score

        print "\t...done."
        self.chrom_score_vector_dict = chrom_score_vector_dict


    def makeNucsDists(self, annotated_nuc_pos_list, chip_sample_name_list=['none']):
        '''
        Add the nucleosome distances for each nucleosome position for each gene in the sample.
        '''
        print "Calculating the annotated nucleosome position distances for sample %s..." % self.sample_name

        for annotated_nuc_pos in annotated_nuc_pos_list:
            if annotated_nuc_pos not in self.annotated_nuc_dists_dict:
                self.annotated_nuc_dists_dict[annotated_nuc_pos] = {}
            for chip_sample_name in chip_sample_name_list:
                self.annotated_nuc_dists_dict[annotated_nuc_pos][chip_sample_name] = np.array( [self[gene_name][annotated_nuc_pos]['pos']
                                                                                                for gene_name in self.keys()
                                                                                                if annotated_nuc_pos in self[gene_name]
                                                                                                and self[gene_name]['in_set'][chip_sample_name] ])
    def annotateAllNucsTSS(self):
        '''
        Calls function annotateNucsByGene for all genes in the verified genes file.
        The strategy is to first identify the +1 nuc, then add all the pos nucs,
        then add the neg nucs, and finally identify the NFR.
        '''

        ##### top level function ####
        #
        #
        print "\tAnnotating all nucleosome positions for sample %s from the file %s..." % (self.sample_name, self.nps_fn)

        # nucs from NPS processed bed file
        self.readNPSbed()
        # dict of vectors of nuc positions keyed by chromosome
        self.buildNucsVectors()

        count = 0
        for gene in self:
            count += 1
            self.annotatePosNucsByGene(gene)
            self.annotateNegNucsByGene(gene)
            self.annotateNFRByGene(gene)

        print "\t...done annotating nucleosome positions for %s genes." % (count)

    def readNPSbed(self):
        '''
        Reads the processed NPS bed file with the centers of the nucleosomes annotated.
        Returns a dict of the rows from the bed file keyed by chromosome.
        '''

        try:
            npsBedfh = open(self.nps_fn, 'r')
        except:
            pdb.set_trace()
            print "Error opening NPS bed file:  %s..." % (self.nps_fn)
            return

        # we'll store the NPS data in self.nps_data_dict = {}
        # This is a dict keyed by chrom

        for line in npsBedfh:
            line = line.strip().split('\t')
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            name = line[3]
            score = float(line[4])
            strand = line[5]
            center1 = int(line[6])
            center2 = int(line[7])
            if not chrom in self.nps_data_dict:
                self.nps_data_dict[chrom] = []
                self.nps_data_dict[chrom].append( [chrom, start, stop, name, score, strand, center1, center2] )
            else:
                self.nps_data_dict[chrom].append( [chrom, start, stop, name, score, strand, center1, center2] )


    def findClosestNuc(self, diffVector, gene):
        '''
        diffVector - the difference vector between the genomic reference feature (e.g. TSS or TTS)
                     and the vector of nuc position centers.
        gene - the NucSampleDict key, used to assess strand orientation of the gene
        '''
        closestNuc = {}
        # finds the index and value of the closest nuc center to the TSS
        closestNuc['idx'] = abs(diffVector).argmin()
        closestNuc['pos'] = diffVector[closestNuc['idx']]
        closestNuc['coord'] = self.nucs_vector_dict[gene['chrom']][closestNuc['idx']]

        # if the gene is on Crick, then invert the sign of the distance
        if gene['Strand'] == '-':
            closestNuc['pos'] = closestNuc['pos'] * -1

        return closestNuc

    def annotatePosNucsByGene(self, gene_name, maxDist = 200):
        '''
        Annotates all nucleosomes associated with all genes
        in the specified set from +1 to the end of the TSS.

        Accepts the dictionary of genes from Pugh and the
        dictionary of nucleosomes by chromosome from NPS.
        Returns nothing but modifies the genes dictionary in place by adding the
        nuc annotations for the NPS dataset being considered.
        '''
        # find the closest position to the TSS
        # this difference vector is used to find the closest nucs with findClosestNuc function
        diffVector = self.nucs_vector_dict[self[gene_name]['chrom']] - float(self[gene_name]['TSS'])

        # first we have to find the +1 nucleosome, then we can find the others relative
        # to the +1, if any of them exist

        while 'AnP1' not in self[gene_name]:
            closestNuc = self.findClosestNuc(diffVector, self[gene_name])
            # if the closest nucleosome is more than maxDist away, then no annotation is made
            if abs(closestNuc['pos']) > maxDist: break
            # Now that we know the nuc is the closest, check to see if it qualifies as -1,0, or +1 position
            nucAnno = self.checkProximalNucpos(closestNuc['pos'])
            #print closestNuc, "###", nucAnno

            ######
            # if the closest nuc is a +1, then annotate it as such and continue annotation with
            #    annotateMorePosNucs;
            # if it is not, then effectively remove it from the search for the +1 by setting
            #    its position to the largest distance in the diffVector, then continue in the
            #    while loop looking for the +1; the diffVector is recalculated in the annotateMorePosNucs
            #    function, so no need to worry about the original value for the overwritten nuc distance
            ######
            if nucAnno == '1':
                self[gene_name]['AnP1'] = { 'pos' : closestNuc['pos'],
                                 'idx' : closestNuc['idx'],
                                 'coord' : closestNuc['coord'],
                                 }

                # the next function will attempt to fill in all positive nuc positions
                # using the +1 as a reference
                self.annotateMorePosNucs(gene_name)
            else:
                diffVector[closestNuc['idx']] = diffVector.argmax()

    def annotateMorePosNucs(self, gene_name):
        '''
        Once the +1 nuc position is known for a gene, annotate all subsequent
        nucs in the locus. The idea is to observe where a nucleosome could have
        been displaced, and then replace it with a hypothetical nucleosome, so that we
        can properly annotate a +3 in the abscence of a +2, for example. This is more
        of a problem deeper into loci, and is an idea borrowed from Pugh 2009, though
        it is modeled more realistically here (Pugh just distributes nucleosomes evenly
        throughout lengthy spacer regions).

        gene - the locus of interest
        nucs_vector_dict - the vector of nuc positions, by chromosome
        model - which model to use for inserting a hypothetical nucleosome; see code
        '''

        # the list of nuc keys for the annotatedNucs dict
        # assumes no more than 1000 positive nucs for a locus
        poslist = ['AnP' + str(i) for i in range(2,1000)]

        chrom = self[gene_name]['chrom']
        # recalculate the diffVector to include all nucs original values
        diffVector = self.nucs_vector_dict[chrom] - self[gene_name]['TSS']

        # the initial positions are all the +1 nucleosome
        # recall the pos is the position relative to the TSS in nucleotides,
        # the idx is the index in the vector of nuc positions, and
        # the coord is the genomic coordinate of the center of the nuc
        prevpos = self[gene_name]['AnP1']['pos']
        previdx = self[gene_name]['AnP1']['idx']
        prevcoord = self[gene_name]['AnP1']['coord']

        # farBound is the farthest the next nuc can be without needing to be skipped for a hypothetical
        # spacer nucleosome (see Pugh 2009)
        while True:
            # print '#### ', sample[gene_name]['GeneID'], ' ####'

            ##
            # There are four models for how much space intervenes in the case of a skipped nucleosome; default is Model Four.
            # Each considers the 73bp to move from the center to the edge of the nuc,
            # then a given max linker, plus another 73 to get to
            # the middle of the next nucleosome
            # Model One:  73 + 147; assumes center of nucleosome should not be more than half a wrap away
            # Model Two: 73 + 147 + 73; assumes next nuc does not leave space for 147bp to be wrapped around another nuc
            # Model Three: 73 + 20 + 147 + 20 + 73; assumes that nothing should change if a nucleosome is removed,
            #              including that there should still be approximately 20bp linkers on each side of the missing
            #              nuc
            # Model Four: 73 + 0 + 147 + 0 + 73; assumes that nothing should change if a nucleosome is removed,
            #              except that there should not necessarily be 20bp linkers on each side of the missing nuc
            ##
            if self.annotation_model == 1:
                farBound = prevpos + (73 + 147)
            elif self.annotation_model == 2:
                farBound = prevpos + (73 + 147 + 73)
            elif self.annotation_model == 3:
                farBound = prevpos + (73 + 20 + 147 + 20 + 73)
            elif self.annotation_model == 4:
                farBound = prevpos + (73 + 0 + 147 + 0 + 73)
            elif self.annotation_model == 5:
                farBound = prevpos + (73 + 0 + 147 + 0 + 73) - 40
            else:
                print "!!! ERROR: Annotation Model set to invalid value: %s" % self.annotation_model

            # strand orientation affects the arithmetic
            if self[gene_name]['Strand'] == '+':
                # go until we get past the TTS
                if ( (farBound + self[gene_name]['TSS']) > self[gene_name]['TTS']):
                    #print "breaking", farBound, gene['TSS'], gene['TTS']
                    break

                # compare to the next nuc in the nucVector, if we're at the end of the
                # chromosome, then break
                idx = previdx + 1
                if idx == len(diffVector):
                    break
                # the next nuc position...
                pos = diffVector[idx]
            elif self[gene_name]['Strand'] == '-':
                # go until we get past the TTS
                if ( (self[gene_name]['TSS'] - farBound) < self[gene_name]['TTS']):
                    #print "breaking", farBound, gene['TSS'], gene['TTS']
                    break
                # strand logic is reversed
                idx = previdx - 1
                if idx < 0:
                    break
                pos = diffVector[idx] * -1

            # if the position was within the farBound, then there's no need to skip,
            # so annotate the nuc with the next position in the list of nucs for this locus
            if pos <= farBound:
                self[gene_name][poslist.pop(0)] = {'pos' : pos,
                                                   'idx' : idx,
                                                   'coord' : self.nucs_vector_dict[chrom][idx],
                                                   }
                prevpos = pos
                previdx = idx
                #print "real nuc at", pos, idx, prevpos, farBound
                self.real_nuc_count += 1

            else:
                # we're inserting a hypothetical nuc, but we don't give it an actual position;
                # the place holder allows us to continue the numbering scheme properly
                # but to calculate the next distance, we want to pretend that hypo nuc
                # really exists and had a center, so that we can calculate the next distance
                # for each model, we position the hypo nuc in the middle of the skipped linker
                self[gene_name][poslist.pop(0)] = { 'pos' : np.nan,
                                                    'idx' : idx,
                                                    'coord' : self.nucs_vector_dict[chrom][idx],
                                                    }
                # for clarity, create a hypothetic position, then set prevpos
                if self.annotation_model == 1:
                    hypopos = prevpos + int((73 + 147) / 2.0)
                    prevpos = hypopos
                elif self.annotation_model == 2:
                    hypopos = prevpos + int((73 + 147 + 73) / 2.0)
                    prevpos = hypopos
                elif self.annotation_model == 3:
                    hypopos = prevpos + int((73 + 20 + 147 + 20 + 73) / 2.0)
                    prevpos = hypopos
                elif self.annotation_model == 4:
                    hypopos = prevpos + int((73 + 0 + 147 + 0 + 73) / 2.0)
                    prevpos = hypopos
                elif self.annotation_model == 5:
                    hypopos = prevpos + int((73 + 0 + 147 + 0 + 73 - 40) / 2.0)
                    prevpos = hypopos
                else:
                    print "!!! ERROR: Annotation Model set to invalid value: %s" % self.annotation_model
                previdx = idx
                #print "fake nuc at", pos, idx, prevpos, farBound
                self.fake_nuc_count += 1

    def annotateNFRByGene(self, gene_name, maxDist = 200):
        '''
        Annotates a nuc in the NFR region of a locus, if it exists.
        Accepts the dictionary of genes from Pugh and the
        dictionary of nucleosomes by chromosome from NPS.

        Returns nothing but modifies the genes dictionary in place by adding the
        nuc annotations for the NPS dataset being considered.
        '''

        # find the closest position to the TSS
        diffVector = self.nucs_vector_dict[self[gene_name]['chrom']] - float(self[gene_name]['TSS'])

        while 'AnP0' not in self[gene_name]:
            closestNuc = self.findClosestNuc(diffVector, self[gene_name])
            if abs(closestNuc['pos']) > maxDist: break
            nucAnno = self.checkProximalNucpos(closestNuc['pos'])
            #print closestNuc, "###", nucAnno

            if nucAnno == '0':
                self[gene_name]['AnP0'] = {'pos' : closestNuc['pos'],
                                'idx' : closestNuc['idx'],
                                'coord' : closestNuc['coord'],
                                }
                self.real_nuc_count += 1
            else:
                diffVector[closestNuc['idx']] = diffVector.argmax()
    def annotateNegNucsByGene(self, gene_name, maxDist=310):
        '''
        Annotates the -1 nucleosome.
        Accepts the dictionary of genes from Pugh and the
        dictionary of nucleosomes by chromosome from NPS.
        Returns nothing but modifies the genes dictionary in place by adding the
        nuc annotation for the NPS dataset being considered.
        Could be modified to retun additional negative nucleosomes, but not implemented.
        '''

        # find the closest position to the TSS
        diffVector = self.nucs_vector_dict[self[gene_name]['chrom']] - float(self[gene_name]['TSS'])

        while 'AnN1' not in self[gene_name]:
            closestNuc = self.findClosestNuc(diffVector, self[gene_name])
            if abs(closestNuc['pos']) > maxDist: break
            nucAnno = self.checkProximalNucpos(closestNuc['pos'])
            #print closestNuc, "###", nucAnno

            if nucAnno == '-1':
                self[gene_name]['AnN1'] = {'pos' : closestNuc['pos'],
                                'idx' : closestNuc['idx'],
                                'coord' : closestNuc['coord'],
                                }
                self.real_nuc_count += 1
            else:
                diffVector[closestNuc['idx']] = diffVector.argmax()
    def buildNucsVectors(self):
        '''
        Construct dictionary of arrays of nucleosomes by chromsome.
        '''
        chroms = []
        for gene in self:
            chroms.append(self[gene]['chrom'])
        chroms = set(chroms)

        # dictionary for the array of centers of nucleosomes keyed by chromosome
        # self.nucs_vector_dict
        for chrom in chroms:
            # get an array of values for the chromosome of interest

            self.nucs_vector_dict[chrom] = np.array([self.nps_data_dict[chrom][i][6] for i in range(len(self.nps_data_dict[chrom]))])

    def checkProximalNucpos(self, closestpos):
        '''
        Checks to see if the nucleosome is a -1,0,1 nucleosome.
        These distances are from Pugh 2009 Genome Biol.
        '''
        if closestpos <= -111 and closestpos >= -307:
            return '-1'
        elif closestpos <= -6 and closestpos >= -110:
            return '0'
        elif closestpos <= 144 and closestpos >= -5:
            return '1'
        else:
            # return False if the closest nuc doesn't qualify as any of -1,0,1
            return False

    def annotateFPKM(self, fpkm_fn='../data/snyderRNASeq/original_dt/cufflinks/genes.fpkm_tracking',
                     top_percent_cutoff=90.0, bottom_percent_cutoff=10.0):
        '''
        '''
        fpkm_fh = open(fpkm_fn)

        gene_fpkm_dict = {}
        for line in fpkm_fh:
            if line.startswith('tracking_id'):
                    continue
            line = line.strip().split('\t')
            gene = line[3]
            fpkm = float(line[9])
            gene_fpkm_dict[gene] = fpkm

        for gene in self:
            if gene in gene_fpkm_dict:
                self[gene]['fpkm'] = gene_fpkm_dict[gene]
            else:
                self[gene]['fpkm'] = np.nan

        fpkm_vector = np.array([self[gene]['fpkm'] for gene in self])
        print "Loading gene expresion data for %s genes from %s" % (len(fpkm_vector), fpkm_fn)
        fpkm_vector = fpkm_vector[(np.isnan(fpkm_vector) == False)]
        print "A total of %s loci have measureable expression." % (len(fpkm_vector))

        top_cut_val = scoreatpercentile(fpkm_vector, top_percent_cutoff)
        print "using a top FPKM cutoff of %s, %s (np.log2)" % (top_cut_val, np.log2(top_cut_val))
        bottom_cut_val = scoreatpercentile(fpkm_vector, bottom_percent_cutoff)
        print "using a bottom FPKM cutoff of %s, %s (np.log2)" % (bottom_cut_val, np.log2(bottom_cut_val))

        for gene in self:
            if 'in_set' not in self[gene]:
                self[gene]['in_set'] = {}
            if self[gene]['fpkm'] > top_cut_val:
                self[gene]['in_set']['fpkm_top'] = True
            else:
                self[gene]['in_set']['fpkm_top'] = False
            if self[gene]['fpkm'] < bottom_cut_val:
                self[gene]['in_set']['fpkm_bottom'] = True
            else:
                self[gene]['in_set']['fpkm_bottom'] = False

    def annotateMatchedFPKMSet(self,):
        pass


    def annotateNoneChipBinding(self,):
        '''
        '''
        for gene in self:
            if not 'in_set' in self[gene]:
                self[gene]['in_set'] = {}
                self[gene]['in_set']['none'] = True

    def annotateChipBinding(self, chipDataSamplesContainer, feature='TSS', dist5p=500, dist3p=300):
        '''
        For each nucleosome annotation dataset, this function annotates each gene in the
        dataset as bound or unbound (True | False), given there is or is not a ChIP summit
        within the specified distance of the TSS.
        '''
        for chip_data_sample_name in chipDataSamplesContainer:
            print "Annotating %s with ChIP binding from ChIP sample %s..." % (self.sample_name, chip_data_sample_name)

            chipDataSamplesContainer[chip_data_sample_name].loadChromScoreVectorDict()
            chromScoreVectorDict = chipDataSamplesContainer[chip_data_sample_name].chrom_score_vector_dict

            # list to track the summits included
            summits_included = []

            for gene in self:
                chrom = self[gene]['chrom']

                if not 'chip_enrichment' in self[gene]:
                    self[gene]['chip_enrichment'] = {}
                if not chip_data_sample_name in self[gene]['chip_enrichment']:
                    self[gene]['chip_enrichment'][chip_data_sample_name] = np.nan
                strand = self[gene]['Strand']
                if strand == '+':
                    start = self[gene][feature] - dist5p
                    stop = self[gene][feature] + dist3p
                elif strand == '-':
                    start = self[gene][feature] - dist3p
                    stop = self[gene][feature] + dist5p
                self[gene]['chip_enrichment'][chip_data_sample_name] = np.mean(chromScoreVectorDict[chrom][start:stop])

                self[gene]['in_set'][chip_data_sample_name] = False

                '''
                summitNp.Array = np.array([i[1] for i in chipDataSamples[chipDataSample].summits if i[0] == chrom])
                testVal = np.array([i[1] for i in chipDataSamples[chipDataSample].summits if i[0] == chrom]) - self[gene][feature]
                # summit_dists.append(min(abs(testVal)))
                summitValues = np.array([i[2] for i in chipDataSamples[chipDataSample].summits if i[0] == chrom])
                #if gene == 'YBR009C':
                #if gene == 'YGL007C-A':
                #test = (testVal < 0) & (testVal > (-1*distance))
                locus_length = abs(self[gene]['TTS'] - self[gene]['TSS'])
                if strand == '+':
                    #test = ((testVal) > -1*dist5p) & ((testVal) < dist3p) & ((testVal) < locus_length)
                    test = (testVal > -1*dist5p) & (testVal < dist3p)
                elif strand == '-':
                    testVal = testVal * -1
                    #test = ((testVal) > -1*dist5p) & ((testVal) < dist3p) & ((testVal) < locus_length)
                    test = (testVal > -1*dist5p) & (testVal < dist3p)
                #test = abs(testVal) < distance
                '''
                summit_array = np.array([i[1] for i in chipDataSamplesContainer[chip_data_sample_name].summits if i[0] == chrom])
                gene_feature_pos = self[gene][feature]
                dist_to_feature_array = summit_array - gene_feature_pos
                gene_length = abs(self[gene]['TTS'] - self[gene]['TSS'])
                if strand == '+':
                    #nearby_summit_idx = (dist_to_feature_array > (-1 * dist5p)) & (dist_to_feature_array < dist3p)
                    nearby_summit_idx = (dist_to_feature_array > (-1 * dist5p)) & (dist_to_feature_array < dist3p) & (dist_to_feature_array < gene_length)
                elif strand == '-':
                    #nearby_summit_idx = (dist_to_feature_array < dist5p) & (dist_to_feature_array > (-1 * dist3p))
                    nearby_summit_idx = (dist_to_feature_array < dist5p) & (dist_to_feature_array > (-1 * dist3p)) & (dist_to_feature_array > (-1 * gene_length))
                if nearby_summit_idx.any():
                    nearby_summit_poss = summit_array[nearby_summit_idx]
                    for summit_pos in nearby_summit_poss:
                        included_summit = chrom + '\t' + str(summit_pos)
                        summits_included.append(included_summit)

                        summit_values = np.array([i[2] for i in chipDataSamplesContainer[chip_data_sample_name].summits if i[0] == chrom])
                        self[gene]['in_set'][chip_data_sample_name] = max(summit_values[nearby_summit_idx])
                    self[gene]['in_set']['none'] = False

            summits_included = list(set(summits_included))
            # for i in chipDataSamples[chipDataSample].summit_tss_dists:
            #    if i.replace(':', '\t') not in summits_included:
            #        print i, chipDataSamples[chipDataSample].summit_tss_dists[i], "not in summits included!"
            #    # else:
            #        #print i, "is in summits included!"

            # for i in summits_included: print chipDataSample, i
            print "%s summits from sample %s were included in the plot" % (len(summits_included), chip_data_sample_name)

def loadReferenceData(pugh_reference_genes_fn, yeast_genes_reference_fn='../reference/yeast_genes.tab', restrict_to_verified=False):
    '''
    Reads the TSS and nucleosome data from reference data file (e.g. Pugh 2009 file)
    and loads returns it in a NucSampleDict.
    '''
    try:
        ref_fh = open(pugh_reference_genes_fn, 'r')
    except:
        print "Error opening the reference gene annotation file:  %s..." % (pugh_reference_genes_fn)
        return

    # Irritating hash table to replace chr01-style chrom names
    chrom_lookup_dict = { 'chr01' : 'chrI',
                    'chr02' : 'chrII',
                    'chr03' : 'chrIII',
                    'chr04' : 'chrIV',
                    'chr05' : 'chrV',
                    'chr06' : 'chrVI',
                    'chr07' : 'chrVII',
                    'chr08' : 'chrVIII',
                    'chr09' : 'chrIX',
                    'chr10' : 'chrX',
                    'chr11' : 'chrXI',
                    'chr12' : 'chrXII',
                    'chr13' : 'chrXIII',
                    'chr14' : 'chrXIV',
                    'chr15' : 'chrXV',
                    'chr16' : 'chrXVI',
                    }
    # get the header for the dictionary keys
    header = ref_fh.readline().strip().split('\t')

    # get the subset of genes we're interested in
    verified_gene_list = []

    with open(yeast_genes_reference_fn) as ref_gene_fh:
        for line in ref_gene_fh:
            if line.startswith('ORF'):
                continue
            line = line.strip().split('\t')
            verified_gene_list.append(line[0])

    genes = MnaseSampleDict(pugh_reference_genes_fn)
    genes.restrict_to_verified = restrict_to_verified

    while True:
        line = ref_fh.readline()
        if line:
            line = line.strip().split('\t')
        else:
            break
        gene = line[0]
        #if (gene not in verified_gene_list) and (restrict_to_verified_gene_list == True) and gene.startswith('Y'):
        if (gene not in verified_gene_list) and (restrict_to_verified == True):
            continue
        genes[gene] = {}
        # unusual flat file read b/c rows are not all the same length
        for i in range(len(line)):
            key = header[i]
            # irritating chromNumber replacement
            if key == 'chrom':
                if line[i] in chrom_lookup_dict:
                    line[i] = chrom_lookup_dict[line[i]]
            if key == 'TSS' or key == 'TTS':
                genes[gene][key] = int(line[i])
            else:
                genes[gene][key] = str(line[i])

    return genes

class ChipDataSamplesContainer(dict):
    '''
    Container for ChIP sample datasets.
    '''
    def __init__(self, top_sample_dir, chip_sample_name_list, bed_fn_suffix='_subpeaks.summits.bed'):
        '''
        - top_sample_dir is the top-level sample dir for the all ChIP experiments

        - chip_sample_name_list is a list of the chip samples to be loaded
              (we do not always  want them all).
        '''
        self.top_sample_dir = top_sample_dir
        sample_dir_list = glob(self.top_sample_dir + '/Sample*/')
        for sample_dir in sample_dir_list:
            sample_metadata_dict = parseSampleMetaData(sample_dir)
            if sample_metadata_dict['sample_name'] in chip_sample_name_list:
                self[ sample_metadata_dict['sample_name'] ] = ChipDataSample(sample_dir, bed_fn_suffix)

    def _checkParams():
        assert os.path.isdir(top_sample_dir)

    def loadAllSummits(self):
        '''
        Loads all the summits for the chipSamples in the container.
        '''
        for chipDataSample in self:
            self[chipDataSample].loadSummits()

    def calcSummitVennStats(self):
        '''
        Calculate enrichment statistics between each set of ChIP summits to
        determine the overlap between the bound regions, as opposed to the
        overlap between the loci bound.
        '''
        pass

    def plotSummitToTSSDists(self, mnaseDataSamplesContainer, dist_to_feature3p, dist_to_feature5p, out_fn=''):
        '''
        '''
        # gets the list of genes from the first of the mnaseDataSamples
        genes = dict(mnaseDataSamplesContainer.values()[0].items())

        for chip_data_sample_name in self:
            good_count = 0
            bad_count = 0
            summit_tss_dists = {}
            for summit in self[chip_data_sample_name].summits:
                summit_name = summit[0] + ':' + str(summit[1])
                chrom = summit[0]
                if chrom == 'chrM':
                    continue
                summit_pos = summit[1]
                tss_pos_array = np.array([genes[gene]['TSS'] for gene in genes if genes[gene]['chrom'] == chrom])

                diff_array = tss_pos_array - summit_pos
                # min_summit_tss_dist = min(diff_array)
                min_abs_summit_tss_dist = min(abs(diff_array))
                # if the min distance is positive, we take it with a bias
                # I dunno a better heuristic that doesn't involve random
                # choice
                # if the dist is in the diff_array, then it's downstream of the TSS,
                # so as long as it's less that the dist3p and the gene_length, we count it
                if min_abs_summit_tss_dist in diff_array:
                    summit_tss_dists[summit_name] = min_abs_summit_tss_dist
                    # if the dist is less that the the dist3p, then
                    if (min_abs_summit_tss_dist < dist_to_feature3p):
                        #print sample, "good:", summit_name
                        good_count += 1
                    else:
                        #print sample, 'bad', summit_name
                        bad_count += 1
                # if the dist is not in the diff_array, it's upstream of the TSS,
                # so, as long as it's greater than the -1*dist5p, then it's
                # we can match to the -1*diff_array
                elif min_abs_summit_tss_dist in diff_array * -1:
                    summit_tss_dists[summit_name] = -1 * min_abs_summit_tss_dist
                    if (-1*min_abs_summit_tss_dist) > (-1*dist_to_feature5p):
                        #print sample, "good:", summit_name
                        good_count += 1
                    else:
                        #print sample, 'bad', summit_name
                        bad_count += 1
            print chip_data_sample_name, "good count is", good_count
            print chip_data_sample_name, "bad count is", bad_count
            self[chip_data_sample_name].summit_tss_dists = summit_tss_dists

            plt.figure()
            plt.hist(summit_tss_dists.values(), bins=len(summit_tss_dists) - 1)
            plt.xlim((-1000,1000))
            out_fn = 'output/summit_dist_to_tss_histogram_' + chip_data_sample_name + '_' + str(dist_to_feature3p) + '_' + str(dist_to_feature5p) + '.tiff'
            plt.savefig(out_fn)

class ChipDataSample(dict):
    '''
    ChIP data object.
    '''
    def __init__(self, sample_dir, bed_fn_suffix):
        self.summits = []
        self.peaks = {}
        self.regions = []
        self.sample_dir = sample_dir
        self.bed_fn_suffix = bed_fn_suffix
        self.peaks_bed_fn = ''
        self.bedGraph_fn = glob(self.sample_dir + '/FASTQ_single_sacCer3/bwa*/*.bedGraph')[0]
        self.loadMetadata()
        self.chromosome_size_fn = 'reference/Saccharomyces_cerevisiae/UCSC/sacCer3/Annotation/Genes/ChromInfo.txt'

    def loadMetadata(self):
        '''
        Loads metadata from SampleID file to object attributes.
        '''
        sample_metadata_dict = parseSampleMetaData(self.sample_dir)
        for key, value in sample_metadata_dict.items():
            setattr(self, key, value)

    def loadPeaks(self, peaks_bed_fn_suffix):
        print "Loading ChIP Peaks for sample %s" % (self.sample_name)

        tmp = glob(self.sample_dir + '/FASTQ*/bwa/macs/*' + peaks_bed_fn_suffix)
        if len(tmp) != 1:
            print "Error: more than one bed file with this suffix %s" % (peaks_bed_fn_suffix)
            return
        
        peaks_bed_fh = open(tmp[0])
        for line in peaks_bed_fh:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            chrom = line[0]
            start = int(line[1])
            stop = int(line[2])
            name = line[3]
            score = np.float16(line[4])

            self.peaks[name] = {'chrom' : chrom, 
                                'start' : start, 
                                'stop' : stop, 
                                'score' : score,
                                }

    def loadSummits(self):
        '''
        Loads summits data into a list.
        '''
        print "Loaded ChIP summits data for sample %s..." % (self.sample_name)
        tmp = glob(self.sample_dir + '/FASTQ*/bwa/macs/*' + self.bed_fn_suffix)
        #tmp = glob(self.sample_dir + '/macs/*' + self.bed_fn_suffix)

        if tmp == []:
            print "\tNo summits file for ChIP Sample %s... none loaded." % (self.sample_name)
            return
        else:
            summit_fn = tmp[0]

        try:
            summit_fh = open(summit_fn)
        except:
            print "Error opening summit file %s: " % (summit_fn)

        for line in summit_fh:
            if line.startswith('Chromosome'):
                continue
            line = line.strip().split('\t')
            chrom = line[0]
            if self.bed_fn_suffix == '_summits.bed':
                summit = int(line[1])
                score = float(line[4])
            elif self.bed_fn_suffix == '_subpeaks.summits.bed':
                summit = int(line[1])
                score = float(line[4])
            else:
                import pdb
                pdb.set_trace()
            self.summits.append([chrom, summit, score])
        print "\tdone."


    def loadChromScoreVectorDict(self):
        '''
        Accepts a bedGraph filename and loads all chromosomes into float vectors
        of the scores across the bedGraph regions.
        '''

        def _setChromosomeSizes():
            '''
            '''
            self.chromosome_size_dict = {}
            chromosome_size_fh = open(self.chromosome_size_fn)
            for line in chromosome_size_fh:
                line = line.strip().split('\t')
                if line == ['']:
                    continue
                else:
                    print "...%s, %s" % (line[0], line[1]) 
                    self.chromosome_size_dict[line[0]] = int(line[1])

            print "...done getting chromosome sizes."

        try:
            bg_fh = open(self.bedGraph_fn, 'r')
        except:
            print "Error opening bedGraph file %s." % (self.bedGraph_fn)

        print "Loading bedGraph scores from %s..." % (self.bedGraph_fn)

        _setChromosomeSizes()
        chrom_score_vector_dict = {}

        for chrom in self.chromosome_size_dict:
            chrom_score_vector_dict[chrom] = np.zeros(self.chromosome_size_dict[chrom], dtype=float)

        for line in bg_fh:
            if line.startswith('track'):
                continue
            line = line.strip('\n').split('\t')
            chrom = line[0]
            start = int(line[1]) - 1  # bed format is 1-indexed
            stop = int(line[2]) - 1
            score = float(line[3])
            if chrom in chrom_score_vector_dict:
                chrom_score_vector_dict[chrom][start:stop] = score

        print "\t...done."
        self.chrom_score_vector_dict = chrom_score_vector_dict

def loadAnnotatedNucsData(annotated_nuc_pos_list, sample_name_list, mnaseDataSamplesContainer=None,
                          pugh_reference_genes_fn='reference/justGenesPugh.tab', annotation_model=5,
                          mnase_data_top_dir='../data/mnaseData/', restrict_to_verified=False):
    '''
    Accepts:
         mnaseDataSamplesContainer - A NucDataSamplesContainer  object to load nucData into;
                                if this is unspecified, the function instantiates an empty
                                container to load
         genesfn - the gene annotation file (from Pugh 2009)
         annotated_nuc_pos_list - a list of nucleosome positions to annotate
         sample_name_list - a list of samples (e.g. experiments) to annotate

    Returns:
         sample - the sample data dictionary keyed by genes from the Pugh 2009 file
         nucsData - a dictionary of the samples with nuc position
                    annotations by gene
    '''
    ###
    # Load the sample data for the datasets requested.
    ###

    def _loadSampleMetadata(sample_dir):
        '''
        Load sample metadata into the mnaseDataSample.
        '''
        def _loadBedGraph():
            bedGraph_fn_list = glob(sample_dir + '/FASTQ*/bwa/*scaled.bedGraph')
            #bedGraph_fn_list = glob(sample_dir + '/bwa*/*.bedGraph')
            assert len(bedGraph_fn_list) == 1
            mnaseDataSample.bedGraph_fn = bedGraph_fn_list[0]

        mnaseDataSample = loadReferenceData(pugh_reference_genes_fn, restrict_to_verified=restrict_to_verified)
        mnaseDataSample.loadMetadata(sample_dir)
        if mnaseDataSample.sample_name in sample_name_list:
            _loadBedGraph()
            mnaseDataSample.annotation_model = annotation_model
            mnaseDataSamplesContainer.annotateSample(mnaseDataSample)
            # need to annotate the fake ChIP binding for the unmasked data vectors
            mnaseDataSample.annotateNoneChipBinding()
            # calculates the annotated nucleosome distances from the TSSes
            mnaseDataSample.makeNucsDists(annotated_nuc_pos_list)

    #####
    ## Main body starts here
    #####

    ### why this none check? can this function update an existing object?
    if mnaseDataSamplesContainer == None:
        mnaseDataSamplesContainer = MnaseDataSamplesContainer()

    mnaseDataSamplesContainer.annotated_nuc_pos_list = annotated_nuc_pos_list

    mnase_data_dir_list = glob(mnase_data_top_dir + '/Sample*')
    if mnase_data_dir_list == []:
        print "Error finding MNase data in %s."  % ( mnase_data_top_dir )
        print "Exiting."
        sys.exit()
    else:
        print "Loading MNase-Seq data from these directories:"
        for mnase_data_dir in mnase_data_dir_list:
            print "\t%s" % (mnase_data_dir)
            _loadSampleMetadata(mnase_data_dir)

    return mnaseDataSamplesContainer

def main():
    '''
    '''
    print "module not intended to be run directly"
    return

if __name__ == '__main__':

    main()
