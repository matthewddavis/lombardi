import sys
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from glob import glob
from dataObjects import *
from analysisTools import *
from scipy.stats import hypergeom, ttest_1samp, ks_2samp, ttest_ind, scoreatpercentile, f_oneway, mannwhitneyu

class MetaFeatureGraph(object):
    '''
    This whole nightmarish object class should have been merged into the graphingLib file,
    but it rests here in its completed hackness.
    '''
    def __init__(self, mnaseDataSamplesContainer, density_type, reference_dataset, chipDataSamplesContainer=['none'],
                 feature='TSS', dist5p=500, dist3p=500, set_fig2D_colors_flag=False, fig_fn='MetaFeature.tiff'):
        '''
        reference_dataset can be one of the chip experiment names, or 'none'; or an fpkm cutoff gorup
        '''
        self.mnaseDataSamplesContainer = mnaseDataSamplesContainer
        self.chipDataSamplesContainer = chipDataSamplesContainer
        self.reference_dataset = reference_dataset
        self.density_type = density_type
        self.score_vectors_list = []
        self.feature = feature
        self.dist5p = dist5p
        self.dist3p = dist3p
        self.fig_fn = fig_fn
        self.stats_fn = fig_fn.split('.')[0] + '.stats'
        self.set_fig2D_colors_flag = set_fig2D_colors_flag

    class ScoreVectors(list):
        '''
        '''
        def __init__(self, density_type, sample_name, reference_dataset, vector_limit_type):
            '''
            '''
            self.vector_limit_type = vector_limit_type
            self.density_type = density_type
            self.sample_name = sample_name
            self.reference_dataset = reference_dataset
            self.plot_color = (0,0,0)

    def writeBoundTSSTable(self):
        '''
        Write out a table of the ChIP samples and the genes with a bound
        TSS for each sample in the figure.
        '''
        out_fh = open(self.fig_fn + '.boundTSS.tab', 'w')

        histone_gene_ref_fh = open('reference/histone_genes.tab')
        histone_genes = {}
        for line in histone_gene_ref_fh:
            line = line.strip().split('\t')
            histone_genes[line[1]] = line[0]

        # gets the list of genes from the first of the mnaseDataSamples
        genes = dict(self.mnaseDataSamplesContainer.values()[0].items())

        for chip_sample in self.chipDataSamplesContainer:
            out_fh.write('%s\n' % (chip_sample))

            for gene in genes:
                if genes[gene]['in_set'][chip_sample] != False:
                    if gene in histone_genes:
                        out_fh.write("\t%s\t%s\n" % (gene, histone_genes[gene]))
                    else:
                        out_fh.write("\t%s\n" % (gene))
        out_fh.close()

    def calculateBoundScoreVectorsMNase(self, mnase_sample_name):
        '''
        Calculates the mean read density over the area spanned on each side of the locus feature (e.g. TSS)
        from -dist5p to +dist3p.

        If unbound is True, modifies vectors for both the bound and unbound regions in the specified ChIP sample,
        otherwise, just the bound vector is returned.

        density_type can either be "chip" or "mnase" and the read densities that correspond are returned in the vectors.

        cases can be "both" "bound" or "unbound"
        '''
        span = self.dist3p - (-1 * self.dist5p)

        self.mnaseDataSamplesContainer[mnase_sample_name].chromScoreVectorDict = \
            loadChromScoreVectorDict(self.mnaseDataSamplesContainer[mnase_sample_name].bedGraph_fn)

        if 'fpkm' in self.reference_dataset:
            posScoreVectors = self.ScoreVectors( self.density_type, self.mnaseDataSamplesContainer[mnase_sample_name].sample_name,
                                                 self.reference_dataset, vector_limit_type='exp_top')
            negScoreVectors = self.ScoreVectors( self.density_type, self.mnaseDataSamplesContainer[mnase_sample_name].sample_name,
                                                 self.reference_dataset, vector_limit_type='exp_bottom')

        elif 'IP' in self.reference_dataset:
            posScoreVectors = self.ScoreVectors( self.density_type, self.mnaseDataSamplesContainer[mnase_sample_name].sample_name,
                                                 self.reference_dataset, vector_limit_type='chip_pos')
            negScoreVectors = self.ScoreVectors( self.density_type, self.mnaseDataSamplesContainer[mnase_sample_name].sample_name,
                                                 self.reference_dataset, vector_limit_type='chip_neg')
        elif self.reference_dataset == 'none':
            scoreVectors = self.ScoreVectors(self.density_type, self.mnaseDataSamplesContainer[mnase_sample_name].sample_name,
                                             self.reference_dataset, vector_limit_type='raw_data')

        for gene in self.mnaseDataSamplesContainer[mnase_sample_name]:
            chrom = self.mnaseDataSamplesContainer[mnase_sample_name][gene]['chrom']
            strand = self.mnaseDataSamplesContainer[mnase_sample_name][gene]['Strand']
            tssPos = self.mnaseDataSamplesContainer[mnase_sample_name][gene][self.feature]
            if strand == '+':
                start = tssPos - self.dist5p
                stop = tssPos + self.dist3p
                # watch out for running off the chromosome
                scoreVector = np.zeros(span, dtype=float)
                vectorLen = len(self.mnaseDataSamplesContainer[mnase_sample_name].chromScoreVectorDict[chrom][start:stop])
                scoreVector[0:vectorLen] = self.mnaseDataSamplesContainer[mnase_sample_name].chromScoreVectorDict[chrom][start:stop]
            elif strand == '-':
                start = tssPos + self.dist5p
                stop = tssPos - self.dist3p
                # watch out for running off the chromosome
                scoreVector = np.zeros(span, dtype=float)
                vectorLen = len(self.mnaseDataSamplesContainer[mnase_sample_name].chromScoreVectorDict[chrom][stop:start][::-1])
                scoreVector[0:vectorLen] = self.mnaseDataSamplesContainer[mnase_sample_name].chromScoreVectorDict[chrom][stop:start][::-1]
            else:
                print 'Error: strand "%s" invalid for gene %s' % (strand, gene)

            if 'fpkm' in self.reference_dataset:
                if self.mnaseDataSamplesContainer[mnase_sample_name][gene]['in_set']['fpkm_top'] == True:
                    posScoreVectors.append(scoreVector)
                elif self.mnaseDataSamplesContainer[mnase_sample_name][gene]['in_set']['fpkm_bottom'] == True:
                    negScoreVectors.append(scoreVector)
                continue


            if self.reference_dataset == 'none':
                scoreVectors.append(scoreVector)
            elif self.mnaseDataSamplesContainer[mnase_sample_name][gene]['in_set'][self.reference_dataset]:
                posScoreVectors.append(scoreVector)
            else:
                negScoreVectors.append(scoreVector)

        if self.reference_dataset == 'none':
            self.score_vectors_list.append(scoreVectors)
        else:
            self.score_vectors_list.append(posScoreVectors)
            self.score_vectors_list.append(negScoreVectors)

    def calculateBoundScoreVectorsChIP(self, chip_sample_name):
        '''
        Calculates the mean read density over the area spanned on each side of the locus feature (e.g. TSS)
        from -dist5p to +dist3p.

        If unbound is True, modifies vectors for both the bound and unbound regions in the specified ChIP sample,
        otherwise, just the bound vector is returned.
        '''
        span = self.dist3p - (-1 * self.dist5p)

        if 'fpkm' in self.reference_dataset:
            posScoreVectors = self.ScoreVectors( self.density_type, self.chipDataSamplesContainer[chip_sample_name].sample_name,
                                                 self.reference_dataset, vector_limit_type='exp_top')
            negScoreVectors = self.ScoreVectors( self.density_type, self.chipDataSamplesContainer[chip_sample_name].sample_name,
                                                 self.reference_dataset, vector_limit_type='exp_bottom')
        elif ('IP' in self.reference_dataset) or (self.reference_dataset == 'none'):
            posScoreVectors = self.ScoreVectors( self.density_type, self.chipDataSamplesContainer[chip_sample_name].sample_name,
                                                 self.reference_dataset, vector_limit_type='chip_pos')
            negScoreVectors = self.ScoreVectors( self.density_type, self.chipDataSamplesContainer[chip_sample_name].sample_name,
                                                 self.reference_dataset, vector_limit_type='chip_neg')

        assert self.density_type == 'chip'
        self.chipDataSamplesContainer[chip_sample_name].loadChromScoreVectorDict()


        # gets the list of genes from the first of the mnaseDataSamples
        genes = dict(self.mnaseDataSamplesContainer.values()[0].items())

        for gene in genes:
            chrom = genes[gene]['chrom']
            strand = genes[gene]['Strand']
            tssPos = genes[gene][self.feature]
            if strand == '+':
                start = tssPos - self.dist5p
                stop = tssPos + self.dist3p
                # watch out for running off the chromosome
                # in practice, nothing is removed from the dataset
                # of verified genes with distances < 1000bp
                if start < 0:
                    continue
                if stop > len(self.chipDataSamplesContainer[chip_sample_name].chrom_score_vector_dict[chrom]):
                    continue
                scoreVector = np.zeros(span, dtype=float)
                scoreVector = self.chipDataSamplesContainer[chip_sample_name].chrom_score_vector_dict[chrom][start:stop]

            #elif strand == '-':
            #    start = tssPos - self.dist3p
            #    stop = tssPos + self.dist5p
            #    # watch out for running off the chromosome
            #    # in practice, nothing is removed from the dataset
            #    # of verified genes with distances < 1000bp
            #    if start < 0:
            #        continue
            #    if stop > len(self.chipDataSamplesContainer[chip_sample_name].chrom_score_vector_dict[chrom]):
            #        continue
            #    scoreVector = np.zeros(span, dtype=float)
            #    scoreVector = self.chipDataSamplesContainer[chip_sample_name].chrom_score_vector_dict[chrom][start:stop]

            elif strand == '-':
                start = tssPos - self.dist3p
                stop = tssPos + self.dist5p
                # watch out for running off the chromosome
                # in practice, nothing is removed from the dataset
                # of verified genes with distances < 1000bp
                if start < 0:
                    continue
                if stop > len(self.chipDataSamplesContainer[chip_sample_name].chrom_score_vector_dict[chrom]):
                    continue

                scoreVector = np.zeros(span, dtype=float)
                vectorLen = len(self.chipDataSamplesContainer[chip_sample_name].chrom_score_vector_dict[chrom][start:stop][::-1])
                if vectorLen != len(scoreVector):
                    pdb.set_trace()
                scoreVector[0:vectorLen] = self.chipDataSamplesContainer[chip_sample_name].chrom_score_vector_dict[chrom][start:stop][::-1]

            else:
                print 'Error: strand "%s" invalid for gene %s' % (strand, gene)

            if 'fpkm' in self.reference_dataset:
                if genes[gene]['in_set']['fpkm_top'] == True:
                    posScoreVectors.append(scoreVector)
                elif genes[gene]['in_set']['fpkm_bottom'] == True:
                    negScoreVectors.append(scoreVector)
                continue
            if genes[gene]['in_set'][chip_sample_name]:
                if chip_sample_name != 'none':
                    if (self.reference_dataset == 'none') or (genes[gene]['in_set'][self.reference_dataset]):
                        posScoreVectors.append(scoreVector)
                    else:
                        negScoreVectors.append(scoreVector)
                else:
                    negScoreVectors.append(scoreVector)
            else:
                negScoreVectors.append(scoreVector)

        #if chip_sample_name != 'none':
        self.score_vectors_list.append(negScoreVectors)
        #else:
        self.score_vectors_list.append(posScoreVectors)

    def graph(self, cases='both', showControl=False):
        '''
        '''
        # load chip ScoreVectors if required
        if self.density_type == 'chip':
            for chip_data_sample_name in self.chipDataSamplesContainer:
                self.calculateBoundScoreVectorsChIP(chip_data_sample_name)
        if self.density_type == 'mnase':
            for mnase_data_sample_name in self.mnaseDataSamplesContainer:
                self.calculateBoundScoreVectorsMNase(mnase_data_sample_name)

        def _setSampleColors():
            sample_color_dict = {}

            #if cases == 'neg':
            if (self.density_type == 'mnase') and (cases != 'both') and ('fpkm' not in self.reference_dataset):
                colors = getColors(len(self.score_vectors_list))
            #elif (self.density_type == 'chip') and cases=='both':
            elif cases == 'both':
                colors = getColors(len(self.score_vectors_list)/2)
            #elif cases == 'both' and ('fpkm' in self.reference_dataset):
            #    colors = getColors(len(self.score_vectors_list)/2)
            for scoreVector in sorted(self.score_vectors_list, key=lambda x: x.sample_name):
                if scoreVector.sample_name not in sample_color_dict:
                    sample_color_dict[scoreVector.sample_name] = colors.pop()
                    scoreVector.plot_color = sample_color_dict[scoreVector.sample_name]
                elif scoreVector.sample_name in sample_color_dict:
                    scoreVector.plot_color = sample_color_dict[scoreVector.sample_name]

        def _setFigure2Colors():
            colors = getColors(len(self.score_vectors_list))
            for scoreVector in self.score_vectors_list:
                if 'OE' in scoreVector.sample_name:
                    scoreVector.plot_color = (0.0, 0.0, 1.0)
                elif 'GAL_YTA7' in scoreVector.sample_name:
                    scoreVector.plot_color = (0.0, 1.0, 0.0)
                if 'DEX_YTA7' in scoreVector.sample_name:
                    scoreVector.plot_color = (1.0, 0.0, 0.0)

        if self.set_fig2D_colors_flag:
            _setFigure2Colors()
        else:
            _setSampleColors()

        fig = plt.figure(dpi=300, figsize=(12.8, 10.24))
        ax = plt.subplot(111)

        x = np.arange( (-1 * self.dist5p), self.dist3p)

        statsfh =  open(self.stats_fn, 'w')
        statsfh.write('SampleID\tRestrictionCriteria\t#Loci\n')

        for scoreVector in sorted(self.score_vectors_list, key=lambda x: x.reference_dataset):
            if 'neg' in scoreVector.vector_limit_type:
                statsfh.write('%s\t%s\t%s\n' % (scoreVector.sample_name.decode('unicode-escape').encode('UTF-8'), 'none', len(scoreVector)))
            else:
                statsfh.write('%s\t%s\t%s\n' % (scoreVector.sample_name.decode('unicode-escape').encode('UTF-8'), scoreVector.vector_limit_type, len(scoreVector)))
            y = np.mean(scoreVector, axis=0)
            if self.density_type == 'mnase':
                labelText = self.mnaseDataSamplesContainer[scoreVector.sample_name].common_id.decode('unicode-escape')
            elif self.density_type == 'chip':
                if scoreVector.sample_name == 'none':
                    labelText = ''
                else:
                    labelText = self.chipDataSamplesContainer[scoreVector.sample_name].common_id.decode('unicode-escape')
                    #labelText = self.chipDataSamplesContainer[scoreVector.sample_name].sample_name
            if cases == 'both':
                # Top 20% FPKM
                #if (scoreVector.vector_limit_type == 'exp_top') and (self.reference_dataset == 'fpkm_top'):
                if (scoreVector.vector_limit_type == 'exp_top'):
                    labelText = labelText + ' Top 20% FPKM'
                    lineStyle = '-'
                # elif (scoreVector.vector_limit_type == 'exp_bottom') and (self.reference_dataset == 'fpkm_top'):
                elif (scoreVector.vector_limit_type == 'exp_bottom'):
                    labelText = labelText + ' Bottom 20% FPKM'
                    lineStyle = '--'
                #
                # Bottom 20% FPKM
                #elif (scoreVector.vector_limit_type == 'exp_pos') and (self.reference_dataset == 'fpkm_bottom'):
                #    labelText = labelText + ' Bottom 20% FPKM'
                #    lineStyle = '-'
                #elif (scoreVector.vector_limit_type == 'exp_neg') and (self.reference_dataset == 'fpkm_bottom'):
                #    labelText = labelText + ' Top 80% FPKM'
                #    lineStyle = '--'
                # bound in a ChIP Experiment
                elif (scoreVector.vector_limit_type == 'chip_pos'):
                    labelText = labelText + ' Bound'
                    lineStyle = '-'
                elif (scoreVector.vector_limit_type == 'chip_neg'):
                    labelText = labelText + ' Unbound'
                    lineStyle = '--'
                else:
                    continue
                ax.plot(x, y, linewidth=3, alpha=0.5, c=scoreVector.plot_color, label=labelText, linestyle=lineStyle)
            # Just the positive vectors
            elif cases == 'pos':
                # Top 20% FPKM
                if (scoreVector.vector_limit_type == 'exp_pos') and (self.reference_dataset == 'fpkm_top'):
                    labelText = labelText + ' Top 20% FPKM'
                    lineStyle = '-'
                # Bottom 20% FPKM
                elif (scoreVector.vector_limit_type == 'exp_pos') and (self.reference_dataset == 'fpkm_bottom'):
                    labelText = labelText + ' Bottom 20% FPKM'
                    lineStyle = '-'
                # bound in a ChIP Experiment
                elif (scoreVector.vector_limit_type == 'chip_pos'):
                    labelText = labelText + ' Bound'
                    lineStyle = '-'
                else:
                    continue
                ax.plot(x, y, linewidth=3, alpha=0.5, c=scoreVector.plot_color, label=labelText, linestyle=lineStyle)
            # Just the negative vectors -- used for plotting raw MNase data
            elif cases == 'neg':
                # raw MNase
                if (scoreVector.vector_limit_type == 'raw_data'):
                    labelText = labelText + ''
                    lineStyle = '-'
                else:
                    continue
                raw_fn = self.stats_fn.replace('.stats', '.raw.') + scoreVector.sample_name
                raw_fh = open(raw_fn, 'w')
                for i,datum in enumerate(x):
                    raw_fh.write('%s\t%s\n' % (x[i], y[i]) )
                ax.plot(x, y, linewidth=3, alpha=0.5, c=scoreVector.plot_color, label=labelText, linestyle=lineStyle)
            '''
            elif cases == 'bound':
                if scoreVector.bound:
                    labelText = labelText + ' Bound'
                    lineStyle = '-'
                    ax.plot(x, y, linewidth=3, alpha=0.5, c=scoreVector.plot_color, label=labelText, linestyle=lineStyle)
            elif cases == 'unbound':
                if scoreVector.unbound:
                    lineStyle = '-'
                    ax.plot(x, y, linewidth=3, alpha=0.5, c=scoreVector.plot_color, label=labelText, linestyle=lineStyle)
            elif cases == 'top20fpkm':
                if scoreVector.bound:
                    labelText = labelText + ' > Top 20% FPKM'
                    lineStyle = '-'
                elif scoreVector.unbound:
                    labelText = labelText + ' < Top 20% FPKM'
                    lineStyle = '--'
                ax.plot(x, y, linewidth=3, alpha=0.5, c=scoreVector.plot_color, label=labelText, linestyle=lineStyle)
            elif cases == 'bottom20fpkm':
                if scoreVector.bound:
                    labelText = labelText + ' Bottom 20% FPKM'
                    lineStyle = '-'
                elif scoreVector.unbound:
                    labelText = labelText + ' Other FPKM'
                    lineStyle = '--'
                ax.plot(x, y, linewidth=3, alpha=0.5, c=scoreVector.plot_color, label=labelText, linestyle=lineStyle)
                '''
        plt.xlim((-1*self.dist5p, self.dist3p))

        box = ax.get_position()
        #ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])

        fontP = FontProperties()
        fontP.set_size('medium')
        #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, prop=fontP)
        lines, labels = ax.get_legend_handles_labels()
        label_dict = dict(zip(labels, lines))
        sorted_keys = sorted(label_dict.keys())

        if self.set_fig2D_colors_flag:

            labels = [u'Gal Yta7 OE IP Top 20% FPKM',
                      u'Gal Yta7 OE IP Bottom 20% FPKM',
                      u'Gal Yta7 IP Top 20% FPKM',
                      u'Gal Yta7 IP Bottom 20% FPKM',
                      u'Glu Yta7 IP Top 20% FPKM',
                      u'Glu Yta7 IP Bottom 20% FPKM',
                      ]

            lines = [label_dict[k] for k in labels]
            ax.legend(lines, labels, loc='best', frameon=False, prop=fontP)
        else:
            ax.legend(loc='best', frameon=False, prop=fontP)
        #ax.legend( [label_dict[k] for k in sorted_keys], sorted_keys, loc='best', frameon=False, prop=fontP)
        #ax.legend(loc='best', bbox_to_anchor=(1, 0.5), prop=fontP)
        #ax.legend(loc='upper left', prop=fontP)
        #ax.legend(loc=0)
        #plt.title(chipDataSample.sample_name)
        if self.feature == 'TSS':
            plt.xlabel('Distance from Transcription Start Site (TSS)')
        #elif:
        #    pass
        #    #plt.xlabel('Distance from ' + self.feature)

        plt.tick_params(pad=10, labeltop='off', labelright='off')

        if self.density_type == 'chip':
            plt.ylabel('Scaled ChIP Read Density')
        elif self.density_type == 'mnase':
            plt.ylabel('Scaled MNase Read Count')
        if not self.fig_fn.endswith('.tiff'):
            self.fig_fn += '.tiff'
        fig.savefig(self.fig_fn)

def calcChIPVennDiagramStats(mnase_sample_name_list, annotated_nuc_pos_list,
                             chip_sample_name_list=['none'], reference_dataset='none',
                             chip_top_sample_dir='../data/chipData', feature='TSS',
                             dist_to_feature5p=450, dist_to_feature3p=1000,
                             out_fn='output/VennChipStats.txt',
                             restrict_to_verified=False):
    '''
    '''
    mnaseDataSamplesContainer = loadAnnotatedNucsData(annotated_nuc_pos_list=annotated_nuc_pos_list,
                                                      sample_name_list=mnase_sample_name_list,
                                                      annotation_model=5,
                                                      restrict_to_verified=restrict_to_verified)

    chipDataSamplesContainer = ChipDataSamplesContainer(top_sample_dir=chip_top_sample_dir,
                                                        chip_sample_name_list=chip_sample_name_list,
                                                        #bed_fn_suffix='_summits.bed',
                                                        bed_fn_suffix='_subpeaks.summits.bed',
                                                        )
    chipDataSamplesContainer.loadAllSummits()
    mnaseDataSamplesContainer.annotateChipBindingAllSamples(chipDataSamplesContainer, feature=feature, dist5p=dist_to_feature5p, dist3p=dist_to_feature3p)

    out_fh = open(out_fn, 'w')
    mnase_sample_name = mnase_sample_name_list[0]
    ref_chip_sample_name = reference_dataset
    for test_chip_sample_name in chip_sample_name_list:
        ref_vector = np.array([mnaseDataSamplesContainer[mnase_sample_name][gene]['in_set'][ref_chip_sample_name]
                               for gene in mnaseDataSamplesContainer[mnase_sample_name]], dtype=bool)
        test_vector = np.array([mnaseDataSamplesContainer[mnase_sample_name][gene]['in_set'][test_chip_sample_name]
                                for gene in mnaseDataSamplesContainer[mnase_sample_name]], dtype=bool)
        out_fh.write("ref sample is %s\n" % (ref_chip_sample_name))
        out_fh.write("test sample is %s\n" % (test_chip_sample_name))

        typeIobjs = sum(ref_vector)
        draws = sum(test_vector)
        hits = sum(ref_vector & test_vector)
        popSize = len(ref_vector)
        out_fh.write("draws / hits / targets is %s / %s / %s\n" % (draws, hits, typeIobjs))
        rv = hypergeom(popSize, typeIobjs, draws)
        p_val = rv.pmf(hits)
        out_fh.write("p_val is %0.e\n" % (p_val))

def graphMetaFeature(mnase_sample_name_list, annotated_nuc_pos_list, density_type,
                     chip_sample_name_list=[], reference_dataset='none',
                     chip_top_data_dir='../data/chipData/', cases='both',
                     feature='TSS', dist5p=500, dist3p=700, dist_to_feature5p=450,
                     dist_to_feature3p=1000,
                     annotate_expression=False,
                     #bed_fn_suffix='_summits.bed',
                     bed_fn_suffix='_subpeaks.summits.bed',
                     set_fig2D_colors_flag=False,
                     fig_fn='output/MetaFeature.tiff',
                     restrict_to_verified=False):
    '''
    This function wraps the metaFeatureGraph object and prepares the data sent to it.

    It can be used to plot MNase-Seq data, ChIP-Seq data, or both, considering conditions
    of either type of data in the graph.
    '''
    # Load the annotated nucleosome position data for the specified mnase_sample_name_list
    # Use annotation model 5 to annotate the nuc positions.
    mnaseDataSamplesContainer = loadAnnotatedNucsData(annotated_nuc_pos_list, mnase_sample_name_list, annotation_model=5,
                                                      restrict_to_verified=restrict_to_verified)

    # load the chip data, if it's specified.
    # If chip_sample_name_list is 'none', the data structure is still constructed
    # but without limits on binding data
    chipDataSamplesContainer = ChipDataSamplesContainer(chip_top_data_dir,
                                                        chip_sample_name_list,
                                                        # bedfn_str='/FASTQ_single_sacCer3/bwa/*subpeaks.oldmacs_peaks',
                                                        bed_fn_suffix=bed_fn_suffix,
                                                        )
    chipDataSamplesContainer.loadAllSummits()
    chipDataSamplesContainer.plotSummitToTSSDists(mnaseDataSamplesContainer, dist_to_feature3p, dist_to_feature5p)

    mnaseDataSamplesContainer.annotateChipBindingAllSamples(chipDataSamplesContainer, feature=feature, dist5p=dist_to_feature5p, dist3p=dist_to_feature3p)

    if annotate_expression == True:
        for mnase_data_sample_name in mnaseDataSamplesContainer:
            mnaseDataSamplesContainer[mnase_data_sample_name].annotateFPKM()

    metaFeatureGraph = MetaFeatureGraph(mnaseDataSamplesContainer, density_type, reference_dataset, chipDataSamplesContainer,
                                        feature, dist5p, dist3p, set_fig2D_colors_flag, fig_fn)

    if density_type == 'chip':
        metaFeatureGraph.writeBoundTSSTable()
    metaFeatureGraph.graph(cases=cases)

def makeTSSUnsmoothedMasked(mnaseDataSamplesContainer, chip_sample_name_list=['none'], fig_fn=False, use_common_id_flag=True):
    '''
    Draws figure with unsmoothened nucleosome densities across meta-TSS.
    '''
    print "\tDrawing the unsmoothened nucleosome distributions around the meta-TSS..."

    colors = getColors(len(mnaseDataSamplesContainer))
    fig = plt.figure(dpi=300, figsize=(10.24,7.68))

    out_fn = fig_fn + '.count_table'
    out_fh = open(out_fn, 'w')
    out_fh.write('Sample\tNucleosome Position\tBinding Restrictions\t# Nucleosomes\n')

    ymax = 0
    for chipSample in chip_sample_name_list:
        for mnaseDataSample in mnaseDataSamplesContainer:
            col = colors.pop()
            if use_common_id_flag:
                dataLabel = mnaseDataSamplesContainer[mnaseDataSample].common_id
            else:
                dataLabel = mnaseDataSamplesContainer[mnaseDataSample].sample_name.decode('unicode-escape')
            for annotated_nuc_pos in mnaseDataSamplesContainer.annotated_nuc_pos_list:
                idx = (np.isnan(mnaseDataSamplesContainer[mnaseDataSample].annotated_nuc_dists_dict[annotated_nuc_pos][chipSample]) == False)
                out_fh.write("%s\t%s\t%s\t%s\n" % (mnaseDataSample, annotated_nuc_pos, chipSample, len(idx)))
                if dataLabel:
                    xhist = plt.hist(mnaseDataSamplesContainer[mnaseDataSample].annotated_nuc_dists_dict[annotated_nuc_pos][chipSample][idx],
                                     bins=100, histtype='step', color=col, label=dataLabel, alpha=0.5)
                    dataLabel = False
                else:
                    xhist = plt.hist(mnaseDataSamplesContainer[mnaseDataSample].annotated_nuc_dists_dict[annotated_nuc_pos][chipSample][idx],
                                     bins=100, histtype='step', color=col, alpha=0.5)
                if max(xhist[0]) > ymax:
                    ymax = max(xhist[0])

    plt.ylim((0,ymax*1.1))
    #plt.yticks(())
    plt.ylabel('Scaled Nucleosome Count')
    plt.xlabel('Distance from TSS to Nucleosome Center')
    plt.legend(loc=2, frameon=False)

    plt.tick_params(pad=10, labeltop='off', labelright='off')

    print "\tFigure is made...."
    # The figure is returned either way.
    if fig_fn != False:
        if not fig_fn.endswith('.tiff'):
            fig_fn += '.tiff'
        fig.savefig(fig_fn)
    print "\tFigure is saved..."
    return fig

def makeTSSSmoothedMasked(mnaseDataSamplesContainer, chip_sample_name_list=['none'], smoothing_type='hanning', bin_denominator=10, fig_fn=False, use_common_id_flag=True):
    '''
    Draws figure with smoothened nucleosome densities across meta-TSS.

    You can mask the figure with required binding in one of the chip samples,
    if they have been added to the locus annotations in the mnaseDataSamplesContainer
    structure.
    '''
    out_fn = fig_fn + '.count_table'
    out_fh = open(out_fn, 'w')
    out_fh.write('Sample\tNucleosome Position\tBinding Restrictions\t# Nucleosomes\n')

    print "\tDrawing the smoothened nucleosome distributions around the meta-TSS..."
    colors = getColors(len(mnaseDataSamplesContainer) * len(chip_sample_name_list))
    fig = plt.figure(dpi=300, figsize=(10.24,7.68))
    ymax = 0
    for chip_sample_name in chip_sample_name_list:
        for mnaseDataSample in mnaseDataSamplesContainer:
            col = colors.pop()
            if chip_sample_name == 'none':
                if use_common_id_flag:
                    data_label = mnaseDataSamplesContainer[mnaseDataSample].common_id.decode('unicode-escape')
                else:
                    data_label = mnaseDataSamplesContainer[mnaseDataSample].sample_name.decode('unicode-escape')
            else:
                if use_common_id_flag:
                    data_label = mnaseDataSamplesContainer[mnaseDataSample].common_id.decode('unicode-escape') + ' ' + chip_sample_name
                else:
                    data_label = mnaseDataSamplesContainer[mnaseDataSample].common_id.decode('unicode-escape') + ' ' + chip_sample_name

            for annotated_nuc_pos in mnaseDataSamplesContainer.annotated_nuc_pos_list:
                # the number of nucleosomes being counted in each distribution is
                # reflected here as the length of the following idx
                idx = (np.isnan(mnaseDataSamplesContainer[mnaseDataSample].annotated_nuc_dists_dict[annotated_nuc_pos][chip_sample_name]) == False)
                out_fh.write("%s\t%s\t%s\t%s\n" % (mnaseDataSample, annotated_nuc_pos, chip_sample_name, len(idx)))

                # this if statement is a trick to only plot the data_label once per sample, as we're
                # drawing each annotated nuc distribution separately on the plot for each sample
                if not data_label == '':
                    xhist = plt.hist(mnaseDataSamplesContainer[mnaseDataSample].annotated_nuc_dists_dict[annotated_nuc_pos][chip_sample_name][idx],
                                     bins=len(mnaseDataSamplesContainer[mnaseDataSample].annotated_nuc_dists_dict[annotated_nuc_pos][chip_sample_name][idx]),
                                     histtype='step', color=col, visible=False)
                    windowsize = len(xhist[1]) / bin_denominator
                    x = smooth(xhist[1], windowsize, smoothing_type)[:-1]
                    y = smooth(xhist[0], windowsize, smoothing_type)
                    plt.plot(x,y, label=data_label, linewidth=3, c=col, alpha=0.5)
                    data_label = ''
                else:
                    xhist = plt.hist(mnaseDataSamplesContainer[mnaseDataSample].annotated_nuc_dists_dict[annotated_nuc_pos][chip_sample_name][idx],
                                     bins=len(mnaseDataSamplesContainer[mnaseDataSample].annotated_nuc_dists_dict[annotated_nuc_pos][chip_sample_name][idx]),
                                     histtype='step', color=col, visible=False)
                    windowsize = len(xhist[1]) / bin_denominator
                    x = smooth(xhist[1], windowsize, smoothing_type)[:-1]
                    y = smooth(xhist[0], windowsize, smoothing_type)
                    # plt.plot(x,y, label=data_label, linewidth=3, c=col, alpha=0.5)
                    plt.plot(x,y, linewidth=3, c=col, alpha=0.5)
                if max(y) > ymax:
                    ymax = max(y)

    plt.ylim((0,ymax*1.1))
    #plt.yticks(())
    plt.ylabel('Scaled Nucleosome Count')
    plt.xlabel('Distance from TSS to Nucleosome Center')
    plt.legend(loc=2, frameon=False)

    plt.tick_params(pad=10, labeltop='off', labelright='off')

    print "\tFigure is made...."
    # The figure is returned either way.
    if fig_fn != False:
        if not fig_fn.endswith('.tiff'):
            fig_fn += '.tiff'
        fig.savefig(fig_fn)
    print "\tFigure is saved..."
    return fig


def writeNucDataStatistics(mnaseDataSamplesContainer,
                           ref_nuc_data_sample_name,
                           chip_sample_name='none',
                           test_type='ks',
                           out_fn=False,
                           use_common_id_flag=True,
                           flip_unbound_flag=False):
    '''
    Uses either t-test or K-S test to determine the signicance of the differences
    of the distributions of nucleosomes between samples (paired at each position).

    Writes a simple table out to a file.
    '''
    try:
        out_fh = open(out_fn, 'w')
    except:
        print "Error opening the out file for table statistics %s" % (out_fn)

    out_fh.write("Sample\t ")
    # a little hash to pretty-up the Nuc Pos names
    nuc_pos_name_dict = {'AnN1' : '-1',
                         'AnP1' : '+1',
                         'AnP2' : '+2',
                         'AnP3' : '+3',
                         'AnP4' : '+4',
                         'AnP5' : '+5',
                         'AnP6' : '+6',
                         }
    for nuc_pos in mnaseDataSamplesContainer.annotated_nuc_pos_list:
        out_fh.write("%s\t" % (nuc_pos_name_dict[nuc_pos]))
    out_fh.write('\n')

    for mnaseDataSample in mnaseDataSamplesContainer:
        #if mnaseDataSample == ref_nuc_data_sample_name:
        #    continue
        if use_common_id_flag:
            out_fh.write("%s\t" % (mnaseDataSamplesContainer[mnaseDataSample].common_id.decode('unicode-escape') ).encode('UTF-8'))
        else:
            out_fh.write("%s\t" % (mnaseDataSamplesContainer[mnaseDataSample].sample_name.decode('unicode-escape') ).encode('UTF-8'))

        for nuc_pos in mnaseDataSamplesContainer.annotated_nuc_pos_list:
            if test_type == 'ks':
                alpha_val, p_val, test_mean = calcNucDataStatistics(mnaseDataSamplesContainer, mnaseDataSample,
                                                                    ref_nuc_data_sample_name, nuc_pos,
                                                                    chip_sample_name=chip_sample_name, test_type='ks',
                                                                    flip_unbound_flag=flip_unbound_flag)
                out_fh.write("%3.2f\t" % (round(test_mean)))
            elif test_type == 't':
                t_val, p_val, test_mean = calcNucDataStatistics(mnaseDataSamplesContainer, mnaseDataSample,
                                                                ref_nuc_data_sample_name, nuc_pos,
                                                                chip_sample_name=chip_sample_name, test_type='t',
                                                                #flip_unbound_flag=flip_unbound_flag
                                                                )
                out_fh.write("%3.0f\t" % (round(test_mean)))

            elif test_type == 'f':
                t_val, p_val, test_mean = calcNucDataStatistics(mnaseDataSamplesContainer, mnaseDataSample,
                                                                ref_nuc_data_sample_name, nuc_pos,
                                                                chip_sample_name=chip_sample_name,
                                                                )
                out_fh.write("%s\t" % (test_mean))

        out_fh.write('\n')

        for nuc_pos in mnaseDataSamplesContainer.annotated_nuc_pos_list:
            if test_type == 'ks':
                alpha_val, p_val, test_mean = calcNucDataStatistics(mnaseDataSamplesContainer, mnaseDataSample,
                                                                    ref_nuc_data_sample_name, nuc_pos,
                                                                    chip_sample_name=chip_sample_name, test_type='ks',
                                                                    flip_unbound_flag=flip_unbound_flag)
                out_fh.write("\t(%1.0e)" % (p_val))
            elif test_type == 't':
                t_val, p_val, test_max = calcNucDataStatistics(mnaseDataSamplesContainer, mnaseDataSample,
                                                               ref_nuc_data_sample_name, nuc_pos,
                                                               chip_sample_name=chip_sample_name, test_type='t',
                                                               #flip_unbound_flag=flip_unbound_flag
                                                               )
                out_fh.write("\t(%1.3f)" % (p_val))

            elif test_type == 'f':
                t_val, p_val, test_mean = calcNucDataStatistics(mnaseDataSamplesContainer, mnaseDataSample,
                                                                ref_nuc_data_sample_name, nuc_pos,
                                                                chip_sample_name=chip_sample_name,
                                                                )
                out_fh.write("%s\t" % (test_mean) )

        out_fh.write('\n')
    out_fh.close()

def calcNucDataStatistics(mnaseDataSamplesContainer,
                          test_nuc_data_sample_name,
                          ref_nuc_data_sample_name,
                          nuc_pos,
                          chip_sample_name='none',
                          test_type='t'
                          ):
    '''
    Caluclates the t-test and ks-test values for the test_nuc_data_sample with
    reference to the ref_nuc_data_sample. First a vector of non-NaN entries in each
    dataset are constructed, and then these are used to calculate the test statistics.
    '''

    def _getNucPosVector(nuc_sample_name):
        bound_values = []
        unbound_values = []

        for locus in mnaseDataSamplesContainer[nuc_sample_name]:
            if nuc_pos in mnaseDataSamplesContainer[nuc_sample_name][locus]:
                # bound vector first
                if mnaseDataSamplesContainer[nuc_sample_name][locus]['in_set'][chip_sample_name] != False:
                    bound_values.append(mnaseDataSamplesContainer[nuc_sample_name][locus][nuc_pos]['pos'])
                elif mnaseDataSamplesContainer[nuc_sample_name][locus]['in_set'][chip_sample_name] == False:
                    unbound_values.append(mnaseDataSamplesContainer[nuc_sample_name][locus][nuc_pos]['pos'])

        return (np.array(bound_values), np.array(unbound_values))

    def _removeNaN(_vector):
        # remove NaN
        _vector_idx = np.isnan(_vector) == False
        _vector = _vector[_vector_idx]

        return _vector

    def _initVectorDict(nuc_sample_name):
        bound_vector, unbound_vector = _getNucPosVector(nuc_sample_name)
        bound_vector = _removeNaN(bound_vector)
        unbound_vector = _removeNaN(unbound_vector)

        vector_dict = {'bound_v' : bound_vector,
                       'unbound_v' : unbound_vector,
                       }

        return vector_dict

    test_vector_dict = _initVectorDict(test_nuc_data_sample_name)
    ref_vector_dict = _initVectorDict(ref_nuc_data_sample_name)

    ttest_v = test_vector_dict['bound_v'] - np.mean(ref_vector_dict['bound_v'])
    tref_v = test_vector_dict['unbound_v'] - np.mean(ref_vector_dict['unbound_v'])

    counter = 0
    iter_counter = 0
    test_sample_len = round(len(ttest_v) * .9 )
    ref_sample_len = round(len(tref_v) * .9)

    for i in xrange(1000):
        iter_counter += 1
        if iter_counter % 100 == 0:
            print "iter_counter @ %s" % iter_counter
        np.random.shuffle(ttest_v)
        test_mean = np.mean(ttest_v[0:test_sample_len])
    
        np.random.shuffle(tref_v)
        ref_mean = np.mean(tref_v[0:ref_sample_len])

        if test_mean < ref_mean:
            counter += 1
    p_val = 1 - (float(counter) / 1000)
    t_val = 1.0

    #(t_val, p_val) = ttest_1samp(ttest_v, np.mean(tref_v))
    #(t_val, p_val) = ttest_ind(ttest_v, tref_v)
    #t_val, p_val = mannwhitneyu(ttest_v, tref_v)
    #np.random.shuffle(tref_v)
    #t_val, p_val = mannwhitneyu(ttest_v, tref_v[0:len(ttest_v)])
    test_mean = np.mean(ttest_v)
    #if test_mean >= 0:
    return (t_val, p_val, test_mean)
    #elif test_mean < 0:
    #return (t_val, 1.0 - p_val, test_mean)

def calcNucDataStatistics_not_so_old(mnaseDataSamplesContainer,
                          test_nuc_data_sample_name,
                          ref_nuc_data_sample_name,
                          nuc_pos,
                          chip_sample_name='none',
                          test_type='t'
                          ):
    '''
    Caluclates the t-test and ks-test values for the test_nuc_data_sample with
    reference to the ref_nuc_data_sample. First a vector of non-NaN entries in each
    dataset are constructed, and then these are used to calculate the test statistics.
    '''

    def _getNucPosVector(nuc_sample_name):
        bound_values = []
        unbound_values = []

        for locus in mnaseDataSamplesContainer[nuc_sample_name]:
            if nuc_pos in mnaseDataSamplesContainer[nuc_sample_name][locus]:
                # bound vector first
                if mnaseDataSamplesContainer[nuc_sample_name][locus]['in_set'][chip_sample_name] != False:
                    bound_values.append(mnaseDataSamplesContainer[nuc_sample_name][locus][nuc_pos]['pos'])
                elif mnaseDataSamplesContainer[nuc_sample_name][locus]['in_set'][chip_sample_name] == False:
                    unbound_values.append(mnaseDataSamplesContainer[nuc_sample_name][locus][nuc_pos]['pos'])

        return (np.array(bound_values), np.array(unbound_values))

    def _removeNaN(_vector):
        # remove NaN
        _vector_idx = np.isnan(_vector) == False
        _vector = _vector[_vector_idx]

        return _vector

    def _initVectorDict(nuc_sample_name):
        bound_vector, unbound_vector = _getNucPosVector(nuc_sample_name)
        bound_vector = _removeNaN(bound_vector)
        unbound_vector = _removeNaN(unbound_vector)

        vector_dict = {'bound_v' : bound_vector,
                       'unbound_v' : unbound_vector,
                       }

        return vector_dict

    test_vector_dict = _initVectorDict(test_nuc_data_sample_name)
    ref_vector_dict = _initVectorDict(ref_nuc_data_sample_name)

    #if test_nuc_data_sample_name == 'DEX_Yta7_minus_alpha' and ref_nuc_data_sample_name == 'DEX_WT_alpha':
    #    import pdb
    #    pdb.set_trace()

    ttest_v = test_vector_dict['bound_v'] - np.mean(ref_vector_dict['bound_v'])
    tref_v = test_vector_dict['unbound_v'] - np.mean(ref_vector_dict['unbound_v'])

    #(t_val, p_val) = ttest_1samp(ttest_v, np.mean(tref_v))
    #(t_val, p_val) = ttest_ind(ttest_v, tref_v)
    t_val, p_val = mannwhitneyu(ttest_v, tref_v)
    #np.random.shuffle(tref_v)
    #t_val, p_val = mannwhitneyu(ttest_v, tref_v[0:len(ttest_v)])
    test_mean = np.mean(ttest_v)
    return (t_val, p_val/2.0, test_mean)


def calcNucDataStatistics_ftest(mnaseDataSamplesContainer,
                          test_nuc_data_sample_name,
                          ref_nuc_data_sample_name,
                          nuc_pos,
                          chip_sample_name='none',
                          test_type='f',
                          ):
    '''
    Caluclates the t-test and ks-test values for the test_nuc_data_sample with
    reference to the ref_nuc_data_sample. First a vector of non-NaN entries in each
    dataset are constructed, and then these are used to calculate the test statistics.
    '''

    def _getNucPosVector(nuc_sample_name):
        bound_values = []
        unbound_values = []

        for locus in mnaseDataSamplesContainer[nuc_sample_name]:
            if nuc_pos in mnaseDataSamplesContainer[nuc_sample_name][locus]:
                # bound vector first
                if mnaseDataSamplesContainer[nuc_sample_name][locus]['in_set'][chip_sample_name] != False:
                    bound_values.append(mnaseDataSamplesContainer[nuc_sample_name][locus][nuc_pos]['pos'])
                elif mnaseDataSamplesContainer[nuc_sample_name][locus]['in_set'][chip_sample_name] == False:
                    unbound_values.append(mnaseDataSamplesContainer[nuc_sample_name][locus][nuc_pos]['pos'])

        return (np.array(bound_values), np.array(unbound_values))

    def _removeNaN(_vector):
        # remove NaN
        _vector_idx = np.isnan(_vector) == False
        _vector = _vector[_vector_idx]

        return _vector

    def _initVectorDict(nuc_sample_name):
        bound_vector, unbound_vector = _getNucPosVector(nuc_sample_name)
        bound_vector = _removeNaN(bound_vector)
        unbound_vector = _removeNaN(unbound_vector)

        vector_dict = {'bound_v' : bound_vector,
                       'unbound_v' : unbound_vector,
                       }

        return vector_dict

    test_vector_dict = _initVectorDict(test_nuc_data_sample_name)
    ref_vector_dict = _initVectorDict(ref_nuc_data_sample_name)

    f_val, p_val = f_oneway(ref_vector_dict['unbound_v'],
                            ref_vector_dict['bound_v'],
                            test_vector_dict['unbound_v'],
                            test_vector_dict['bound_v'])
    means_str = ','.join([ str(np.mean(ref_vector_dict['unbound_v'])),
                           str(np.mean(ref_vector_dict['bound_v'])),
                           str(np.mean(test_vector_dict['unbound_v'])),
                           str(np.mean(test_vector_dict['bound_v'])),
                           ] )
    #means_dict = { 'ref_unbound_mean' : np.mean(ref_vector_dict['unbound_v']),
    #               'ref_bound_mean' : np.mean(ref_vector_dict['bound_v']),
    #               'test_unbound_mean' : np.mean(test_vector_dict['unbound_v']),
    #               'test_bound_mean' : np.mean(test_vector_dict['bound_v']),
    #               }

    return (f_val, p_val, means_str)




def calcNucDataStatistics_old(mnaseDataSamplesContainer,
                          test_nuc_data_sample_name,
                          ref_nuc_data_sample_name,
                          nuc_pos,
                          chip_sample_name='none',
                          test_type='ks', flip_unbound_flag=False):
    '''
    Caluclates the t-test and ks-test values for the test_nuc_data_sample with
    reference to the ref_nuc_data_sample. First a vector of non-NaN entries in each
    dataset are constructed, and then these are used to calculate the test statistics.
    '''

    def _getNucPosVector(mnaseDataSamplesContainer, sample_name, nuc_pos, chip_sample_name='none'):
        values = []

        for locus in mnaseDataSamplesContainer[sample_name]:
            if nuc_pos in mnaseDataSamplesContainer[sample_name][locus]:
                if flip_unbound_flag:
                    if mnaseDataSamplesContainer[sample_name][locus]['in_set'][chip_sample_name] == False:
                        values.append(mnaseDataSamplesContainer[sample_name][locus][nuc_pos]['pos'])
                else:
                    if mnaseDataSamplesContainer[sample_name][locus]['in_set'][chip_sample_name] != False:
                        values.append(mnaseDataSamplesContainer[sample_name][locus][nuc_pos]['pos'])

        return np.array(values)

    test_vector = _getNucPosVector(mnaseDataSamplesContainer, test_nuc_data_sample_name, nuc_pos, chip_sample_name=chip_sample_name)

    idx = np.isnan(test_vector) == False

    test_vector = test_vector[idx]
    #print test_vector, "test"
    ref_vector = _getNucPosVector(mnaseDataSamplesContainer, ref_nuc_data_sample_name, nuc_pos, chip_sample_name=chip_sample_name)

    idx = np.isnan(ref_vector) == False
    ref_vector = ref_vector[idx]
    #print ref_vector, "ref"
    if test_type == 't':
        (t_val, p_val) = ttest_1samp(test_vector, np.mean(ref_vector))
        #print "t-test %s: %0.20f, %0.20f, %s, %s" % (nuc_pos, t_val, -np.log10(p_val), len(test_vector), len(ref_vector))
        test_mean = np.mean(test_vector)
        return (t_val, p_val, test_mean)
    elif test_type == 'ks':
        (alpha_val, p_val) = ks_2samp(test_vector, ref_vector)
        #print "ks-test %s: %0.20f, %0.20f, %s, %s" % (nuc_pos, alpha_val, -np.log10(p_val), len(test_vector), len(ref_vector))
        test_mean = np.mean(test_vector)
        return (alpha_val, p_val, test_mean)
    else:
        print "Warning: invalid test type!"
        return (None, None, None)

class ReadDensityPlot(object):
    def __init__(self, ylabel='', xlabel=''):
        self.fig = plt.figure(dpi=300, figsize=(12.8, 10.24))
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)

    def addDensity(self, chromScoreVector, left_bound, right_bound, labelText, scaleVectors=False, lineStyle='-', col='blue'):
        scoreVector = chromScoreVector
        if scaleVectors:
            scoreVector = scoreVector / max(scoreVector) *100
        plt.plot(scoreVector, linewidth=3, label=labelText, linestyle=lineStyle, color=col, alpha=0.5)
        plt.xlim((left_bound,right_bound))
        plt.legend(loc='upper left', frameon=False)
        #plt.yticks(())
        plt.tick_params(labeltop='off', labelbottom='off')

    def addLocus(self, locus_start, locus_end, locus_text, locus_color):
        #locus_text = r"\textit{%s}" % (locus_text)
        locus_text = "%s" % (locus_text)
        locus_width = locus_end - locus_start
        text_center = locus_width*.3 + locus_start
        plt.barh(0, locus_width, 100, locus_start, color=locus_color, alpha=0.3)
        plt.text(text_center, 10, locus_text, style='italic')

def scatterChIPvsFPKM(mnase_data_sample_name, chip_data_sample_name, chip_top_data_dir='chipData',
                      per_lower=.20, per_upper=.80,
                      annotated_nuc_pos_list = ['AnN1', 'AnP1', 'AnP2', 'AnP3', 'AnP4', 'AnP5', 'AnP6'],
                      dist5p=450,
                      dist3p=1000,
                      fig_fn='output/chip_fpkm_scatter.tiff'):

    mnaseDataSamplesContainer = loadAnnotatedNucsData(annotated_nuc_pos_list,
                                                    sample_name_list=[mnase_data_sample_name],
                                                    pugh_reference_genes_fn='reference/justGenesPugh.tab',
                                                    annotation_model=5,
                                                    restrict_to_verified=True)

    chipDataSamplesContainer = ChipDataSamplesContainer(top_sample_dir=chip_top_data_dir,
                                                   chip_sample_name_list=[chip_data_sample_name],
                                                   #bed_fn_suffix='_summits.bed')
                                                   bed_fn_suffix='_subpeaks.summits.bed')

    chipDataSamplesContainer.loadAllSummits()
    mnaseDataSamplesContainer.annotateChipBindingAllSamples(chipDataSamplesContainer, feature='TSS', dist5p=dist5p, dist3p=dist3p)

    mnaseDataSamplesContainer[mnase_data_sample_name].annotateFPKM()

    x = np.array([mnaseDataSamplesContainer[mnase_data_sample_name][locus]['chip_enrichment'][chip_data_sample_name]
               for locus in mnaseDataSamplesContainer[mnase_data_sample_name]])
    y = np.array([mnaseDataSamplesContainer[mnase_data_sample_name][i]['fpkm']
               for i in mnaseDataSamplesContainer[mnase_data_sample_name].keys()])

    fig = plt.figure(dpi=300, figsize=(12.8,10.24))
    ax = plt.subplot(111)
    box = ax.get_position()
    #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    plt.scatter(np.log2(x), np.log2(y), color='grey', s=96, alpha=.3, edgecolor='', label='All Expressing Loci')

    bound_idx = np.array([mnaseDataSamplesContainer[mnase_data_sample_name][locus]['in_set'][chip_data_sample_name]
                       for locus in mnaseDataSamplesContainer[mnase_data_sample_name]]) > 0
    plt.scatter(np.log2(x[bound_idx]), np.log2(y[bound_idx]),
                color='orange', s=96, alpha=0.6, edgecolor='',
                #label='Bound in %s' % (chipDataSamplesContainer[chip_data_sample_name].common_id) )
                label='Bound by Yta7')

    print np.corrcoef(np.log2(x[bound_idx]), np.log2(y[bound_idx]))

    top_fpkm_idx = np.log2(y) > scoreatpercentile(np.log2(y), per_upper)
    plt.scatter(np.log2(x[bound_idx&top_fpkm_idx]),
                np.log2(y[bound_idx&top_fpkm_idx]),
                color='red', s=96, alpha=0.6,
                edgecolor='',
                #label='Bound in %s\nand in Top %s%% of FPKM' %
                #(chipDataSamplesContainer[chip_data_sample_name].common_id, int(100-per_upper)) )
                label='Bound by Yta7 -- Top %s%% of FPKM' % (int(100-per_upper)) )

    print np.corrcoef(np.log2(x[bound_idx&top_fpkm_idx]), np.log2(y[bound_idx&top_fpkm_idx]))

    bottom_fpkm_idx = np.log2(y) < scoreatpercentile(np.log2(y), per_lower)
    plt.scatter(np.log2(x[bound_idx&bottom_fpkm_idx]),
                np.log2(y[bound_idx&bottom_fpkm_idx]),
                color='blue', s=96, alpha=0.6,
                edgecolor='',
                #label='Bound in %s\nand in Bottom %s%% of FPKM' %
                #(chipDataSamplesContainer[chip_data_sample_name].common_id, int(per_lower)) )
                label='Bound by Yta7 -- Bottom %s%% of FPKM' % (int(per_lower)) )

    print np.corrcoef(np.log2(x[bound_idx&bottom_fpkm_idx]), np.log2(y[bound_idx&bottom_fpkm_idx]))

    plt.xlim((5,11))
    plt.ylim((0,16))

    plt.xlabel('Yta7 ChIP Enrichment')
    plt.ylabel('FPKM')

    plt.tick_params(pad=10, labeltop='off', labelright='off')

    fontP = FontProperties()
    fontP.set_size('small')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), prop=fontP, scatterpoints=1, frameon=False, ncol=2)

    plt.savefig(fig_fn)

def makeBoxPlot(mnase_data_sample_name='DEX_WT_alpha', chip_data_sample_name='DEX_YTA7_IP_alpha', dist5p=450, dist3p=1000, fig_fn='output/Plt.Figure2C.tiff'):
    '''
    '''

    chipDataSamplesContainer = ChipDataSamplesContainer(top_sample_dir='chipData',
                                                        chip_sample_name_list=[chip_data_sample_name],
                                                        #bed_fn_suffix='_summits.bed')
                                                        bed_fn_suffix='_subpeaks.summits.bed')

    chipScoreVectorDict = loadChromScoreVectorDict(chipDataSamplesContainer[chip_data_sample_name].bedGraph_fn)

    annotated_nuc_pos_list = ['AnN1', 'AnP1', 'AnP2', 'AnP3', 'AnP4', 'AnP5', 'AnP6']
    mnaseDataSamplesContainer = loadAnnotatedNucsData(annotated_nuc_pos_list,
                                                    sample_name_list=[mnase_data_sample_name],
                                                    pugh_reference_genes_fn='reference/justGenesPugh.tab',
                                                    restrict_to_verified=True,
                                                    annotation_model=5)

    mnaseScoreVectorDict = loadChromScoreVectorDict(mnaseDataSamplesContainer[mnase_data_sample_name].bedGraph_fn)

    chipDataSamplesContainer.loadAllSummits()

    mnaseDataSamplesContainer.annotateChipBindingAllSamples(chipDataSamplesContainer, feature='TSS', dist5p=dist5p, dist3p=dist3p)

    for mnase_data_sample_name in mnaseDataSamplesContainer:
        mnaseDataSamplesContainer[mnase_data_sample_name].annotateFPKM()

    bound_idx = np.array([mnaseDataSamplesContainer[mnase_data_sample_name][locus]['in_set'][chip_data_sample_name] != False
                       for locus in mnaseDataSamplesContainer[mnase_data_sample_name]])
    unbound_idx = np.array([mnaseDataSamplesContainer[mnase_data_sample_name][locus]['in_set'][chip_data_sample_name] == False
                         for locus in mnaseDataSamplesContainer[mnase_data_sample_name]])
    expr_idx = np.array(np.isnan(np.array([mnaseDataSamplesContainer[mnase_data_sample_name][locus]['fpkm']
                                  for locus in mnaseDataSamplesContainer[mnase_data_sample_name]])) == False)

    # add .01 for log calculation to avoid log(0);
    fpkms = np.log2(np.array([mnaseDataSamplesContainer[mnase_data_sample_name][locus]['fpkm']
                        for locus in mnaseDataSamplesContainer[mnase_data_sample_name]])+.01)
    bound_x = fpkms[bound_idx & expr_idx]
    unbound_y = fpkms[unbound_idx & expr_idx]
    #allexpr_y = fpkms[expr_idx]

    plt.figure(dpi=300, figsize=((1*10.24, 1*7.68)))
    plt.boxplot([bound_x, unbound_y], vert=1, notch=1, bootstrap=10000, widths=.2, sym='')
    #plt.boxplot([bound_x, allexpr_y], vert=1, notch=1, bootstrap=10000, widths=.2, sym='')
    plt.xticks(np.arange(4), ('', 'YTA7-Bound Loci', 'Unbound Loci', ''))
    #plt.xticks(np.arange(4), ('', 'YTA7 Bound Loci', 'All Expressing Loci', ''))
    plt.ylabel('FPKM')

    plt.tick_params(pad=10, labeltop='off', labelright='off')

    if not fig_fn.startswith('output/'):
        fig_fn = 'output/' + fig_fn
    if not fig_fn.endswith('.tiff'):
        fig_fn = fig_fn + '.tiff'
    plt.savefig(fig_fn)

    (t_val, p_val) = ttest_ind(bound_x, unbound_y)
    with open(fig_fn + '.stats', 'w') as out_fh:
        out_fh.write('# Bound Loci\t%s\n' % (len(bound_x)))
        out_fh.write('# Unbound Loci\t%s\n' % (len(unbound_y)))
        out_fh.write('t-test score (independent samples)\t%s\n' % (t_val))
        out_fh.write('p-val\t< %s\n' % (p_val))
