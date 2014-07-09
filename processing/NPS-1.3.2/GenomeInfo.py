#!/home/liulab/shin/python/bin/python
# zy@jimmy.harvard.edu

genome = {'hg18' : 'Human', 'mm9': 'Mouse', 'dm2': 'Fly' }
          	
chrlength = {}
validlength = {}
controlgenome = {}

# hg18
chrlength['hg18'] = {'chr1' : 247249719,  'chr2' : 242951149,  'chr3' : 199501827,
                        'chr4' : 191273063,  'chr5' : 180857866,  'chr6' : 170899992,
                        'chr7' : 158821424,  'chr8' : 146274826,  'chr9' : 140273252,
                        'chr10' : 135374737, 'chr11' : 134452384, 'chr12' : 132349534,
                        'chr13' : 114142980, 'chr14' : 106368585, 'chr15' : 100338915,
                        'chr16' : 88827254,  'chr17' : 78774742,  'chr18' : 76117153,
                        'chr19' : 63811651,  'chr20' : 62435964,  'chr21' : 46944323,
                        'chr22' : 49691432,  'chrX' : 154913754,  'chrY' : 57772954 }

validlength['hg18'] = {'chr1' : 224999719,  'chr2' : 237709794,  'chr3' : 194704822,
                        'chr4' : 187297063,  'chr5' : 177702766,  'chr6' : 167273991,
                        'chr7' : 154952424,  'chr8' : 142612826,  'chr9' : 120143252,
                        'chr10' : 131624728, 'chr11' : 131130753, 'chr12' : 130303032,
                        'chr13' : 95559980, 'chr14' : 88290585, 'chr15' : 81341915,
                        'chr16' : 78884752,  'chr17' : 77800220,  'chr18' : 74656155,
                        'chr19' : 55785651,  'chr20' : 59505253,  'chr21' : 34170106,
                        'chr22' : 34851311,  'chrX' : 151058754,  'chrY' : 25652954 }

controlgenome['hg18'] = 2400000000

# hg19
chrlength['hg19'] = {"chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,\
                     "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,\
                     "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,\
                     "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,\
                     "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,\
                     "chr21": 48129895, "chr22": 51304566, "chrX": 155270560, "chrY": 59373566, "chrM": 16571}

# D.melanogaster genome info
chrlength['dm2']={'chr2L':22407834, 'chr2R':20766785, 'chr2h':1694122, 'chr3L':23771897,\
                  'chr3R':27905053, 'chr3h':2955737, 'chr4':1281640, 'chr4h':88110,\
                  'chrX':22224390, 'chrXh':359526, 'chrYh':396896, 'chrU':8724946,'chrM':19517}

# yeast genome from SGD
chrlength['sacCer3'] = {'chrI':230218, 'chrII':813184, 'chrIII':316620, 'chrIV':1531933, \
                    'chrV':576874, 'chrVI':270161, 'chrVII':1090940, 'chrVIII':562643, \
                    'chrIX':439888, 'chrX':745751, 'chrXI':666816, 'chrXII':1078177, \
                    'chrXIII':924431, 'chrXIV':784333, 'chrXV':1091291, 'chrXVI':948066, \
                    'chrM':85779} 

chrlength['sc'] = {'chr1':230208, 'chr2':813178, 'chr3':316617, 'chr4':1531919, \
                    'chr5':576869, 'chr6':270148, 'chr7':1090947, 'chr8':562643, \
                    'chr9':439885, 'chr10':745741, 'chr11':666454, 'chr12':1078175, \
                    'chr13':924429, 'chr14':784333, 'chr15':1091289, 'chr16':948062, \
                    'chrM':85779} 

chrlength['mm9'] = {'chr1':197195432, 'chr2':181748087, 'chr3':159599783, 'chr4':155630120,\
                    'chr5':152537259, 'chr6':149517037, 'chr7':152524553, 'chr8':131738871,\
                    'chr9':124076172, 'chr10':129993255, 'chr11':121843856, 'chr12':121257530,\
                    'chr13':120284312, 'chr14':125194864, 'chr15':103494974, 'chr16':98319150,\
                    'chr17':95272651, 'chr18':90772031, 'chr19':61342430, 'chrX':166650296,\
                    'chrY':15902555, 'chrM':16299}          

chrlength['mm8'] = {'chr1':197069962, 'chr2':181976762, 'chr3':159872112, 'chr4':155029701,\
                    'chr5':152003063, 'chr6':149525685, 'chr7':145134094, 'chr8':132085098,\
                    'chr9':124000669,  'chr10':129959148,  'chr11':121798632, 'chr12':120463159,\
                    'chr13':120614378, 'chr14':123978870, 'chr15':103492577, 'chr16':98252459,\
                    'chr17':95177420, 'chr18':90736837, 'chr19':61321190, 'chrX':165556469,\
                    'chrY':16029404, 'chrM':16299}        


class GenomeInfo(object):
    def __init__(self, genomeversion):
        self.genomeversion = genomeversion
        if (not genome.has_key(self.genomeversion)):
            raise 'The genome version %s is not supported' %(self.genomeversion)
        
    def OrganismName(self):
        if genome.has_key(self.genomeversion):
            self.organismname = genome[self.genomeversion]
        else:
            raise 'The genome version %s is not supported' %(self.genomeversion)
        
        return self.organismname
    
    def GetChrLength(self, type = 'total'):
        if type == 'total':
            if chrlength.has_key(self.genomeversion):
                return chrlength[self.genomeversion]
            else:
                raise 'The genome version %s is not supported' %(self.genomeversion)
        elif type == 'valid':
            if validlength.has_key(self.genomeversion):
                return validlength[self.genomeversion]
            else:
                raise 'The genome version %s is not supported' %(self.genomeversion)
        else:
            raise 'The type %s is invalide (only \'total\' and \'valid\' are valid type)' %(self.genomeversion)

    def GetBgLength(self):
        if controlgenome.has_key(self.genomeversion):
            return controlgenome[self.genomeversion]
        else:
            raise 'The genome version %s is not supported' %(self.genomeversion)

if __name__ == '__main__':
    a = GenomeInfo('hg18')
    print a.OrganismName()
    print a.GetChrLength()
    print a.GetBgLength()
    
