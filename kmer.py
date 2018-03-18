from Bio import SeqIO
import pandas as pd
import argparse
import os


class Kmer:
    '''
    kmer creater

    :param sequence:   whole sequence (1 string from multi-fasta)
    :param kmer_name:  kmer name
    :param locus_name: record id from input fasta (all that after '>' symbol)

    :return:  <Kmer_object>
    '''

    counter = 0
    sequence = ''
    kmer_name = ''
    chr = ''
    def __init__(self, sequence, kmer_name, locus_name):
        self.sequence = sequence
        self.kmer_name = kmer_name
        self.locus_name = locus_name
        self.length = len(str(self))
        self.chords = []

    def increase(self):
        self.counter += 1

    def locate(self, locus):
        self.chords.append(locus)

    def infogram(self):
        self.full_descr = pd.DataFrame({'k_mer_name': self.kmer_name,
                                        'Locus': ';'.join(['{0}-{1}'.format(i[0], i[1]) for i in self.chords]),
                                        'sequence': self.sequence,
                                        'Chromosome': [self.locus_name] * len(self.chords),
                                        'Count': self.counter}).head(1)
        return self.full_descr


def wheel(loc, seq, k):
    '''
    dictionary writer

    :param
    :param loc:   locus name
    :param seq:   sequence name
    :param k:     k-mer length

    :return kmer_dict:  {sequence: {kmer_name: <Kmer_object>}
    '''

    seq_lng = len(seq)
    kmer_dict = {}
    for index in range(seq_lng - k + 1):
        current_kmer = seq[index:(index + k)]
        if current_kmer in kmer_dict:
            kmer_dict[current_kmer].increase()
        else:
            kmer_dict[current_kmer] = Kmer(seq, current_kmer, loc)
            kmer_dict[current_kmer].increase()
        kmer_dict[current_kmer].locate([index, (index + k)])
    return kmer_dict

if __name__ == '__main__':

    def str2bool(v):
        '''
        string variable converter

        :param v: string value, that should be boolean

        :return:  True/False
        '''

        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = argparse.ArgumentParser(description='''You can use this script for k-mer mapping.
                                               Output format is simple pandas dataframe.
                                               If you want you can get only subset from whole dataframe''')

    parser.add_argument('-i', '--input', help='Paste your path to input file here.', type=str)
    parser.add_argument('-o', '--output', help='Paste your path to output file here.',
                       default=os.getcwd()+'/output', type=str)
    parser.add_argument('-k', '--kmer', help='k-mer length',
                       default=23, type=int)
    parser.add_argument('-s', '--subset', help='''You can choose k-mer name (for example \'AATGC\'),
    which will be taken as whole dataframe's subset''', type=str)
    parser.add_argument('-mf', '--mfreq', help='''Return dataframe with most popular kmer''',
                        default=True, type=str2bool)
    parser.add_argument('-bd', '--big_data', help='If you not want get whole dataframe, you should use 0',
                        type=str2bool)

    args = parser.parse_args()


    i = args.input
    o = args.output
    k = args.kmer
    s = args.subset
    mf = args.mfreq
    bd = args.big_data

    if i is None:
        raise argparse.ArgumentTypeError('You forgot your path to file (try -i)')
    i = os.path.abspath(i)
    if 'k' not in globals():
        raise argparse.ArgumentTypeError('You forgot paste -k')

    kmer_dict, sequence_dict = {}, {}
    for record in SeqIO.parse(i, 'fasta'):
        sequence_dict[record.id] = str(record.seq)
        for key, value in sequence_dict.items():
            kmer_dict[key] = wheel(key, value, k)

    if mf == True:
        mf_dict = {}
        for kk, vv in kmer_dict.items():
            list_i = [('A', 0)]
            for kkk, vvv in vv.items():
                if vvv.counter < list_i[0][1]:
                    continue
                elif vvv.counter == list_i[0][1]:
                    list_i.append((kkk, vvv.counter))
                elif vvv.counter > list_i[0][1]:
                    list_i = [(kkk, vvv.counter)]
            mf_dict[kk] = list_i

        with open('{}_mf.csv'.format(o), 'a') as wr:
            for kk, vv in mf_dict.items():
                for index in vv:
                    wr.write('{},{},{}\n'.format(kk, index[0], index[1]))

    if s is not None:
        list_d = []
        for value in kmer_dict.values():
            for val in value.values():
                if val.kmer_name == s:
                    list_d.append(val.infogram())
                    break
        summ_df = pd.concat(list_d)
        summ_df.to_csv('{}_selected.csv'.format(o))

    if bd == True:
        list_of_df = []
        for value in kmer_dict.values():
            for val in value.values():
                list_of_df.append(val.infogram())
        summ_df = pd.concat(list_of_df)
        summ_df.to_csv('{}_bd.csv'.format(o))

