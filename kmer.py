from Bio import SeqIO
import pandas as pd
import argparse

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='''You can use this script for k-mer mapping.
                                               Output format is simple pandas dataframe.
                                               If you want you can get only subset from whole dataframe''')

   parser.add_argument('-i', '--input', help='Paste your path to input file here.', default=None, type=str)
   parser.add_argument('-o', '--output', help='Paste your path to output file here.', default='output.txt', type=str)
   parser.add_argument('-k', '--kmer', help='k-mer length', default=None, type=int)
   parser.add_argument('-s', '--subset', help='''You can choose k-mer name (for example \'AATGC\'),
                                           which will be taken as whole dataframe's subset''', type=str)
   parser.add_argument('-mf', '--mfreq', help='''Return dataframe with most popular kmer''', default=0, type=bool)
   parser.add_argument('-bd', '--big_data', help='''If you want get whole dataframe, you schould use 1''', default=0)

   args = parser.parse_args()

   i = args.input
   o = args.output
   k = args.kmer
   s = args.subset
   mf = args.mfreq
   bd = args.big_data

if i is None:
    raise InputError('You forgot your path to file (try -i)')
if k is None:
    raise KmerError('You forgot paste -k')



##print(i)
class Kmer:
    counter = 0
    sequence = ''
    chr = ''
    def __init__(self, kmer_name, locus_name):
        self.sequence = kmer_name
        self.locus_name = locus_name
        self.length = len(str(self))
        self.loci_list = []

    def increase(self):
        self.counter += 1

    def locate(self, locus):
        self.loci_list.append(locus)

    def infogram(self):
        self.info = pd.DataFrame({'k_mer_name': self.sequence, 'Locus': self.loci_list,
                                  'Chromosome': [self.locus_name] * len(self.loci_list)})
        return self.info


def wheel(loc, seq, k):
    seq_lng = len(seq)
    kmer_dict = {}
    for index in range(seq_lng - k + 1):
        current_kmer = seq[index:(index + k)]
        if current_kmer in kmer_dict:
            kmer_dict[current_kmer].increase()
        else:
            kmer_dict[current_kmer] = Kmer(current_kmer, loc)
            kmer_dict[current_kmer].increase()
        kmer_dict[current_kmer].locate([index, (index + k)])
    return kmer_dict



sequence_list = []
for record in SeqIO.parse(i, 'fasta'):
    sequence_list.append((str(record.id), str(record.seq).  upper()))
    kmer_list = wheel(sequence_list[0][0], sequence_list[0][1], int(k))


if mf != 0:
    i = [1]
    for key, value in kmer_list.items():
        if kmer_list[key].counter < i[0]:
            continue
        else:
            i = (kmer_list[key].counter, key)
            with open(o+'most_freq.txt', 'w') as w:
                w.write(str(i))

if s != None:
    inf=kmer_list[s].infogram()
    print(inf)
    inf.to_csv(o+'target_seq'+'.csv')

if bd != 0:
    list_of_df = []
    for val in kmer_list.values():
        list_of_df.append(val.infogram())
    summ_df = pd.concat(list_of_df)
    summ_df.to_csv(o+'whole_df.csv')

