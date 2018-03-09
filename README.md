### kmer_finder
You can clone whole repo by standart command:
```
git clone https://github.com/rostkick/kmer_finder.git
```
Also if you need to save it **without commit-log** you schould use:
```
git clone —depth=1 https://github.com/rostkick/kmer_finder.git
```
If you want to save **only one** file (for example without saving testing_data) you should:
  1) copy link of favorite file
  2) go to https://minhaskamal.github.io/DownGit/#/home, paste it, and create Download Link.

Keys description:  
`-i, --input`     Paste your path to fasta file here.  
`-o, --output`    Paste your path to output file here.  
`-k, --kmer`      Specify kmer length.  
`-s, --subset`    Paste kmer sequence here, which will be taken as subset.  
`-mf, --mfreq`    Use 1 if you want to get frequency and name of most frequence.   
`-bd, --big_data` Use 1 if you want to get whole dataframe with all kmers.

Simple usage:  
1. Getting most freq kmer name.
```
$ python3 kmer.py -i /home/usr/example.fasta -o /home/usr/out -k 23 -mf 1
```
2. Getting subset.  
```
$ python kmer.py -i /home/usr/example.fasta -o /home/usr/out -k 23 -s 'ATGTAGT..TAGT'
```
3. Getting whole dataframe.  
```
$ python kmer.py -i /home/usr/example.fasta -o /home/usr/out -k 23 -bd 1
```
