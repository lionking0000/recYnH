[What is recYnH?]

This program is taking recYnH sequencing files and generating recYnH interaction score matrix to correponding genes.

We assumed that the sequencing format is the same as we described in Nature Communication 2018 Paper.

# After setting proper parameter in recYnH.py file
# User can run the program as following command

usage: recYnH.py [-h] COMMAND ...

recYnH program

positional arguments:
  COMMAND     sub-command help
    align     a help for align
    merge     a help for merge

optional arguments:
  -h, --help  show this help message and exit

Commands:
  align      Align the FASTQ sequencing files into bait and prey sequences to generate interaction matrix
  merge      Merge two interaction matries to generate an interaction score matrix

Run 'recYnH.py COMMAND --help' for more information on a command.

---------------------------------------------------------------------------------------------------------------------

# python recYnH.py align -i1 [Input fasta file] -f1 [Input fastq file1] -f2 [Input fastq file2] -o [Output folder]

sh Y2H_Blastn.sh 2017-11-17_MiSeq S1_WD_R1 S1_WD_R2 ../data/A463-MGj69.RBP-MAP.-150 S1_W > qjobs/qjob_Sebastian_2017-11-17_MiSeq_S1_W1.sh

# Typical output

[ Starting recY2H Analysis ]
[ Read No-selection Sequence File ] ./data/XXXX.fastq.gz
[ Read Selection Sequence File ] ./data/YYYY.fastq.gz
[ Calculating recY2H Scores ]
[ Finishing recY2H Analysis ]


# It will gives recY2H Scores for all the sequences


[Output Files]


# How to use Docker

[Build image]

docker build -t recynh .

or

sh docker_build.sh

[Test code]

sh docker_align.sh -i1 /share/db/A463-MGj69.RBP-MAP.-150.fa -f1 /share/fastq/S1_WD_R1.300000.fastq.gz -f2 /share/fastq/S1_WD_R2.300000.fastq.gz -o /share/output/2017-11-17_MiSeq -n S1_WD
sh docker_align.sh -i1 /share/db/A463-MGj69.RBP-MAP.-150.fa -f1 /share/fastq/S2_QD_R1.300000.fastq.gz -f2 /share/fastq/S2_QD_R2.300000.fastq.gz -o /share/output/2017-11-17_MiSeq -n S2_QD

sh docker_merge.sh -m1 /share/output/2017-11-17_MiSeq/S1_WD.ppi.txt -m2 /share/output/2017-11-17_MiSeq/S2_QD.ppi.txt

or

sh docker_merge.sh -m1 /share/output/example/S1_W.ppi.txt -m2 /share/output/example/S2_Q.ppi.txt


[Run image]

docker run -d -v ~/myshared:/share --name myynh recynh tail -f /dev/null

docker run -d -v /Users/jyang/Dropbox/Code/recYnH/share:/share --name myynh lionking0000/recynh:0.1 tail -f /dev/null

docker exec myynh2 Y2H_Blastn.sh 2017-11-17_MiSeq S1_WD_R1.10000 S1_WD_R2.10000 ./db/A463-MGj69.RBP-MAP.-150 S1_W


/Users/jyang/Dropbox/Code/recYnH/share is a shared volume used for convenience for placing input and output files

The command above keeps container running under the name myynh


# running interactive mode
docker run -v /Users/jyang/Dropbox/Code/recYnH/share:/share -ti lionking0000/recynh:0.1 /bin/bash

* example
Y2H_Blastn.sh 2017-11-17_MiSeq S1_WD_R1.10000 S1_WD_R2.10000 ./db/A463-MGj69.RBP-MAP.-150 S1_W


[Excute command]

docker exec myynh rec-Ynh.py -i /share/reads/SRR493366.fastq -o /share/out/test


# References & Related information

https://github.com/toniher/ELMSeq/blob/master/Dockerfile



sh Y2H_Blastn.sh 2017-11-17_MiSeq S1_WD_R1 S1_WD_R2 ../data/A463-MGj69.RBP-MAP.-150 S1_W > qjobs/qjob_Sebastian_2017-11-17_MiSeq_S1_W1.sh


## in the docker

cd /share
mkdir output/2017-11-17_MiSeq/
cutadapt -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o output/2017-11-17_MiSeq//S1_WD_R1.fastq -p output/2017-11-17_MiSeq//S1_WD_R2.fastq -m 15 --discard-untrimmed ./fastq/S1_WD_R1.10000.fastq.gz ./fastq/S1_WD_R2.10000.fastq.gz

main.py fastq_to_fasta output/2017-11-17_MiSeq//S1_WD_R1.fastq > output/2017-11-17_MiSeq//S1_WD_R1.fa
main.py fastq_to_fasta output/2017-11-17_MiSeq//S1_WD_R2.fastq > output/2017-11-17_MiSeq//S1_WD_R2.fa
blastn -db ./db/A463-MGj69.RBP-MAP.-150.fa  -query output/2017-11-17_MiSeq//S1_WD_R1.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/2017-11-17_MiSeq//S1_WD_R1.blastn
blastn -db ./db/A463-MGj69.RBP-MAP.-150.fa  -query output/2017-11-17_MiSeq//S1_WD_R2.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/2017-11-17_MiSeq//S1_WD_R2.blastn
main.py BLASTN ./db/A463-MGj69.RBP-MAP.-150.fa output/2017-11-17_MiSeq//S1_WD_R1.blastn output/2017-11-17_MiSeq//S1_WD_R2.blastn > output/2017-11-17_MiSeq//S1_W.ppi.txt

