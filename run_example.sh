
# Run EXP2 with 25k reads
#python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.-150.fa -f1 ./example/fastq/EXP2_W_R1.25000.fastq -f2 ./example/fastq/EXP2_W_R2.25000.fastq -o ./example/output/ -n EXP2_W.25k
#python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.-150.fa -f1 ./example/fastq/EXP2_Q_R1.25000.fastq -f2 ./example/fastq/EXP2_Q_R2.25000.fastq -o ./example/output/ -n EXP2_Q.25k
#python recYnH.py score -m1 ./example/output/EXP2_W.25k -m2 ./example/output/EXP2_Q.25k -o ./example/output/ -n EXP2.25k

# Run EXP3 with 25k reads
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.-150.fa -f1 ./example/fastq/EXP3_W_R1.25000.fastq -f2 ./example/fastq/EXP3_W_R2.25000.fastq -o ./example/output/ -n EXP3_W.25k
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.-150.fa -f1 ./example/fastq/EXP3_Q_R1.25000.fastq -f2 ./example/fastq/EXP3_Q_R2.25000.fastq -o ./example/output/ -n EXP3_Q.25k
python recYnH.py score -m1 ./example/output/EXP3_W.25k -m2 ./example/output/EXP3_Q.25k -o ./example/output/ -n EXP3.25k
