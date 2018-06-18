# Run EXP1 with 25k reads
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP1_W_R1.25000.fastq -f2 ./example/fastq/EXP1_W_R2.25000.fastq -o ./example/output/ -n EXP1_W.25k
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP1_Q_R1.25000.fastq -f2 ./example/fastq/EXP1_Q_R2.25000.fastq -o ./example/output/ -n EXP1_Q.25k
python recYnH.py score -m1 ./example/output/EXP1_W.25k -m2 ./example/output/EXP1_Q.25k -o ./example/output/ -n EXP1.25k

# Run EXP2 with 25k reads
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP2_W_R1.25000.fastq -f2 ./example/fastq/EXP2_W_R2.25000.fastq -o ./example/output/ -n EXP2_W.25k
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP2_Q_R1.25000.fastq -f2 ./example/fastq/EXP2_Q_R2.25000.fastq -o ./example/output/ -n EXP2_Q.25k
python recYnH.py score -m1 ./example/output/EXP2_W.25k -m2 ./example/output/EXP2_Q.25k -o ./example/output/ -n EXP2.25k

# Run EXP3 with 25k reads
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP3_W_R1.25000.fastq -f2 ./example/fastq/EXP3_W_R2.25000.fastq -o ./example/output/ -n EXP3_W.25k
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP3_Q_R1.25000.fastq -f2 ./example/fastq/EXP3_Q_R2.25000.fastq -o ./example/output/ -n EXP3_Q.25k
python recYnH.py score -m1 ./example/output/EXP3_W.25k -m2 ./example/output/EXP3_Q.25k -o ./example/output/ -n EXP3.25k

# Merge 5k reads
python recYnH.py merge -i ./example/output/EXP1.25k.nis.txt ./example/output/EXP2.25k.nis.txt ./example/output/EXP3.25k.nis.txt -o ./example/output/ -n EXP.25k.avgIS



# Run EXP1 with 5k reads
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP1_W_R1.5000.fastq -f2 ./example/fastq/EXP1_W_R2.5000.fastq -o ./example/output/ -n EXP1_W.5k
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP1_Q_R1.5000.fastq -f2 ./example/fastq/EXP1_Q_R2.5000.fastq -o ./example/output/ -n EXP1_Q.5k
python recYnH.py score -m1 ./example/output/EXP1_W.5k -m2 ./example/output/EXP1_Q.5k -o ./example/output/ -n EXP1.5k

# Run EXP2 with 5k reads
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP2_W_R1.5000.fastq -f2 ./example/fastq/EXP2_W_R2.5000.fastq -o ./example/output/ -n EXP2_W.5k
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP2_Q_R1.5000.fastq -f2 ./example/fastq/EXP2_Q_R2.5000.fastq -o ./example/output/ -n EXP2_Q.5k
python recYnH.py score -m1 ./example/output/EXP2_W.5k -m2 ./example/output/EXP2_Q.5k -o ./example/output/ -n EXP2.5k

# Run EXP3 with 5k reads
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP3_W_R1.5000.fastq -f2 ./example/fastq/EXP3_W_R2.5000.fastq -o ./example/output/ -n EXP3_W.5k
python recYnH.py align -i1 ./example/db/roth2016_control_set_plus_control.fa -f1 ./example/fastq/EXP3_Q_R1.5000.fastq -f2 ./example/fastq/EXP3_Q_R2.5000.fastq -o ./example/output/ -n EXP3_Q.5k
python recYnH.py score -m1 ./example/output/EXP3_W.5k -m2 ./example/output/EXP3_Q.5k -o ./example/output/ -n EXP3.5k

# Merge 5k reads
python recYnH.py merge -i ./example/output/EXP1.5k.nis.txt ./example/output/EXP2.5k.nis.txt ./example/output/EXP3.5k.nis.txt -o ./example/output/ -n EXP.5k.avgIS

