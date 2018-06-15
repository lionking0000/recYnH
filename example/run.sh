#!/bin/bash

#===============================================================================
# Y2H_Blastn.sh
# created 2017.07.03
# written by Jae-Seong Yang
#
# Blastn short-read based ppi map
#===============================================================================

# READ: https://www.ncbi.nlm.nih.gov/books/NBK279668/

# $1 : directory
# $2 : File1 (without .fastq.gz)
# $3 : File2 (without .fastq.gz)
# $4 : Fasta (without .fa)
# $5 : Outname (prefix only)
# $6 : Output Folder

#===============================================================================
# Need to make blastn database (only -100 bp containing) as following instructions before running this script
#
# python main.py generate_last_nts ../data/P170_4_library_MGj5615-Jun-2017122014.fa 100 > ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa
# makeblastdb -in P170_4_library_MGj5615-Jun-2017122014.-100.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/E375pLT-20170706.fa 100 > ../data/E375pLT-20170706.-100.fa 
# makeblastdb -in ../data/E375pLT-20170706.-100.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/roth2016_control_set_plus_control.fa 100 > ../data/roth2016_control_set_plus_control.-100.fa 
# makeblastdb -in ../data/roth2016_control_set_plus_control.-100.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/EGFR_entry_clones.fa  100 > ../data/EGFR_entry_clones.-100.fa
# makeblastdb -in EGFR_entry_clones.-100.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/roth2016_control_set_plus_control.fa 150 > ../data/roth2016_control_set_plus_control.-150.fa 
# makeblastdb -in ../data/roth2016_control_set_plus_control.-150.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/A463-MGj69.RBP-MAP.fa 150 > ../data/A463-MGj69.RBP-MAP.-150.fa 
# makeblastdb -in ../data/A463-MGj69.RBP-MAP.-150.fa -dbtype nucl
#===============================================================================


#===============================================================================
# JOB_LIST (what have been done)
#
### output/2016-12-22_MiSeq/Friedrich/
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17543_S1_R1 17543_S1_R2 ../data/roth2016_control_set_plus_control.-100 17543_S1 > qjobs/qjob_2016-12-22_MiSeq_S1.sh # ==> 1968532 (Roth -W)
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17544_S2_R1 17544_S2_R2 ../data/roth2016_control_set_plus_control.-100 17544_S2 > qjobs/qjob_2016-12-22_MiSeq_S2.sh # ==> 1541502 (Roth -A) from 2418986 (trimed seq) from xxx
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17545_S3_R1 17545_S3_R2 ../data/roth2016_control_set_plus_control.-100 17545_S3 > qjobs/qjob_2016-12-22_MiSeq_S3.sh # ==> 2190445 (Roth Seaprep -W)
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17546_S4_R1 17546_S4_R2 ../data/roth2016_control_set_plus_control.-100 17546_S4 > qjobs/qjob_2016-12-22_MiSeq_S4.sh # ==> 1877953  (Roth Seaprep -A)
#
### output/2017-02-04_MiSeq/Friedrich/ ## (not use because of contamination)
#
#
### output/2017-02-22_MiSeq/Friedrich/
#
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S1_R1 S1_R2 ../data/roth2016_control_set_plus_control.-100 S1 > qjobs/qjob_2017-02-22_MiSeq_S1.sh # (Roth -W) # ==> 960245
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S2_R1 S2_R2 ../data/roth2016_control_set_plus_control.-100 S2 > qjobs/qjob_2017-02-22_MiSeq_S2.sh # (Roth -A) # ==> 2024754
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S3_R1 S3_R2 ../data/roth2016_control_set_plus_control.-100 S3 > qjobs/qjob_2017-02-22_MiSeq_S3.sh # (Roth -Q) # ==> 813912
#
### output/2017-03-03_MiSeq/Friedrich (Roth)
#
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S1_R1 S1_R2 ../data/roth2016_control_set_plus_control.-100 S1 > qjobs/qjob_2017-03-03_MiSeq_S1.sh # (Roth -W) # ==> 1414435
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S2_R1 S2_R2 ../data/roth2016_control_set_plus_control.-100 S2 > qjobs/qjob_2017-03-03_MiSeq_S2.sh # (Roth -A) # ==> 969159
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S3_R1 S3_R2 ../data/roth2016_control_set_plus_control.-100 S3 > qjobs/qjob_2017-03-03_MiSeq_S3.sh # (Roth -Q) # ==> 1255002
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S4_R1 S4_R2 ../data/roth2016_control_set_plus_control.-100 S4 > qjobs/qjob_2017-03-03_MiSeq_S4.sh # (Roth Aonly -W) to check toxicity # ==> 23
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S5_R1 S5_R2 ../data/roth2016_control_set_plus_control.-100 S5 > qjobs/qjob_2017-03-03_MiSeq_S5.sh # (Roth Bonly -W) to check toxicity # ==> 169
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S6_R1 S6_R2 ../data/roth2016_control_set_plus_control.-100 S6 > qjobs/qjob_2017-03-03_MiSeq_S6.sh # (Roth Bonly -A) to check auto-activation # ==> 1296
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S7_R1 S7_R2 ../data/roth2016_control_set_plus_control.-100 S7 > qjobs/qjob_2017-03-03_MiSeq_S7.sh # (Roth Seaprep -Q) # ==> 734766
# 
# 2017-06-08_MiSeq (P170 toxicity test and Roth Seaprep)
#  
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S49_R1 S49_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S49 > qjobs/qjob_2017-06-08_MiSeq_S49.sh # (P170 Bonly -W) to check toxicity           # ==> 9299 
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S50_R1 S50_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S50 > qjobs/qjob_2017-06-08_MiSeq_S50.sh # (P170 Bonly -W-A) to check auto-activation  # ==> 7991
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S51_R1 S51_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S51 > qjobs/qjob_2017-06-08_MiSeq_S51.sh # (P170 Aonly -W) to check toxicity           # ==> 3965
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S52_R1 S52_R2 ../data/roth2016_control_set_plus_control.-100 S52 > qjobs/qjob_2017-06-08_MiSeq_S52.sh # (Roth Seaprep -W) to check complexcity          # ==> 2801344
#
#
# 2017-06-12_MiSeq (P170)
# sh Y2H_Blastn.sh 2017-06-12_MiSeq S53_R1 S53_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S53 > qjobs/qjob_2017-06-08_MiSeq_S53.sh # (P170 -W) no selection   # ==> 3250858
# sh Y2H_Blastn.sh 2017-06-12_MiSeq S54_R1 S54_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S54 > qjobs/qjob_2017-06-08_MiSeq_S54.sh # (P170 -W-A) selection    # ==> 3553164
#
### 2017-07-03_MiSeq : 60
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 60_W_R1 60_W_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 60_W # ==> 2629509
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 60_Q_R1 60_Q_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 60_Q # ==> 2507307
#
### 2017-07-03_MiSeq : 58
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 58_AW_R1 58_AW_R2 ../data/E375pLT-20170706.-100 58_AW # ==> 219
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 58_BW_R1 58_BW_R2 ../data/E375pLT-20170706.-100 58_BW # ==> 219
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 58_BQ_R1 58_BQ_R2 ../data/E375pLT-20170706.-100 58_BQ" # ==> 306
#
# 2017-08-15_MiSeq (Test for Seaprep again)
# sh Y2H_Blastn.sh 2017-08-15_MiSeq S61_SWD_R1 S61_SWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S61_SWD > qjobs/qjob_2017-08-15_MiSeq_S61_SWD.sh # SWD (P170 Seaprep -W)       # ==> 2909222 out of 6012622 (48.4%)
# sh Y2H_Blastn.sh 2017-08-15_MiSeq S61_PWD_R1 S61_PWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S61_PWD > qjobs/qjob_2017-08-15_MiSeq_S61_PWD.sh # SWD (P170 Plate -W)         # ==> 1956131 out of 4473498 (43.7%)
# sh Y2H_Blastn.sh 2017-08-15_MiSeq S61_PQD_R1 S61_PQD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S61_PQD > qjobs/qjob_2017-08-15_MiSeq_S61_PQD.sh # SWD (P170 Plate -Q)         # ==> 1684777 out of 3553019 (47.4%)
#
# 2017-08-22_MiSeq (Roth - Seaprep with selection)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SWD_S1_R1 S64-SWD_S1_R2 ../data/roth2016_control_set_plus_control.-100 S64_SWD > qjobs/qjob_2017-08-22_MiSeq_S64_SWD.sh # SWD (Roth Seaprep -W)  # ==> 1491017 out of 2562522 (58.2%)           ==> 1450453 (For New)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SA4D_S2_R1 S64-SA4D_S2_R2 ../data/roth2016_control_set_plus_control.-100 S64_SA4D > qjobs/qjob_2017-08-22_MiSeq_S64_SA4D.sh # SWD (Roth Seaprep +A1/4) # ==> 1177404 out of 2321645 (50.7%)     ==> 1145754 (For New)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SA8D_S3_R1 S64-SA8D_S3_R2 ../data/roth2016_control_set_plus_control.-100 S64_SA8D > qjobs/qjob_2017-08-22_MiSeq_S64_SA8D.sh # SWD (Roth Seaprep +A1/8) # ==> 1125083 out of 2076630 (54.2%)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SQD_S4_R1 S64-SQD_S4_R2 ../data/roth2016_control_set_plus_control.-100 S64_SQD > qjobs/qjob_2017-08-22_MiSeq_S64_SQD.sh # SWD (Roth Seaprep -Q)    # ==> 618177 out of 1437540 (43.0%)
#
# 2017-08-28_MiSeq (P170 - Seaprep)
# sh Y2H_Blastn.sh 2017-08-28_MiSeq S68_SWD_R1 S68_SWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S68_SWD > qjobs/qjob_2017-08-28_MiSeq_S68_SWD.sh # SWD (P170 Seaprep -W)            # 1999636
# sh Y2H_Blastn.sh 2017-08-28_MiSeq S68_SA4D_R1 S68_SA4D_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S68_SA4D > qjobs/qjob_2017-08-28_MiSeq_S68_SA4D.sh # SA4D (P170 Seaprep +A1/4)    # 1943064
# sh Y2H_Blastn.sh 2017-08-28_MiSeq S68_SA8D_R1 S68_SA8D_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S68_SA8D > qjobs/qjob_2017-08-28_MiSeq_S68_SA8D.sh # SA8D (P170 Seaprep +A1/8)    # 1372704
# sh Y2H_Blastn.sh 2017-08-28_MiSeq S68_SQD_R1 S68_SQD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S68_SQD > qjobs/qjob_2017-08-28_MiSeq_S68_SQD.sh # SQD (P170 Seaprep -Q)            # 1247668
#
# 2017-08-30_MiSeq (P170 - Seaprep)
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SWD_R1 S65_SWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SWD > qjobs/qjob_2017-08-30_MiSeq_S65_SWD.sh # SWD (P170 Seaprep -W)            # 1836282
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SAD_R1 S65_SAD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SAD > qjobs/qjob_2017-08-30_MiSeq_S65_SAD.sh # SAD (P170 Seaprep +A?)           # 2108754
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SQD_R1 S65_SQD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SQD > qjobs/qjob_2017-08-30_MiSeq_S65_SQD.sh # SQD (P170 Seaprep -Q)            # 1870105
#
# 2017-10-09_MiSeq (EGFR - 1)
# sh Y2H_Blastn.sh 2017-10-09_MiSeq S1_WD_R1 S1_WD_R2 ../data/EGFR_entry_clones.-100 S1_WD > qjobs/qjob_2017-10-09_MiSeq_S1_WD.sh # WD
# sh Y2H_Blastn.sh 2017-10-09_MiSeq S2_WD_R1 S2_WD_R2 ../data/EGFR_entry_clones.-100 S2_WD > qjobs/qjob_2017-10-09_MiSeq_S2_WD.sh
# sh Y2H_Blastn.sh 2017-10-09_MiSeq S3_A8D_R1 S3_A8D_R2 ../data/EGFR_entry_clones.-100 S3_A8D > qjobs/qjob_2017-10-09_MiSeq_S3_A8D.sh
# sh Y2H_Blastn.sh 2017-10-09_MiSeq S4_QD_R1 S4_QD_R2 ../data/EGFR_entry_clones.-100 S4_QD > qjobs/qjob_2017-10-09_MiSeq_S4_QD.sh
#
# sh Y2H_Blastn.sh 2017-10-16_MiSeq S1_WD_R1 S1_WD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S1_WD > qjobs/qjob_2017-10-16_MiSeq_S1_WD.sh # WD (P170 Seaprep -W) 
# sh Y2H_Blastn.sh 2017-10-16_MiSeq S2_WD_R1 S2_WD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S2_WD > qjobs/qjob_2017-10-16_MiSeq_S2_WD.sh # WD (P170 Seaprep -W) 
# sh Y2H_Blastn.sh 2017-10-16_MiSeq S3_AD_R1 S3_AD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S3_AD > qjobs/qjob_2017-10-16_MiSeq_S3_AD.sh # WD (P170 Seaprep -W) 
# sh Y2H_Blastn.sh 2017-10-16_MiSeq S4_AD_R1 S4_AD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S4_AD > qjobs/qjob_2017-10-16_MiSeq_S4_AD.sh # WD (P170 Seaprep -W) 
#
# 2017-10-30_MiSeq (R75 auto-activator)
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S1_BWD_R1 S1_BWD_R2 ../data/roth2016_control_set_plus_control.-100 S1_BWD > qjobs/qjob_2017-10-30_MiSeq_S1_BWD.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S2_BA2D_R1 S2_BA2D_R2 ../data/roth2016_control_set_plus_control.-100 S2_BA2D > qjobs/qjob_2017-10-30_MiSeq_S2_BA2D.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S3_BQD_R1 S3_BQD_R2 ../data/roth2016_control_set_plus_control.-100 S3_BQD > qjobs/qjob_2017-10-30_MiSeq_S3_BQD.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S4_AWD_R1 S4_AWD_R2 ../data/roth2016_control_set_plus_control.-100 S4_AWD > qjobs/qjob_2017-10-30_MiSeq_S4_AWD.sh
#
# 2017-11-03_MiSeq (R75 technical repeat)
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S1_W_R1 S1_W_R2 ../data/roth2016_control_set_plus_control.-100 S1_W > qjobs/qjob_2017-11-03_MiSeq_S1_W.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S2_Q_R1 S2_Q_R2 ../data/roth2016_control_set_plus_control.-100 S2_Q > qjobs/qjob_2017-11-03_MiSeq_S2_Q.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S3_W_R1 S3_W_R2 ../data/roth2016_control_set_plus_control.-100 S3_W > qjobs/qjob_2017-11-03_MiSeq_S3_W.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S4_Q_R1 S4_Q_R2 ../data/roth2016_control_set_plus_control.-100 S4_Q > qjobs/qjob_2017-11-03_MiSeq_S4_Q.sh
#
#===============================================================================


#===============================================================================
# Use all reference sequences to map
### 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 ) 
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17543_S1_R1 17543_S1_R2 ../data/roth2016_control_set_plus_control 17543_S1 Blastn_All_Ref > qjobs/all.qjob_2016-12-22_MiSeq_S1.sh # ==> 1968532 (Roth -W)
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17544_S2_R1 17544_S2_R2 ../data/roth2016_control_set_plus_control 17544_S2 Blastn_All_Ref > qjobs/all.qjob_2016-12-22_MiSeq_S2.sh # ==> 1541502 (Roth -A) from 2418986 (trimed seq) from xxx
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17545_S3_R1 17545_S3_R2 ../data/roth2016_control_set_plus_control 17545_S3 Blastn_All_Ref > qjobs/all.qjob_2016-12-22_MiSeq_S3.sh # ==> 2190445 (Roth Seaprep -W)
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17546_S4_R1 17546_S4_R2 ../data/roth2016_control_set_plus_control 17546_S4 Blastn_All_Ref > qjobs/all.qjob_2016-12-22_MiSeq_S4.sh # ==> 1877953  (Roth Seaprep -A)

### 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA ) 
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S1_R1 S1_R2 ../data/roth2016_control_set_plus_control S1 Blastn_All_Ref > qjobs/all.qjob_2017-02-22_MiSeq_S1.sh # (Roth -W) # ==> 960245
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S2_R1 S2_R2 ../data/roth2016_control_set_plus_control S2 Blastn_All_Ref > qjobs/all.qjob_2017-02-22_MiSeq_S2.sh # (Roth -A) # ==> 2024754
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S3_R1 S3_R2 ../data/roth2016_control_set_plus_control S3 Blastn_All_Ref > qjobs/all.qjob_2017-02-22_MiSeq_S3.sh # (Roth -Q) # ==> 813912

### 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S1_R1 S1_R2 ../data/roth2016_control_set_plus_control S1 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S1.sh # (Roth -W) # ==> 1414435
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S2_R1 S2_R2 ../data/roth2016_control_set_plus_control S2 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S2.sh # (Roth -A) # ==> 969159
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S3_R1 S3_R2 ../data/roth2016_control_set_plus_control S3 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S3.sh # (Roth -Q) # ==> 1255002
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S4_R1 S4_R2 ../data/roth2016_control_set_plus_control S4 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S4.sh # (Roth Aonly -W) to check toxicity # ==> 23
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S5_R1 S5_R2 ../data/roth2016_control_set_plus_control S5 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S5.sh # (Roth Bonly -W) to check toxicity # ==> 169
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S6_R1 S6_R2 ../data/roth2016_control_set_plus_control S6 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S6.sh # (Roth Bonly -A) to check auto-activation # ==> 1296
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S7_R1 S7_R2 ../data/roth2016_control_set_plus_control S7 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S7.sh # (Roth Seaprep -Q) # ==> 734766

### 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61] 
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SWD_S1_R1 S64-SWD_S1_R2 ../data/roth2016_control_set_plus_control S64_SWD Blastn_All_Ref > qjobs/all.qjob_2017-08-22_MiSeq_S64_SWD.sh # SWD (Roth Seaprep -W)  # ==> 1491017 out of 2562522 (58.2%)           ==> 1450453 (For New)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SA4D_S2_R1 S64-SA4D_S2_R2 ../data/roth2016_control_set_plus_control S64_SA4D Blastn_All_Ref > qjobs/all.qjob_2017-08-22_MiSeq_S64_SA4D.sh # SWD (Roth Seaprep +A1/4) # ==> 1177404 out of 2321645 (50.7%)     ==> 1145754 (For New)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SA8D_S3_R1 S64-SA8D_S3_R2 ../data/roth2016_control_set_plus_control S64_SA8D Blastn_All_Ref > qjobs/all.qjob_2017-08-22_MiSeq_S64_SA8D.sh # SWD (Roth Seaprep +A1/8) # ==> 1125083 out of 2076630 (54.2%)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SQD_S4_R1 S64-SQD_S4_R2 ../data/roth2016_control_set_plus_control S64_SQD Blastn_All_Ref > qjobs/all.qjob_2017-08-22_MiSeq_S64_SQD.sh # SWD (Roth Seaprep -Q)    # ==> 618177 out of 1437540 (43.0%)

### 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W; 
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S52_R1 S52_R2 ../data/roth2016_control_set_plus_control S52 Blastn_All_Ref > qjobs/all.qjob_2017-06-08_MiSeq_S52.sh # (Roth Seaprep -W) to check complexcity   

### 2017-10-30_MiSeq (R75 auto-activator)
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S1_BWD_R1 S1_BWD_R2 ../data/roth2016_control_set_plus_control S1_BWD Blastn_All_Ref > qjobs/all.qjob_2017-10-30_MiSeq_S1_BWD.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S2_BA2D_R1 S2_BA2D_R2 ../data/roth2016_control_set_plus_control S2_BA2D Blastn_All_Ref > qjobs/all.qjob_2017-10-30_MiSeq_S2_BA2D.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S3_BQD_R1 S3_BQD_R2 ../data/roth2016_control_set_plus_control S3_BQD Blastn_All_Ref > qjobs/all.qjob_2017-10-30_MiSeq_S3_BQD.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S4_AWD_R1 S4_AWD_R2 ../data/roth2016_control_set_plus_control S4_AWD Blastn_All_Ref > qjobs/all.qjob_2017-10-30_MiSeq_S4_AWD.sh

# 2017-11-03_MiSeq; Roth75 - technical repeat  
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S1_W_R1 S1_W_R2 ../data/roth2016_control_set_plus_control S1_W Blastn_All_Ref > qjobs/all.qjob_2017-11-03_MiSeq_S1_W.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S2_Q_R1 S2_Q_R2 ../data/roth2016_control_set_plus_control S2_Q Blastn_All_Ref > qjobs/all.qjob_2017-11-03_MiSeq_S2_Q.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S3_W_R1 S3_W_R2 ../data/roth2016_control_set_plus_control S3_W Blastn_All_Ref > qjobs/all.qjob_2017-11-03_MiSeq_S3_W.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S4_Q_R1 S4_Q_R2 ../data/roth2016_control_set_plus_control S4_Q Blastn_All_Ref > qjobs/all.qjob_2017-11-03_MiSeq_S4_Q.sh   


#2018-03-26_MiSeq - EGFR #2 + Sarah/Hannah micro-exons and mutant
#sh Y2H_Blastn.sh 2018-03-26_MiSeq  S1_W12D_R1 S1_W12D_R2 ../data/EGFR_entry_clones.-100 S1_W12D Blastn_All_Ref > qjobs/qjob_2018-03-26_MiSeq_S1_W12D.sh 
#sh Y2H_Blastn.sh 2018-03-26_MiSeq  S2_WHA12D_R1 S2_WHA12D_R2 ../data/EGFR_entry_clones.-100 S2_WHA12D Blastn_All_Ref > qjobs/qjob_2018-03-26_MiSeq_S2_WHA12D.sh       
#sh Y2H_Blastn.sh 2018-03-26_MiSeq  S3_W9D_R1 S3_W9D_R2 ../data/EGFR_entry_clones.-100 S3_W9D Blastn_All_Ref > qjobs/qjob_2018-03-26_MiSeq_S3_W9D.sh  
#===============================================================================


#===============================================================================
# A463 library (all)
#sh Y2H_Blastn.sh 2017-11-17_MiSeq S1_WD_R1 S1_WD_R2 ../data/A463-MGj69.RBP-MAP.-150 S1_W > qjobs/qjob_Sebastian_2017-11-17_MiSeq_S1_W2.sh
#sh Y2H_Blastn.sh 2017-11-17_MiSeq S2_QD_R1 S2_QD_R2 ../data/A463-MGj69.RBP-MAP.-150 S2_Q > qjobs/qjob_Sebastian_2017-11-17_MiSeq_S2_Q2.sh


# A463 library (w/o auto activator)
#sh Y2H_Blastn.sh 2017-12-07_MiSeq S1_WD_R1 S1_WD_R2 ../data/A463-MGj69.RBP-MAP.-150 S1_W > qjobs/qjob_Sebastian_2017-12-07_MiSeq_S1_W.sh                
#sh Y2H_Blastn.sh 2017-12-07_MiSeq S2_QD_R1 S2_QD_R2 ../data/A463-MGj69.RBP-MAP.-150 S2_Q > qjobs/qjob_Sebastian_2017-12-07_MiSeq_S2_Q.sh                

# A463 library (w/o auto activator)
#sh Y2H_Blastn.sh 2018-02-21_MiSeq S1_WD_R1 S1_WD_R2 ../data/A463-MGj69.RBP-MAP.-150 S1_W > qjobs/qjob_Sebastian_2018-02-21_MiSeq_S1_W.sh                
#sh Y2H_Blastn.sh 2018-02-21_MiSeq S2_QD_R1 S2_QD_R2 ../data/A463-MGj69.RBP-MAP.-150 S2_Q > qjobs/qjob_Sebastian_2018-02-21_MiSeq_S2_Q.sh                
#===============================================================================


# make output folder
echo "cd /share"
cd /share

if [ ! -d output/$1 ]; then
   echo "mkdir output/$1"
   mkdir output/$1
fi

if [ ! -d output/$1/$6 ]; then
   echo "mkdir output/$1/$6"
   mkdir output/$1/$6
fi



# remove 5' end adaptor sequences in both fastq file; remove reads if any of read does not have an adaptor sequence ( Temporarily now using fastq files in Friedrich folder which are already trimmed )
# -m 15 : Discard trimmed reads that are shorter than LENGTH.
# --discard-untrimmed : Discard reads that do not contain the adapter.
#
#	bait      --CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT
#	prey      GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT
#	             ***  ** **** ****************************************
#
echo "cutadapt -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o output/$1/$6/$2.fastq -p output/$1/$6/$3.fastq -m 15 --discard-untrimmed ../fastq/$2.fastq.gz ../fastq/$3.fastq.gz"
cutadapt -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o output/$1/$6/$2.fastq -p output/$1/$6/$3.fastq -m 15 --discard-untrimmed ./fastq/$2.fastq.gz ./fastq/$3.fastq.gz
#cutadapt -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o output/$1/Blastn/$2.fastq -p output/$1/Blastn/$3.fastq -m 15 --discard-untrimmed ../$1/$2.fastq.gz ../$1/$3.fastq.gz

# convert fastq to fasta ( Temporarily now using fastq file generated in Friedrich folder; since it is the same fastq. But We need to change it to Blastn for general cases )
#echo "python main.py fastq_to_fasta output/$1/Friedrich/$2.fastq > output/$1/Blastn/$2.fa"
#echo "python main.py fastq_to_fasta output/$1/Friedrich/$3.fastq > output/$1/Blastn/$3.fa"
echo "main.py fastq_to_fasta output/$1/$6/$2.fastq > output/$1/$6/$2.fa"
main.py fastq_to_fasta output/$1/$6/$2.fastq > output/$1/$6/$2.fa
echo "main.py fastq_to_fasta output/$1/$6/$3.fastq > output/$1/$6/$3.fa"
main.py fastq_to_fasta output/$1/$6/$3.fastq > output/$1/$6/$3.fa
#python main.py fastq_to_fasta output/2017-07-03_MiSeq/Friedrich/60_Q_R1.fastq > output/2017-07-03_MiSeq/Blastn/60_Q_R1.fa
#python main.py fastq_to_fasta output/2017-07-03_MiSeq/Friedrich/60_Q_R2.fastq > output/2017-07-03_MiSeq/Blastn/60_Q_R2.fa


# blastn-short search
echo "blastn -db $4.fa  -query output/$1/$6/$2.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/$1/$6/$2.blastn"
blastn -db $4.fa  -query output/$1/$6/$2.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/$1/$6/$2.blastn
#blastn -db ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa  -query output/2017-07-03_MiSeq/Blastn/60_Q_R1.fa -task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 1e-8 > output/2017-07-03_MiSeq/Blastn/60_Q_R1.blastn
#echo "blastn -db ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa  -query output/2017-07-03_MiSeq/Blastn/60_Q_R2.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn"
echo "blastn -db $4.fa  -query output/$1/$6/$3.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/$1/$6/$3.blastn"
blastn -db $4.fa  -query output/$1/$6/$3.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/$1/$6/$3.blastn
#blastn -db ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa  -query output/2017-07-03_MiSeq/Blastn/60_Q_R2.fa -task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 1e-8 > output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn
#echo "blastn -db ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa  -query output/$1/Friedrich/$2.fastq -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8"


# python parse blastn output and make ppi map
echo "main.py BLASTN $4.fa output/$1/$6/$2.blastn output/$1/$6/$3.blastn > output/$1/$6/$5.ppi.txt"
main.py BLASTN $4.fa output/$1/$6/$2.blastn output/$1/$6/$3.blastn > output/$1/$6/$5.ppi.txt
# cat output/2017-07-03_MiSeq/Blastn/60_Q_R1.blastn | awk '{print $1}' | uniq | wc    ## ==> 3209686
# cat output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn | awk '{print $1}' | uniq | wc    ## ==> 3029067
# python main.py BLASTN ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa output/2017-07-03_MiSeq/Blastn/60_Q_R1.blastn output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn    ##  ==> 2507307



# echo "python main.py BLASTN_BARCODE $4.fa output/$1/$6/$2.blastn output/$1/$6/$3.blastn > output/$1/$6/$5.ppi.txt" # 2018/03/28 







