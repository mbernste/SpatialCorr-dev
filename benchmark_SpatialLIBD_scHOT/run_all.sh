########################################################################################################################################
#   All commands needed to run the scHOT/SpatialCorr benchmark. Note, this should be refactored into a Snakemake pipeline.
########################################################################################################################################

# 1
python run_SpatialDC_SpatialLIBD_one_pair.py PLD3,MYL12B -o results_timing_rand_pairs/SpatialDC/PLD3_MYL12B.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py PLD3,MYL12B -n -o results_timing_rand_pairs/SpatialDC_no_mc/PLD3_MYL12B.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py PLD3,MYL12B -b -o results_timing_rand_pairs/SpatialCorr_br/PLD3_MYL12B.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py PLD3,MYL12B -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/PLD3_MYL12B.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py PLD3,MYL12B -o results_timing_rand_pairs/scHOT/PLD3_MYL12B.scHOT.tsv

# 2
python run_SpatialDC_SpatialLIBD_one_pair.py MDH1,JTB -o results_timing_rand_pairs/SpatialDC/MDH1_JTB.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MDH1,JTB -n -o results_timing_rand_pairs/SpatialDC_no_mc/MDH1_JTB.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MDH1,JTB -b -o results_timing_rand_pairs/SpatialCorr_br/MDH1_JTB.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MDH1,JTB -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/MDH1_JTB.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py MDH1,JTB -o results_timing_rand_pairs/scHOT/MDH1_JTB.scHOT.tsv

# 3
python run_SpatialDC_SpatialLIBD_one_pair.py ANAPC11,ECHS1 -o results_timing_rand_pairs/SpatialDC/ANAPC11_ECHS1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ANAPC11,ECHS1 -n -o results_timing_rand_pairs/SpatialDC_no_mc/ANAPC11_ECHS1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ANAPC11,ECHS1 -b -o results_timing_rand_pairs/SpatialCorr_br/ANAPC11_ECHS1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ANAPC11,ECHS1 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/ANAPC11_ECHS1.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py ANAPC11,ECHS1 -o results_timing_rand_pairs/scHOT/ANAPC11_ECHS1.scHOT.tsv

# 4
python run_SpatialDC_SpatialLIBD_one_pair.py HNRNPA3,ACOT7 -o results_timing_rand_pairs/SpatialDC/HNRNPA3_ACOT7.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py HNRNPA3,ACOT7 -n -o results_timing_rand_pairs/SpatialDC_no_mc/HNRNPA3_ACOT7.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py HNRNPA3,ACOT7 -b -o results_timing_rand_pairs/SpatialCorr_br/HNRNPA3_ACOT7.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py HNRNPA3,ACOT7 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/HNRNPA3_ACOT7.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py HNRNPA3,ACOT7 -o results_timing_rand_pairs/scHOT/HNRNPA3_ACOT7.scHOT.tsv

# 5
python run_SpatialDC_SpatialLIBD_one_pair.py LZTS3,FAU -o results_timing_rand_pairs/SpatialDC/LZTS3_FAU.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py LZTS3,FAU -n -o results_timing_rand_pairs/SpatialDC_no_mc/LZTS3_FAU.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py LZTS3,FAU -b -o results_timing_rand_pairs/SpatialCorr_br/LZTS3_FAU.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py LZTS3,FAU -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/LZTS3_FAU.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py LZTS3,FAU -o results_timing_rand_pairs/scHOT/LZTS3_FAU.scHOT.tsv

# 6
python run_SpatialDC_SpatialLIBD_one_pair.py ALDOC,CHMP3 -o results_timing_rand_pairs/SpatialDC/ALDOC_CHMP3.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ALDOC,CHMP3 -n -o results_timing_rand_pairs/SpatialDC_no_mc/ALDOC_CHMP3.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ALDOC,CHMP3 -b -o results_timing_rand_pairs/SpatialCorr_br/ALDOC_CHMP3.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ALDOC,CHMP3 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/ALDOC_CHMP3.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py ALDOC,CHMP3 -o results_timing_rand_pairs/scHOT/ALDOC_CHMP3.scHOT.tsv

# 7
python run_SpatialDC_SpatialLIBD_one_pair.py SRSF5,CENPX -o results_timing_rand_pairs/SpatialDC/SRSF5_CENPX.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py SRSF5,CENPX -n -o results_timing_rand_pairs/SpatialDC_no_mc/SRSF5_CENPX.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py SRSF5,CENPX -b -o results_timing_rand_pairs/SpatialCorr_br/SRSF5_CENPX.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py SRSF5,CENPX -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/SRSF5_CENPX.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py SRSF5,CENPX -o results_timing_rand_pairs/scHOT/SRSF5_CENPX.scHOT.tsv

# 8
python run_SpatialDC_SpatialLIBD_one_pair.py RPS17,NFE2L1 -o results_timing_rand_pairs/SpatialDC/RPS17_NFE2L1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py RPS17,NFE2L1 -n -o results_timing_rand_pairs/SpatialDC_no_mc/RPS17_NFE2L1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py RPS17,NFE2L1 -b -o results_timing_rand_pairs/SpatialCorr_br/RPS17_NFE2L1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py RPS17,NFE2L1 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/RPS17_NFE2L1.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py RPS17,NFE2L1 -o results_timing_rand_pairs/scHOT/RPS17_NFE2L1.scHOT.tsv

# 9
python run_SpatialDC_SpatialLIBD_one_pair.py MAP2K1,CAMTA2 -o results_timing_rand_pairs/SpatialDC/MAP2K1_CAMTA2.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MAP2K1,CAMTA2 -n -o results_timing_rand_pairs/SpatialDC_no_mc/MAP2K1_CAMTA2.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MAP2K1,CAMTA2 -b -o results_timing_rand_pairs/SpatialCorr_br/MAP2K1_CAMTA2.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MAP2K1,CAMTA2 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/MAP2K1_CAMTA2.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py MAP2K1,CAMTA2 -o results_timing_rand_pairs/scHOT/MAP2K1_CAMTA2.scHOT.tsv

# 10
python run_SpatialDC_SpatialLIBD_one_pair.py ATP1A3,FAM162A -o results_timing_rand_pairs/SpatialDC/ATP1A3_FAM162A.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ATP1A3,FAM162A -n -o results_timing_rand_pairs/SpatialDC_no_mc/ATP1A3_FAM162A.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ATP1A3,FAM162A -b -o results_timing_rand_pairs/SpatialCorr_br/ATP1A3_FAM162A.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py ATP1A3,FAM162A -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/ATP1A3_FAM162A.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py ATP1A3,FAM162A -o results_timing_rand_pairs/scHOT/ATP1A3_FAM162A.scHOT.tsv

# 11
python run_SpatialDC_SpatialLIBD_one_pair.py SARS,QDPR -o results_timing_rand_pairs/SpatialDC/SARS_QDPR.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py SARS,QDPR -n -o results_timing_rand_pairs/SpatialDC_no_mc/SARS_QDPR.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py SARS,QDPR -b -o results_timing_rand_pairs/SpatialCorr_br/SARS_QDPR.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py SARS,QDPR -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/SARS_QDPR.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py SARS,QDPR -o results_timing_rand_pairs/scHOT/SARS_QDPR.scHOT.tsv

# 12
python run_SpatialDC_SpatialLIBD_one_pair.py TMEM130,DAD1 -o results_timing_rand_pairs/SpatialDC/TMEM130_DAD1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py TMEM130,DAD1 -n -o results_timing_rand_pairs/SpatialDC_no_mc/TMEM130_DAD1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py TMEM130,DAD1 -b -o results_timing_rand_pairs/SpatialCorr_br/TMEM130_DAD1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py TMEM130,DAD1 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/TMEM130_DAD1.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py TMEM130,DAD1 -o results_timing_rand_pairs/scHOT/TMEM130_DAD1.scHOT.tsv

# 13
python run_SpatialDC_SpatialLIBD_one_pair.py NOVA1,HNRNPH1 -o results_timing_rand_pairs/SpatialDC/NOVA1_HNRNPH1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py NOVA1,HNRNPH1 -n -o results_timing_rand_pairs/SpatialDC_no_mc/NOVA1_HNRNPH1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py NOVA1,HNRNPH1 -b -o results_timing_rand_pairs/SpatialCorr_br/NOVA1_HNRNPH1.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py NOVA1,HNRNPH1 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/NOVA1_HNRNPH1.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py NOVA1,HNRNPH1 -o results_timing_rand_pairs/scHOT/NOVA1_HNRNPH1.scHOT.tsv

# 14
python run_SpatialDC_SpatialLIBD_one_pair.py CCDC85B,UQCRH -o results_timing_rand_pairs/SpatialDC/CCDC85B_UQCRH.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py CCDC85B,UQCRH -n -o results_timing_rand_pairs/SpatialDC_no_mc/CCDC85B_UQCRH.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py CCDC85B,UQCRH -b -o results_timing_rand_pairs/SpatialCorr_br/CCDC85B_UQCRH.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py CCDC85B,UQCRH -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/CCDC85B_UQCRH.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py CCDC85B,UQCRH -o results_timing_rand_pairs/scHOT/CCDC85B_UQCRH.scHOT.tsv

# 15
python run_SpatialDC_SpatialLIBD_one_pair.py PDIA6,NUDT14 -o results_timing_rand_pairs/SpatialDC/PDIA6_NUDT14.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py PDIA6,NUDT14 -n -o results_timing_rand_pairs/SpatialDC_no_mc/PDIA6_NUDT14.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py PDIA6,NUDT14 -b -o results_timing_rand_pairs/SpatialCorr_br/PDIA6_NUDT14.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py PDIA6,NUDT14 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/PDIA6_NUDT14.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py PDIA6,NUDT14 -o results_timing_rand_pairs/scHOT/PDIA6_NUDT14.scHOT.tsv

# 16
python run_SpatialDC_SpatialLIBD_one_pair.py RTN3,RPL18 -o results_timing_rand_pairs/SpatialDC/RTN3_RPL18.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py RTN3,RPL18 -n -o results_timing_rand_pairs/SpatialDC_no_mc/RTN3_RPL18.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py RTN3,RPL18 -b -o results_timing_rand_pairs/SpatialCorr_br/RTN3_RPL18.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py RTN3,RPL18 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/RTN3_RPL18.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py RTN3,RPL18 -o results_timing_rand_pairs/scHOT/RTN3_RPL18.scHOT.tsv

# 17
python run_SpatialDC_SpatialLIBD_one_pair.py MTURN,TAF10 -o results_timing_rand_pairs/SpatialDC/MTURN_TAF10.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MTURN,TAF10 -n -o results_timing_rand_pairs/SpatialDC_no_mc/MTURN_TAF10.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MTURN,TAF10 -b -o results_timing_rand_pairs/SpatialCorr_br/MTURN_TAF10.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MTURN,TAF10 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/MTURN_TAF10.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py MTURN,TAF10 -o results_timing_rand_pairs/scHOT/MTURN_TAF10.scHOT.tsv

# 18
python run_SpatialDC_SpatialLIBD_one_pair.py CBX6,FAM19A5 -o results_timing_rand_pairs/SpatialDC/CBX6_FAM19A5.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py CBX6,FAM19A5 -n -o results_timing_rand_pairs/SpatialDC_no_mc/CBX6_FAM19A5.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py CBX6,FAM19A5 -b -o results_timing_rand_pairs/SpatialCorr_br/CBX6_FAM19A5.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py CBX6,FAM19A5 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/CBX6_FAM19A5.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py CBX6,FAM19A5 -o results_timing_rand_pairs/scHOT/CBX6_FAM19A5.scHOT.tsv

# 19
python run_SpatialDC_SpatialLIBD_one_pair.py EPS15,SEPT2 -o results_timing_rand_pairs/SpatialDC/EPS15_SEPT2.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py EPS15,SEPT2 -n -o results_timing_rand_pairs/SpatialDC_no_mc/EPS15_SEPT2.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py EPS15,SEPT2 -b -o results_timing_rand_pairs/SpatialCorr_br/EPS15_SEPT2.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py EPS15,SEPT2 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/EPS15_SEPT2.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py EPS15,SEPT2 -o results_timing_rand_pairs/scHOT/EPS15_SEPT2.scHOT.tsv

# 20
python run_SpatialDC_SpatialLIBD_one_pair.py MRPL33,RPL24 -o results_timing_rand_pairs/SpatialDC/MRPL33_RPL24.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MRPL33,RPL24 -n -o results_timing_rand_pairs/SpatialDC_no_mc/MRPL33_RPL24.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MRPL33,RPL24 -b -o results_timing_rand_pairs/SpatialCorr_br/MRPL33_RPL24.SpatialDC.tsv
python run_SpatialDC_SpatialLIBD_one_pair.py MRPL33,RPL24 -n -b -o results_timing_rand_pairs/SpatialCorr_br_no_mc/MRPL33_RPL24.SpatialDC.tsv
python run_scHOT_SpatialLIBD_one_pair.py MRPL33,RPL24 -o results_timing_rand_pairs/scHOT/MRPL33_RPL24.scHOT.tsv



