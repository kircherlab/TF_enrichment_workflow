#!/bin/bash

# blood data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74912/soft/GSE74912_family.soft.gz -O resources/GSE74912_family.soft.gz
python3 scriptsprepare_blood_data.py 

# cancer data
wget https://api.gdc.cancer.gov/data/116ebba2-d284-485b-9121-faf73ce0a4ec -O resources/TCGA-ATAC_PanCancer_PeakSet.txt
python3 scripts/prepare_cancer_data.py 