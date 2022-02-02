import pandas as pd
import numpy as np

pancancer = pd.read_table('input/TCGA-ATAC_PanCancer_PeakSet.txt')

counts = pancancer['score'].astype(bool).sum(axis=0) #count non-zero rows


#normalize the data
df_cancer_normalized = pd.DataFrame(data=pancancer)

totalSum = pancancer['score'].sum()
df_cancer_normalized['score'] = pancancer['score'] / totalSum * 1000000

#Width
peakWidth = np.array(df_cancer_normalized["end"])-np.array(df_cancer_normalized["start"])
print("peak width: max ",peakWidth.max())
print("peak width: min ",peakWidth.min())


df_cancer_normalized["end"] = df_cancer_normalized["end"]-1

cancer_types = []
for cancer in df_cancer_normalized['name'].array: #.split("_")
    temp = cancer.split("_")
    if temp[0] not in cancer_types:
        cancer_types.append(temp[0])
        


#save each cancer type 
for cancer in cancer_types:
    l = len(cancer)
    #(500 bp width)
    cancer_type_normalized = df_cancer_normalized[df_cancer_normalized['name'].str[0:l] == cancer]
    cancer_500_bed = cancer_type_normalized.loc[:,["seqnames","start","end","score"]]
    cancer_500_bed.to_csv("input/pancancer_"+cancer+"_500.bed.gz", index=False, header=False, compression="gzip", sep ='\t')
    print("input/pancancer_"+cancer+"_500.bed.gz file created")

    #(300 bp width)
    normalized_300 = pd.DataFrame(data=cancer_type_normalized)
    normalized_300["start"] = normalized_300["start"]+100
    normalized_300["end"] = normalized_300["end"]-100
    
    cancer_300_bed = normalized_300.loc[:,["seqnames","start","end","score"]]
    cancer_300_bed.to_csv("input/pancancer_"+cancer+"_300.bed.gz", index=False, header=False, compression="gzip", sep ='\t')
    print("input/pancancer_"+cancer+"_300.bed.gz file created")
    
