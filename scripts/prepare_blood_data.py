import pandas as pd
import numpy as np

# read in the data
bloodData = pd.read_table('resources/GSE74912_ATACseq_All_Counts.txt.gz')

# get dataframe with only monocytes
mono_data = [bloodData["Chr"],bloodData["Start"],bloodData["End"],bloodData["4983-7A"],bloodData["4983-7B"],bloodData["6792-7A"],bloodData["6792-7B"],bloodData["Donor7256-7A-141106"],bloodData["Donor7256-7B-141106"]]
mono_header = ["Chr","Start","End","4983_7A","4983_7B","6792_7A","6792_7B","Donor7256_7A_141106","Donor7256_7B_141106"]
mono_df = pd.concat(mono_data, axis=1, keys=mono_header)

# remove rows with all zeros
mono_df = mono_df.loc[(mono_df.loc[:,"4983_7A":"Donor7256_7B_141106"] != 0).any(axis=1)] 


# majority vote
majority_mono = mono_df[mono_df.loc[:,"4983_7A":"Donor7256_7B_141106"].astype(bool).sum(axis=1) >= 4]


# normalize the majority vote
normalized_mono = pd.DataFrame(data=majority_mono)
monocytes = ["4983_7A","4983_7B","6792_7A","6792_7B","Donor7256_7A_141106","Donor7256_7B_141106"]
for col in monocytes:
    totalSum = mono_df[col].sum()
    normalized_mono[col] = mono_df[col] / totalSum * 1000000



# mean / median
final_df_500 = pd.DataFrame(data=normalized_mono.loc[:,"Chr":"End"])

mean_counts = normalized_mono.loc[:,"4983_7A":"Donor7256_7B_141106"].mean(axis=1)
median_counts = normalized_mono.loc[:,"4983_7A":"Donor7256_7B_141106"].median(axis=1)
diff = mean_counts - median_counts

final_df_500["mean_counts"] = mean_counts
#final_df_500["median_counts"] = median_counts


# Width
peakWidth = np.array(bloodData["End"])-np.array(bloodData["Start"])
print("peak width: max ",peakWidth.max())
print("peak width: min ",peakWidth.min())


# save new data (pandas df as gzipped file)
final_df_500.to_csv("input/blood_mono_500.bed.gz", index=False, header=False, compression="gzip", sep ='\t')
print("input/blood_mono_500.bed.gz file created")

# make dataframe with only 300 bp peak width
final_df_300 = pd.DataFrame(data=final_df_500)
final_df_300["Start"] = final_df_300["Start"]+100
final_df_300["End"] = final_df_300["End"]-100

final_df_300.to_csv("input/blood_mono_300.bed.gz", index=False, header=False, compression="gzip", sep ='\t')
print("input/blood_mono_300.bed.gz file created")