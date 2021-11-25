import glob, os.path
import pandas as pd


#### Name of healthy and cancer file ####
HEALTHY = "blood_mono_500"
CANCER = "pancancer_COAD_500"
#########################################

FILES = glob.glob("transcriptionfactors/*.tsv.gz")
TF = [os.path.basename(x)[:-7] for x in FILES]

SAMPLES = [HEALTHY,CANCER]

rule all:
	input:
		expand("intersect_{samples}/{tf}.intersect_single.bed.gz", samples=SAMPLES, tf=TF),
		expand("results/{cancer}_enriched_intersect.txt", cancer=CANCER)


rule tf_annotation:
	input:
		expand("transcriptionfactors/{tf}.tsv.gz", tf=TF)
	output:
		"resources/tf_annotation.txt"
	shell:
		"""
		printf "filename\twidth\tTF\tlines\n" > {output}

		for f in transcriptionfactors/*.tsv.gz
		do
			file=$f
			filename=$(basename $file)
				width=$(zcat $file | head -1 | awk '{{print $3-$2}}')
				tf=$(zcat $file | head -1 | awk '{{print $4}}')
				lines=$(zcat $file | wc -l)
				printf "$filename\t$width\t$tf\t$lines\n" >> {output}
		done
		"""


rule sort_tf:
	input:
		"transcriptionfactors/{tf}.tsv.gz"
	output:
		"transcriptionfactors/sorted/{tf}.tsv.gz"
	shell:
		"zcat {input} | sort -k1,1 -k2,2n | gzip -c > {output}"


rule bedtools_subtract:
	input:
		file="input/{sample}.bed.gz",
		blacklist="resources/hg38-blacklist.v2.bed.gz"
	output:
		"subtract/{sample}.subtract.bed.gz"
	shell:
		"bedtools subtract -A -a {input.file} -b {input.blacklist} | sort -k1,1 -k2,2n | bedtools merge | gzip -c  > {output}"


rule get_score:
	input:
		file="input/{sample}.bed.gz",
		subtract="subtract/{sample}.subtract.bed.gz"
	output:
		"subtract/{sample}_score.subtract.bed.gz"
	shell:
		"bedtools intersect -a {input.file} -b {input.subtract} | sort -k1,1 -k2,2n | gzip -c > {output}"


rule intersect_blood:
	input:
		file="subtract/"+HEALTHY+"_score.subtract.bed.gz",
		factor="transcriptionfactors/sorted/{tf}.tsv.gz"
	output:
		"intersect_"+HEALTHY+"/{tf}.intersect.bed.gz"
	shell:
		"bedtools intersect -a {input.factor} -b {input.file} -f 1 | cut -f1-5,7 | awk '{{print $0 , \"\thg38:\" $1 \":\" $2 \"-\" $3}}' | gzip -c > {output}"


rule intersect_cancer:
	input:
		file="subtract/"+CANCER+"_score.subtract.bed.gz",
		factor="transcriptionfactors/sorted/{tf}.tsv.gz"
	output:
		TF_to_P=temp("intersect_"+CANCER+"/{tf}_1.intersect.bed.gz"),
		P_to_TF=temp("intersect_"+CANCER+"/{tf}_2.intersect.bed.gz"),
		intersect="intersect_"+CANCER+"/{tf}.intersect.bed.gz"
	shell:
		"""
		bedtools intersect -a {input.factor} -b {input.file} -f 1 | cut -f1-5,7 | awk '{{print $0 , \"\thg38:\" $1 \":\" $2 \"-\" $3}}' | gzip -c > {output.TF_to_P}
		bedtools intersect -a {input.file} -b {output.TF_to_P} -F 1 | gzip -c > {output.P_to_TF}
		scripts/cancer_intersect.py -t {output.TF_to_P} -p {output.P_to_TF} -o {output.intersect}
		"""


rule enriched:
	input:
		expand("intersect_{samples}/{tf}.intersect.bed.gz", samples=SAMPLES, tf=TF),
		tf_info="resources/tf_annotation.txt"
	output:
		expand("results/{cancer}_enriched_intersect.txt", cancer=CANCER)
	shell:
		"""
		blood_total=$(zcat input/{HEALTHY}.bed.gz | wc -l)
		cancer_total=$(zcat input/{CANCER}.bed.gz | wc -l)
		
		printf "MA\tname\tblood_lines\tcancer_lines\ttf_lines\tblood_ratio\tcancer_ratio\ttotal_ratio\n" > {output}

		for file in transcriptionfactors/sorted/*.tsv.gz
		do
			filename=$(basename $file)
			echo $filename
			tf_name=$(cat {input.tf_info} | grep $filename | cut -f 3)
			tf_lines=$(zcat $file | wc -l)
			tf=$(echo $filename | cut -d'.' -f 1,2)

			blood_lines=$(zcat intersect_{HEALTHY}/$tf.intersect.bed.gz | wc -l)
			cancer_lines=$(zcat intersect_{CANCER}/$tf.intersect.bed.gz | wc -l)
			blood_ratio=$(gawk "BEGIN {{print $blood_lines / $blood_total;}}")
			cancer_ratio=$(gawk "BEGIN {{print $cancer_lines / $cancer_total;}}") 
			total_ratio=$(gawk "BEGIN {{print $cancer_ratio / $blood_ratio;}}")
			
			printf "$filename\t$tf_name\t$blood_lines\t$cancer_lines\t$tf_lines\t$blood_ratio\t$cancer_ratio\t$total_ratio\n" >> {output}
		done
		"""


rule single_peak:
	input:
		a="input/{samples}.bed.gz",
		b="intersect_{samples}/{tf}.intersect.bed.gz",
		file="intersect_{samples}/{tf}.intersect.bed.gz"
	output:
		wa=temp("intersect_{samples}/{tf}.wa.bed.gz"),
		single="intersect_{samples}/{tf}.intersect_single.bed.gz"
	shell:
		"""
		bedtools intersect -a {input.a} -b {input.b} -wa | sort -k1,1 -k2,2n | gzip -c > {output.wa}
		paste <(zcat {output.wa} | cut -f1-3) <(zcat {input.file}) | sort -u -k1,1 -k2,2n | cut -f4- | sort -k1,1 -k2,2n | gzip -c > {output.single}
		"""