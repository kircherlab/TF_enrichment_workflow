# TF_enrichment_workflow

This repository contains a snakemake workflow to calculate enrichment of transcription factors in cancer open chromatin regions compared to healthy open chromatin regions.

## Authors

- Lea Burkard (@leaburkard)

## Usage

### Open chromatin regions

This workflow was used with the following data:

- blood open chromatin regions: Corces et al (2016) Nat Genet. https://doi.org/10.1038/ng.3646
- cancer open chromatin regions: Corces et al (2018) Science https://doi.org/10.1126/science.aav1898

The data can be downloaded and preprocessed by exectuing the `download_data.sh` script.


### Transcription factors

To use this workflow transcription factor data needs to be downloaded and added to the `transcriptionfactors/` folder. This worklow was designed and only used with the predicted TFBS of the JASPAR Core Collection. An instruction how to download them is written in the directory. 


### Input parameters

In the beginning of the workflow (`snakefile_enriched_tf.smk`) the file names for the healthy and cancer open chromatin regions need to be given without extensions. These files need to be in the `input/` folder.

```bash
#### Name of healthy and cancer file ####
HEALTHY = "blood_mono_500"
CANCER = "pancancer_COAD_500"
#########################################
```
