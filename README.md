# BarleyPangenomeVariantExtract

Repository of scripts used as part of the HMA (heavy metal ATPase) gene family publication.

The scripts are wrapped into a Snakemake pipeline ([Koster, 2019](https://doi.org/10.1093/bioinformatics/bts480); [Molder, 2021](https://doi.org/10.12688/f1000research.29032.2)). It combines a blast search of a candidate gene against the barley pangenome V2 sequences, sequence extraction, mapping against a reference sequence and variant calling.

## Requirements

1. Installation of Snakemake (tested with Snakemake v.9.5.1.)
2. Download of the genome sequences as published in Jayakodi et al. ([Jayakodi, 2024](https://doi.org/10.1038/s41586-024-08187-1)) and add to a folder called Input/assemblies
3. Creation of the corresponding blast databases from the fasta sequences using ``` makeblastdb ```




## Running of the scripts

```
snakemake -s Snakefile.smk --cores 4 --use-conda
```


## References
Jayakodi, M., Lu, Q., Pidon, H., Rabanus-Wallace, M.T., Bayer, M., Lux, T., et al. (2024). Structural variation in the pangenome of wild and domesticated barley. Nature 636(8043), 654-662. doi: 10.1038/s41586-024-08187-1

Koster, J., and Rahmann, S. (2012). Snakemake--a scalable bioinformatics workflow engine. Bioinformatics 28(19), 2520-2522. doi: 10.1093/bioinformatics/bts480

Molder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., et al. (2021). Sustainable data analysis with Snakemake. F1000Res 10, 33. doi: 10.12688/f1000research.29032.2